# --------------------------
# Step 7: Aim 3 data-driven intra-individual analysis
# --------------------------

# ---------------------
# source functions & load libraries  
# ---------------------
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path ))
source("Functions.R")
load_packs()

main_output_dir <- "Files/Output/Individual_graphs" 

# ---------------------
# read in & clean data  
# ---------------------
# WP plots for 50% inclusion data
path <- "Files/20221013_50pct_inclusion/to_model_dsm_row2.rds"
dsm_speed <- readRDS(path)
model_df <- dsm_speed$data
poly_coefs <- poly(model_df$MinLag_0_WP, 2, raw = FALSE)[,2]
model_df <- model_df %>% mutate(poly_coefs = poly_coefs)
mydf_clean <- read.csv("Files/20221013_50pct_inclusion/mydf_clean_20221014.csv")

test <- "dsm"
L <-
  # clean EMA data - select relevant survey variables and variables to join on
  mydf_clean %>% 
  filter(test == test) %>%
  rename(survey_interruptions = interruptions_anyInterruptions) %>%
  select(user_id, contains("survey"), MinLag_0, value, 
         -survey_glucose_estimate_numeric) %>%
  {
    bind_cols(
      select_if(., is.numeric),
      select_at(., "user_id")) 
  } %>%
  # compute WP deviations
  group_by(user_id) %>%
  mutate_at(vars(contains("survey"), MinLag_0, value),
            funs(WP = (. - mean(., na.rm = TRUE)) / sd(., na.rm = TRUE))) %>%
  ungroup() %>%
  inner_join(model_df %>% unique(.)) %>%
  select(contains("survey"), contains("poly_coefs"), user_id, 
         contains("MinLag_0"), contains("nightly_gl"), contains("value")) %>%
  unique(.) %>%
  group_by(user_id) %>%
  nest() %>%
  ungroup() %>%
  mutate(data = map(.x = data, 
                    ~.x %>% 
                      select_if(~(sum(!is.na(.)) > .5*sum(is.na(.) | !is.na(.)))) %>%
                      select_if(funs(var(., na.rm = TRUE) != 0)))) %>%
  mutate(nrow = map(.x = data, ~nrow(.x))) %>%
  filter(nrow > 10) %>% select(-nrow)

# ---------------------
# combine model + EMA data
# ---------------------
dsm_speed_list <- list(a = dsm_speed)
names(dsm_speed_list) <- path
wrangle_cat_new <- wrangle_func(dsm_speed_list, "50pct")

cat_coefs <- 
  wrangle_cat_new$data[[1]] %>%
  select(user_id = term, condition_mean_quad) %>% 
  unique(.) %>%
  right_join(L) 

graph_distributions <-
  cat_coefs %>%
  unnest(data) %>%
  rename_at(.vars = vars(ends_with("_score")),
            .funs = funs(sub("_score", "", .))) %>%
  rename_at(.vars = vars(ends_with("_emotion")),
            .funs = funs(sub("_emotion", "", .))) %>%
  select(user_id, (!contains("_WP") & contains("survey"))) %>%
  gather(key = key, value = value, -user_id) %>%
  mutate(key = str_remove_all(key, "survey_")) %>%
  ggplot(., aes(x = value)) +
  geom_histogram(alpha = .5, bins = 15, aes(fill = "50pct")) +
  facet_wrap(~ key, scales = "free", ncol = 4) +
  theme_bw() +
  labs(y = "Count", x = "", fill = "EMA\ncompletion")

ggsave("Files/Output/Aim3_distributions.tiff", 
       graph_distributions, 
       units = "in", width = 12, height = 5)

# ---------------------
# plot individual graphs
# ---------------------
# forcing network model to retain MinLag_0 x value path
# other paths retained if in top 10% (excl. diagonal)
legend_labels <-
  cat_coefs %>% unnest(c(data)) %>% 
  select((contains("_WP") & contains("survey")), 
         value_WP, MinLag_0_WP) %>% colnames(.)
color_key <- data.frame(var = legend_labels, 
                        col = brewer.pal(legend_labels %>% length(), "Set3"))

count <- 0
for (i in 1:nrow(cat_coefs)) {
  tryCatch(
    expr = {
      paths_tmp <-
        plot_graph(df = cat_coefs$data[[i]], 
                   decile = 10, 
                   color_key = color_key, 
                   out_dir = paste0(main_output_dir, "/idio_"), 
                   with_names = FALSE,
                   user_id = cat_coefs$user_id[i])
      count <- count + 1
    },
    error = function(e) {
      paths_tmp <- NA_real_
    }
  )
  
  if (i == 1 & is.data.frame(paths_tmp)) {
    paths <- paths_tmp
  } else {
    paths <- bind_rows(paths, paths_tmp)
  }
  print(count)
}

full_color_key <- color_key %>%
  mutate(var = str_remove(var, "survey_"),
         var = str_remove(var, "_WP"),
         var = str_remove(var, "_score"),
         var = str_remove(var, "_numeric"),
         var = if_else(var == "value", "DSM_medRT", var),
         var = if_else(str_detect(var, "Lag"), "glucose", var),
         var = if_else(str_detect(var, "DSM"), "dsm_RT", var))

graphics.off()
png(filename = paste0(main_output_dir, "/legend.png"), 
    width = 4, height = 5, units = "in", res = 300)
plot.new()
legend("center", 
       legend = full_color_key$var %>% str_remove_all(., "_emotion"), 
       pt.bg= full_color_key$col, 
       pch=21, pt.cex = 1, cex = 1, 
       bty = "n", inset = .5)
dev.off()

# ---------------------
# Summarize info
# ---------------------

path_summary <-
  paths %>% 
  filter(!is.na(path) & !is.infinite(path)) %>%
  group_by(key) %>%
  summarise(sum = n(), mean = mean(path)) %>%
  mutate(percent = sum/count*100) %>%
  arrange(desc(percent))
write.csv(path_summary, "Files/Output/MI_path_summary.csv", row.names = FALSE)
write.csv(paths, "Files/Output/MI_path_full.csv", row.names = FALSE)

# count <- 188
MI_summary <-
  paths %>% 
  filter(key != "DSM_medRT" & key != "MinLag_0" & !is.na(path) & !is.infinite(path)) %>%
  group_by(key) %>%
  mutate(sumKey = n()) %>%
  mutate(percentKey = sumKey/count) %>%
  ungroup() %>%
  select(key, percentKey) %>%
  distinct(.)

MI_summary_plot <-
  MI_summary %>%
  ggplot(., aes(x = percentKey, y = reorder(key, percentKey))) +
  geom_bar(stat = "identity", fill = "lightblue4", width = .6) +
  theme_bw() +
  xlab("") +
  labs(x = "Fraction of participants", fill = "Path length", y = "") +
  guides(fill = 'none') +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))

ggsave("Files/Output/MI_summary_plot.tiff", MI_summary_plot, units = "in", width = 10, height = 5)

# ---------------------
# Read in images & combine
# ---------------------

# read in images
image <- list.files(main_output_dir, pattern = "\\.png$", full.names = TRUE)
images <- map(image, image_read)

count <- 0
num <- 40
for (i in 1:ceiling(length(images)/40)) {
  if (count + 40 > length(images)) {
    num <- length(images) - count 
    print(num)
  }
  
  par(mar=rep(0,4)) 
  images_joined <- image_join(images[(1 + count):(num + count)])
  
  # layout the plots into a matrix 
  layout(matrix(1:(ceiling(num/5)*5), ncol=5, byrow=TRUE))
  
  # do the plotting
  for(j in 1:num) {
    plot(NA,xlim=0:1,ylim=0:1,xaxt="n",yaxt="n",bty="n")
    rasterImage(images_joined[j],0,0,1,1)
  }
  
  dev.print(pdf, paste0(main_output_dir, "/comb_plot_", i, ".pdf"),
            height = 8, width = 5)
  
  count <- count + 40
}


