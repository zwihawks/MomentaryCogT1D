# --------------------------
# Step 4: Aim 1 hypothesis-driven analysis
# --------------------------

# ---------------------
# source functions & load libraries  
# ---------------------
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path ))
source("Functions.R")
load_packs()

main_output_dir <- "Files/Output" 
wrangle_for_cat_plots <- FALSE

# ---------------------
# Read in files
# ---------------------

# Load in HB models
files <- list.files(pattern = "to_model", recursive = TRUE)
files <- files[str_detect(files, "pct_inclusion") & 
                 !str_detect(files, "copy") &
                 !str_detect(files, "Archive") &
                 !str_detect(files, "rawTrue") &
                 !str_detect(files, "rawFalse")]
data_files <- files %>% map(readRDS) %>% setNames(files)
models <- data_files[str_detect(names(data_files), "row1.rds") | 
                       str_detect(names(data_files), "row2.rds")] 

# load EMA data
mydf_clean <- read.csv("Files/20221013_50pct_inclusion/mydf_clean_20221014.csv")

# TMB data
TMB_dat <- read.csv("Files/20221013_50pct_inclusion/TMB_dat_20221014.csv")
TMB_tmp <- read.csv("Files/20221013_50pct_inclusion/TMB_tmp_20221014.csv")

# demographic & clinic data - one row per participant
TMB_clean <-
  TMB_tmp %>%
  mutate(demo_gender_male = if_else(demo_gender == "male", 1, 0),
         demo_education_num = if_else(demo_education == "high", 0, NA_real_),
         demo_education_num = if_else(demo_education == "some_college", 1, demo_education_num),
         demo_education_num = if_else(demo_education == "technical", 1, demo_education_num),
         demo_education_num = if_else(demo_education == "college", 2, demo_education_num),
         demo_education_num = if_else(demo_education == "masters", 3, demo_education_num),
         demo_education_num = if_else(demo_education == "grad", 4, demo_education_num),
         demo_ethnicity_africanOrBlack = if_else(demo_ethnicity == "africanOrBlack", 1, 0),
         demo_ethnicity_europeanOrWhite = if_else(demo_ethnicity == "europeanOrWhite", 1, 0),
         demo_hispanic_yes = if_else(demo_hispanic == "yes", 1, 0),
         clinic_cgmusestat = if_else(clinic_CGMUseStat == "Current", 1, 0)) %>%
  {
    bind_cols(
      select_if(., is.numeric),
      select_at(., "user_id")) 
  } %>%
  select(-contains("utc_offset"), -contains("timeOn"), -contains("ump")) %>% 
  unique(.)

# cgm use count (relevant to blinding)
TMB_clean %>% 
  select(user_id, clinic_cgmusestat) %>% 
  distinct(.) %>% 
  summarise(sum = sum(clinic_cgmusestat, na.rm = TRUE)) 


# ---------------------
# Compute bayes R2
# ---------------------

df_R2 <- as.data.frame(matrix(0, length(models), 2))
for (i in 1:length(models)) {
  df_R2[i, 1] <- names(models)[i]
  df_R2[i, 2] <- median(bayes_R2(models[[i]]))
}
write.csv(df_R2, "Files/Output/bayes_R2.csv", row.names = FALSE)
df_R2 %>%
  separate(V1, into = c("omit", "inclusion", "outcome"), sep = "/") %>%
  group_by(outcome) %>%
  summarise(mean(V2))

# ---------------------
# Output model results
# ---------------------

prep_for_table <- 
  models %>%
  map_dfr(broom.mixed::tidy, conf.int = TRUE, 
          effects = c("ran_pars", "fixed"), .id = "Names") %>%
  separate(Names, into = c("omit1", "inclusion", "test"), sep = "/") %>% 
  separate(inclusion, into = c("omit3", "inclusion", "omit4"), sep = "_") %>% 
  separate(test, into = c("omit3", "omit4", "test", "outcome"), sep = "_") %>% 
  mutate(outcome = if_else(str_detect(outcome, "row1.rds"), "accuracy", "RT"),
         pct_test_outcome = paste(inclusion, test, outcome, sep = "_")) %>%
  select(-contains("omit"), -inclusion, -test, -outcome) 

prep_for_table %>%
  select(-pct_test_outcome) %>%
  kbl("html", digits = 3) %>%
  kable_classic(html_font = "Cambria", "basic", "center", full_width = F) %>%
  pack_rows(index = table(prep_for_table$pct_test_outcome)) %>%
  cat(., file = paste0(main_output_dir, "/Aim1_model_statistics.html"))

# ---------------------
# Distributions in 50% data
# ---------------------

dfs <- 
  data_files[str_detect(names(data_files), "[gcpt|dsm].rds")] %>% 
  bind_rows(., .id = "id") %>%
  separate(id, into = c("omit1", "inclusion", "omit2"), sep = "/") %>% 
  separate(inclusion, into = c("omit3", "inclusion", "omit4"), sep = "_") %>% 
  select(-contains("omit")) %>%
  mutate(test_outcome = paste0(test, ": ", outcome)) %>%
  filter((test == "dsm" & outcome == "num_correct") |
           (test == "dsm" & outcome == "medianRTc") |
           (test == "gradcpt" & outcome == "dprime") |
           (test == "gradcpt" & outcome == "medianRTc") ) 

cog_gl_distributions <-
  dfs %>%
  unnest(data) %>%
  select(inclusion, test, outcome, user_id, value, glucose = MinLag_0) %>%
  unite(test_outcome, c(test, outcome), sep = ": ") %>%
  {. ->> for_gl} %>%
  select(inclusion, test_outcome, user_id, value) %>%
  distinct(.) %>%
  bind_rows(for_gl %>% select(-value) %>% 
              mutate(test_outcome = "cgm") %>% 
              rename(value = glucose) %>% 
              distinct(.)) %>%
  ggplot(., aes(x = value)) +
  geom_histogram(aes(fill = inclusion), alpha = .5, 
                 position = "identity", bins = 15) +
  facet_wrap(~ test_outcome, scales = "free", ncol = 5) +
  theme_bw() +
  labs(y = "Count", x = "", fill = "EMA\ncompletion")

ggsave("Files/Output/Aim1_distributions.tiff", 
       cog_gl_distributions, 
       units = "in", width = 12, height = 3)

# ---------------------
# Population-level effects
# ---------------------
# First, plot posterior intervals ("omnibus" test of sorts)
fixef <-
  models %>%
  map(summary, pars = "beta", probs = c(.025, .05, .95, .975)) %>%
  map(as.data.frame) %>%
  map(rownames_to_column) %>%
  bind_rows(.id = "id") %>%
  separate(id, into = c("omit0", "inclusion", "test"), sep = "/") %>% 
  separate(inclusion, into = c("omit1", "inclusion", "omit2"), sep = "_") %>% 
  separate(test, into = c("omit3", "omit4", "test", "outcome"), sep = "_") %>% 
  select(-contains("omit")) %>%
  rename(fixef = rowname) %>%
  mutate(sig = if_else((`5%` > 0 & `95%` > 0) |
                         (`5%` < 0 & `95%` < 0), "^", ""),
         sig = if_else((`2.5%` > 0 & `97.5%` > 0) |
                         (`2.5%` < 0 & `97.5%` < 0), "*", sig))

# refer to posterior_intervals documentation for rationale for using 90% CI
post_plots <-
  fixef %>% 
  mutate(jitter = 0,
         jitter = if_else(inclusion == "50pct", -.02, jitter),
         jitter = if_else(inclusion == "80pct", .02, jitter),
         fixef_num = if_else(fixef == "poly(MinLag_0_WP, 2, raw = FALSE)1", 1, 1.2),
         `90% CI` = if_else(fixef == "poly(MinLag_0_WP, 2, raw = FALSE)1", "Linear", "Quadratic"),
         `95% CI` = if_else(fixef == "poly(MinLag_0_WP, 2, raw = FALSE)1", "Linear", "Quadratic"),
         outcome = if_else(str_detect(outcome, "row1"), "Accuracy", "Speed (ms)"),
         test = toupper(test)) %>%
  posterior_intervals(.) +
  facet_wrap(test ~ outcome, scales = "free") +
  theme(axis.text.y=element_blank(),  #remove y axis labels
        axis.ticks.y=element_blank())  #remove y axis ticks

post_plots_DSMRT <-
  fixef %>% 
  filter(test == "dsm" & str_detect(outcome, "row2")) %>%
  mutate(jitter = 0,
         jitter = if_else(inclusion == "50pct", -.02, jitter),
         jitter = if_else(inclusion == "80pct", .02, jitter),
         fixef_num = if_else(fixef == "poly(MinLag_0_WP, 2, raw = FALSE)1", 1, 1.2),
         `90% CI` = if_else(fixef == "poly(MinLag_0_WP, 2, raw = FALSE)1", "Linear", "Quadratic"),
         `95% CI` = if_else(fixef == "poly(MinLag_0_WP, 2, raw = FALSE)1", "Linear", "Quadratic"),
         outcome = if_else(str_detect(outcome, "row1"), "Accuracy", "Speed (ms)"),
         test = toupper(test)) %>%
  posterior_intervals(.) +
  theme(axis.text.y=element_blank(),  #remove y axis labels
        axis.ticks.y=element_blank(), #remove y axis ticks
        legend.position = "left")  

ggsave("Files/Output/posterior_interval_plots.tiff", post_plots, 
       units = "in", height = 6, width = 9)

# ---------------------
# Individual-level effects
# ---------------------
# caterpillar plots across different inclusion criteria
rm(fixef, data_files)

if (wrangle_for_cat_plots) {
  iterate <- c("dsm_row1", "dsm_row2", "gcpt_row1", "gcpt_row2")
  for (i in 1:length(iterate)) {
    
    test_outcome <- iterate[i]
    wrangle_cat <- wrangle_func(models = models, test_outcome = test_outcome)
    
    if (i == 1) {
      cat_plots <- cat_plot_func(wrangle_cat)
      summary_stats <- summary_stats_func(wrangle_cat)
      percent_strong <- percent_strong_func(wrangle_cat)
      rm(wrangle_cat)
    } else {
      tmp <- cat_plot_func(wrangle_cat)
      summary_stats_tmp <- summary_stats_func(wrangle_cat)
      percent_strong_tmp <- percent_strong_func(wrangle_cat)
      rm(wrangle_cat)
      
      cat_plots <- bind_rows(cat_plots, tmp)
      summary_stats <- bind_rows(summary_stats, summary_stats_tmp)
      percent_strong <- bind_rows(percent_strong, percent_strong_tmp)
      rm(tmp, summary_stats_tmp, percent_strong_tmp)
    }
    print(i)
  }
  
  saveRDS(cat_plots, "Files/Output/for_caterpillar_plots.Rds")
  write.csv(summary_stats, "Files/Output/iLevelEffects_summaryStats.csv", row.names = FALSE)
  write.csv(percent_strong, "Files/Output/iLevelEffects_isStrong.csv", row.names = FALSE)
} else {
  print("Reading in files...")
  cat_plots <- readRDS("Files/Output/for_caterpillar_plots.Rds")
  summary_stats <- read.csv("Files/Output/iLevelEffects_summaryStats.csv")
  percent_strong <- read.csv("Files/Output/iLevelEffects_isStrong.csv")
}

percent_strong %>% 
  group_by(outcome) %>% 
  mutate(total = n()) %>%
  ungroup() %>%
  group_by(outcome, total) %>%
  summarise(sum = sum(color)) %>%
  mutate(percent = sum/total)

# Saving to file
cat_plots_comb <-
  ggpubr::ggarrange(ggpubr::ggarrange(cat_plots$plots[[4]], cat_plots$plots[[3]], nrow = 1, legend = "none"),
                    ggpubr::ggarrange(cat_plots$plots[[2]], cat_plots$plots[[1]], nrow = 1,
                      common.legend = TRUE, legend = "bottom"), 
            nrow = 2, heights = c(.45, .55)) 
ggsave("Files/Output/iLevel_plots_comb.tiff", cat_plots_comb, 
       units = "in", height = 8, width = 16)


quad_plot <-
  quad_plots(models[str_detect(names(models), "dsm_row2")], 
             EMA_cutoffs = c("66pct"), 
             TMB_list = list(NULL),
             lasso = FALSE,
             lasso_results = NULL) 

tmp_legend <- 
  get_legend(
    post_plots_DSMRT + 
      theme(legend.margin = margin(0, 0, 0, 0),
            legend.spacing.y = unit(2, "mm")) ) 
tmp_post_plots_DSMRT <- post_plots_DSMRT + theme(legend.position = 'none')

# ---------------------
# Combine to create Aim 1 figure
# ---------------------
Aim1_comb <-
  cowplot::plot_grid(
  cowplot::plot_grid(
    tmp_legend, tmp_post_plots_DSMRT, quad_plot[[2]][[1]], nrow = 1, 
    rel_widths = c(.1, .35, .45),
    labels = c("A", "", "B")),
  (cat_plots$plots[[2]] + ggtitle("") + 
     facet_wrap(~inclusion, nrow = 1, scales = "free_x") +
     geom_hline(yintercept = 0, lty = 2) +
     theme(legend.position = "left",
           legend.margin = margin(0, 0, 0, 0),
           legend.spacing.y = unit(2, "mm")) ), 
  nrow = 2, align = "h", axis = "r",
  rel_heights = c(.55, .45),
  labels = c("", "C")) 

ggsave("Files/Output/Aim1_comb.tiff", Aim1_comb, 
       units = "in", height = 8, width = 16)

# -----------------------
# Evaluate optimal performance
# -----------------------

epred_summary <- quad_plots(models[str_detect(names(models), "dsm_row2")], 
                            EMA_cutoffs = c("66pct", "80pct", "50pct"), 
                            TMB_list = list(NULL),
                            lasso = FALSE,
                            lasso_results = NULL,
                            n_reps = NA,
                            make_plots = FALSE)

pred_optimal_plotv1 <- plot_min_max(epred_summary, optimal = "min", text = "DSM RT")

ggsave("Files/Output/optimal_performance.tiff", pred_optimal_plotv1[[1]], 
       units = "in", height = 4, width = 7)

# optimal performance occurred 0.78 standard deviations above individuals’ glucose means
# was associated with a 0.66% (6.19 ms) performance gain
summary_stats <-
  pred_optimal_plotv1[[2]] %>% 
  select(contains("deviation"), MinLag_0_WP, user_id, inclusion) %>%
  # group by inclusion to computed weighted avg.
  group_by(inclusion) %>%
  summarise_if(is.numeric, mean, na.rm = TRUE) %>%
  ungroup() %>%
  summarise_if(is.numeric, mean, na.rm = TRUE)

# Performance decrements exceeded 1% (9.76 ms) relative to individuals’ cognitive means 
# when glucose was 1.33 standard deviations below and 2.82 standard deviations 
# above individuals’ glucose means
summary_stats_1pct_dec <-
  epred_summary %>% 
  mutate(percent_deviation = 100*((mean - user_intercept)/user_intercept),
         deviation = mean - user_intercept,
         tmp_deviation = percent_deviation - 1,
         tmp_code = if_else(MinLag_0_WP < 0, "negative", "positive")) %>% 
  group_by(MinLag_0_WP, inclusion, tmp_code) %>% 
  summarise_at(vars(contains("deviation")), mean) %>%
  ungroup() %>%
  group_by(inclusion, tmp_code) %>%
  filter(tmp_deviation > 0) %>%
  filter(tmp_deviation == min(tmp_deviation)) %>%
  ungroup() %>%
  group_by(tmp_code) %>%
  summarise_if(is.numeric, mean) %>%
  ungroup() %>%
  mutate(mean_deviation = mean(deviation))

# Mean & SD glucose
user_gl_sd <- 
  models$`Files/20221013_50pct_inclusion/to_model_dsm_row2.rds`$data %>%
  group_by(user_id) %>%
  summarise(MinLag_0_sd = sd(MinLag_0, na.rm = TRUE),
            MinLag_0_mean = mean(MinLag_0, na.rm = TRUE))

user_gl_sd %>%
  summarise_if(is.numeric, mean, na.rm = TRUE)

# gl deviation in mg/dl
# as well: raw gl value associated with 1% performance decrement
pred_optimal_plotv1[[2]] %>%
  left_join(user_gl_sd) %>%
  mutate(MinLag_0_centered = MinLag_0_WP*MinLag_0_sd,
         MinLag_0_raw = MinLag_0_centered + MinLag_0_mean,
         MinLag_0_SDabove = MinLag_0_raw + 
           (summary_stats_1pct_dec %>% filter(MinLag_0_WP > 0) %>% pull(MinLag_0_WP))*MinLag_0_sd,
         MinLag_0_SDbelow = MinLag_0_raw + 
           (summary_stats_1pct_dec %>% filter(MinLag_0_WP < 0) %>% pull(MinLag_0_WP))*MinLag_0_sd) %>%
  # group by inclusion to compute weighted avg.
  group_by(inclusion) %>%
  summarise(mean_gl_deviation = mean(MinLag_0_centered),
            mean_gl_raw = mean(MinLag_0_raw),
            mean_gl_SDabove = mean(MinLag_0_SDabove),
            mean_gl_SDbelow = mean(MinLag_0_SDbelow)) %>%
  ungroup() %>%
  summarise(mean_gl_deviation = mean(mean_gl_deviation),
            mean_gl_raw = mean(mean_gl_raw),
            mean_gl_SDabove = mean(mean_gl_SDabove),
            mean_gl_SDbelow = mean(mean_gl_SDbelow))

# correlation between age & dsm RT
tmp <-
  TMB_clean %>%
  select(user_id, dsm_medianRTc, demo_age) %>%
  distinct(.) %>%
  group_by(user_id, demo_age) %>%
  summarise(mean_RT = mean(dsm_medianRTc, na.rm = TRUE)) %>%
  ungroup()

cor.test(tmp$mean_RT, tmp$demo_age)




