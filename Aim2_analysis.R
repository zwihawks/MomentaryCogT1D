# --------------------------
# Step 5: Aim 2 data-driven inter-individual difference analysis
# --------------------------

# ---------------------
# source functions & load libraries  
# ---------------------
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path ))
source("Functions.R")
load_packs()

main_output_dir <- "Files/Output" 
run_lasso <- FALSE
lasso_distributions <- FALSE

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
subject_level_data <- 
  read.csv("Files/subject_level_data.csv") %>%
  wrangle_subject_file(., select_SB = FALSE) %>%
  mutate(clinic_STOP_BANG_risk_num = if_else(clinic_STOP_BANG_risk == "low", 0, 1),
         clinic_STOP_BANG_risk_num = if_else(clinic_STOP_BANG_risk == "high", 2, clinic_STOP_BANG_risk_num)) %>%
  select(-clinic_BANG_sum, -clinic_STOP_sum, -clinic_STOP_BANG_risk)
  
TMB_clean <-
  TMB_tmp %>%
  select(-c(clinic_Weight, clinic_Height, clinic_WaistCir, 
            clinic_NeckCir, clinic_SHSeizComaLast12MonthsB)) %>%
  mutate(demo_education_num = if_else(demo_education == "high", 0, NA_real_),
         demo_education_num = if_else(demo_education == "some_college", 1, demo_education_num),
         demo_education_num = if_else(demo_education == "technical", 1, demo_education_num),
         demo_education_num = if_else(demo_education == "college", 2, demo_education_num),
         demo_education_num = if_else(demo_education == "masters", 3, demo_education_num),
         demo_education_num = if_else(demo_education == "grad", 4, demo_education_num),
         demo_ethnicity_africanOrBlack = if_else(demo_ethnicity == "africanOrBlack", 1, 0),
         demo_ethnicity_europeanOrWhite = if_else(demo_ethnicity == "europeanOrWhite", 1, 0),
         demo_hispanic_yes = if_else(demo_hispanic == "yes", 1, 0)) %>%
  {
    bind_cols(
      select_if(., is.numeric),
      select_at(., "user_id"))
  } %>%
  select(-contains("utc_offset"), -contains("timeOn"), -contains("ump")) %>%
  unique(.) %>%
  # left join clinic info
  left_join(subject_level_data)

# ---------------------
# Lasso regression
# ---------------------

dsm_speed <- models[str_detect(names(models), "dsm_row2")]
n_reps <- 1000
if (run_lasso) {
  dsm_model_terms <-
    dsm_speed %>%
    map(spread_draws, `poly(MinLag_0_WP, 2, raw = FALSE)2`, b[group, term], sep = " u") %>%
    bind_rows(.id = "id") %>%
    filter(group == "poly(MinLag_0_WP, 2, raw = FALSE)2") %>%
    mutate(condition_mean = `poly(MinLag_0_WP, 2, raw = FALSE)2` + b) %>%
    group_by(id, term) %>%
    median_qi(condition_mean) %>%
    separate(term, into = c("omit", "user_id"), sep = ":") %>%
    select(-omit) %>% 
    select(id, user_id, condition_mean)
  
  TMB_plus <-
    TMB_clean %>% 
    group_by(user_id) %>%
    mutate(EMA_max = max(EMA_count, na.rm = TRUE)) %>%
    ungroup() %>%
    select(user_id, contains("demo_"), contains("clinic_")) %>%
    {
      bind_cols(
        select_if(., is.numeric),
        select_at(., "user_id")) 
    } %>%
    unique(.) %>%
    # right join so that we have one row per user, inclusion criterion 
    right_join(dsm_model_terms) %>%
    select_if(~(sum(!is.na(.)) > .5*sum(is.na(.) | !is.na(.)))) 
  
  if (lasso_distributions) {
    lasso_distributions <- 
      TMB_plus %>% 
      filter(!is.na(condition_mean)) %>%
      select(-clinic_id, -condition_mean) %>%
      mutate(id_cat = if_else(str_detect(id, "80pct"), "80pct", "50pct"),
             id_cat = if_else(str_detect(id, "66pct"), "66pct", id_cat)) %>%
      gather(key = key, value = value, -user_id, -id, -id_cat) %>%
      mutate(key = str_replace_all(key, "ethnicity", "race")) %>%
      ggplot(., aes(x = value)) +
      geom_histogram(aes(fill = id_cat), alpha = .5, 
                     position = "identity", bins = 15) +
      facet_wrap(~ key, scales = "free", ncol = 4) +
      theme_bw() +
      labs(y = "Count", x = "", fill = "EMA\ncompletion")
    
    ggsave("Files/Output/Aim2_distributions.tiff", 
           lasso_distributions, 
           units = "in", width = 12, height = 15)
  }
  
  for (j in 1:n_reps) {
    for (i in c("80pct", "66pct", "50pct")) {
      for_lasso <- 
        TMB_plus %>% filter(str_detect(id, i)) %>% 
        filter(!is.na(condition_mean)) %>%
        mutate_all( ~ ifelse(is.na(.x), median(.x, na.rm = T), .)) %>%
        mutate_if(is.numeric, funs(z = (. - mean(.)) / sd(.)))
      x_df <- for_lasso %>% select(contains("_z"), -c(clinic_id_z, condition_mean_z))
      cv_model <- cv.glmnet(y = for_lasso$condition_mean, x = data.matrix(x_df), alpha = 1)
      coefs <- coef(cv_model, s = "lambda.min")
      
      if (i == "80pct") {
        lasso_results_inner <-
          data.frame(var = coefs@Dimnames[[1]][coefs@i + 1], 
                     coef = coefs@x,
                     inclusion = i,
                     rep = j) 
        
        cvm <- 
          data.frame(lambda = cv_model$lambda,
                     mse = cv_model$cvm,
                     inclusion = i,
                     rep = j) %>%
          mutate(lambda_min = if_else(lambda == cv_model$lambda.min, 1, 0),
                 lambda_1se = if_else(lambda == cv_model$lambda.1se, 1, 0))
      } else {
        tmp1 <-
          data.frame(var = coefs@Dimnames[[1]][coefs@i + 1], 
                     coef = coefs@x,
                     inclusion = i,
                     rep = j) 
        lasso_results_inner <- bind_rows(lasso_results_inner, tmp1)
        
        tmp2 <- 
          data.frame(lambda = cv_model$lambda,
                     mse = cv_model$cvm,
                     inclusion = i,
                     rep = j) %>%
          mutate(lambda_min = if_else(lambda == cv_model$lambda.min, 1, 0),
                 lambda_1se = if_else(lambda == cv_model$lambda.1se, 1, 0))
        cvm <- bind_rows(cvm, tmp2)
      }
    }
    if (j == 1) {
      lasso_results_outer <- lasso_results_inner
      cvm_outer <- cvm
    } else {
      lasso_results_outer <- bind_rows(lasso_results_outer, lasso_results_inner)
      cvm_outer <- bind_rows(cvm_outer, cvm)
    }
  }
  write.csv(lasso_results_outer, "Files/Output/Lasso_results.csv", row.names = FALSE)
  write.csv(cvm_outer, "Files/Output/cvm_results.csv", row.names = FALSE)
  write.csv((x_df %>% colnames(.) %>% str_remove_all(., "_z") %>% as.data.frame()),
            "Files/Output/Lasso_features.csv")
} else {
  lasso_results_outer <- read.csv("Files/Output/Lasso_results.csv")
  cvm_outer <- read.csv("Files/Output/cvm_results.csv")
}

lasso_results <- 
  lasso_results_outer %>%
  group_by(var, inclusion) %>%
  summarise(coef_mean = mean(coef, na.rm = TRUE),
            coef_sd = sd(coef, na.rm = TRUE),
            coef_n = n(),
            coef_lower = coef_mean - 1.96*coef_sd/sqrt(coef_n),
            coef_upper = coef_mean + 1.96*coef_sd/sqrt(coef_n)) 

mean_RMSE <- 
  cvm_outer %>%
  filter(lambda_min == 1) %>%
  group_by(inclusion) %>%
  summarise(mean_mse = mean(mse),
            sd_mse = sd(mse)) %>%
  mutate(mean_rmse = sqrt(mean_mse),
         sd_rmse = sqrt(sd_mse))

# Require that coefficient is retained in > 50% of reps
lasso_flex <- 
  lasso_results %>%
  filter((coef_n >= n_reps*.5) & var != ("(Intercept)")) %>%
  mutate(`Mean (SD, n)` = paste0(round(coef_mean, 2), 
                                 " (", round(coef_sd, 2), ", ",coef_n, ")"),
         inclusion = paste0(inclusion, ": Mean (SD, n)")) %>%
  select(Variable = var, inclusion, `Mean (SD, n)`) %>%
  pivot_wider(names_from = inclusion, values_from = `Mean (SD, n)`)

lasso_tmp <- lasso_flex
lasso_tmp[1,] <-
  as.list(
    c("mean_RMSE", 
      paste0(as.character(round(mean_RMSE$mean_rmse[1], 2)), " (", 
             as.character(round(mean_RMSE$sd_rmse[1], 2)), ", NA)"), 
      paste0(as.character(round(mean_RMSE$mean_rmse[2], 2)), " (", 
             as.character(round(mean_RMSE$sd_rmse[2], 2)), ", NA)"), 
      paste0(as.character(round(mean_RMSE$mean_rmse[3], 2)), " (", 
             as.character(round(mean_RMSE$sd_rmse[3], 2)), ", NA)")))

lasso_flex <-
  bind_rows(lasso_tmp[1,], lasso_flex) 
options(knitr.kable.NA = '')
lasso_flex %>%
  kbl("html", digits = 3) %>%
  kable_classic(html_font = "Cambria", "basic", "center", full_width = F) %>%
  row_spec(1, background = "lightgray") %>%
  cat(., file = paste0(main_output_dir, "/lasso_reg_table.html"))

lass_plot <-
  ggplot(lasso_results %>% filter(var != "(Intercept)"),
         aes(x = coef_mean, y = var, group = inclusion)) +
  geom_point(aes(fill = inclusion, size = coef_n), shape = 21, alpha = .6) +
  theme_bw() +
  geom_vline(xintercept = 0, lty = 2) +
  labs(y = "", x = "Coefficient", title = "BP predictors of quadratic RE using lasso regression for feature selection", 
       subtitle = "Plotting features retained by at least one model across three EMA inclusion criteria",
       fill = "EMA inclusion\ncut-off",
       size = paste0("Model inclusion\n(n out of ", n_reps, ")"))
ggsave("Files/Output/lasso_plot.tiff", lass_plot, units = "in", height = 5, width = 9)

# Plotting individual and group-average quadratic curves
quad_plot_gradient <-
  quad_plots(model_input = dsm_speed, 
             EMA_cutoffs = c("66pct"), 
             TMB_list = list(NULL, TMB_clean, TMB_clean, TMB_clean, TMB_clean, TMB_clean),
             lasso = TRUE,
             lasso_results = lasso_results,
             n_reps = n_reps*.5,
             make_plots = TRUE) 

quad_plot_gradient_comb <-
  cowplot::plot_grid(
    cowplot::plot_grid(
      quad_plot_gradient[[2]][[2]], 
      quad_plot_gradient[[2]][[3]],
      quad_plot_gradient[[2]][[4]],
      nrow = 1, labels = c("A","B", "C")),
    cowplot::plot_grid(
      NULL,
      quad_plot_gradient[[2]][[5]], 
      quad_plot_gradient[[2]][[6]],
      NULL,
      rel_widths = c(.15, .35, .35, .15),
      nrow = 1, labels = c("", "D","E", "")),
    cowplot::plot_grid(
      quad_plot_gradient[[1]][[2]], 
      quad_plot_gradient[[1]][[3]],
      quad_plot_gradient[[1]][[4]],
      quad_plot_gradient[[1]][[5]],
      quad_plot_gradient[[1]][[6]],
      nrow = 1),
    nrow = 3, align = 'v', axis = 'r', rel_heights = c(.45, .45, .1))
ggsave("Files/Output/quad_plot_gradient_comb.tiff", 
       quad_plot_gradient_comb, units = "in", height = 9, width = 14)



