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
  mutate(clinic_STOP_BANG_risk_num = if_else(clinic_STOP_BANG_risk == "low", 0, 1)) %>%
  mutate(clinic_STOP_BANG_risk_num = if_else(clinic_STOP_BANG_risk == "high", 2, clinic_STOP_BANG_risk_num)) %>%
  select(-clinic_BANG_sum, -clinic_STOP_sum, -clinic_STOP_BANG_risk)

TMB_clean <-
  TMB_tmp %>%
  select(-c(clinic_Weight, clinic_Height, clinic_WaistCir, 
            clinic_NeckCir, clinic_SHSeizComaLast12MonthsB,
            clinic_gluAbsRateChgMean, clinic_gluSDRateChg,
            clinic_gluGeoMean, clinic_gluGeoSD)) %>%
  mutate(demo_education_num = if_else(demo_education == "high", 0, NA_real_)) %>%
  mutate(demo_education_num = if_else(demo_education == "some_college", 1, demo_education_num)) %>%
  mutate(demo_education_num = if_else(demo_education == "technical", 1, demo_education_num)) %>%
  mutate(demo_education_num = if_else(demo_education == "college", 2, demo_education_num)) %>%
  mutate(demo_education_num = if_else(demo_education == "masters", 3, demo_education_num)) %>%
  mutate(demo_education_num = if_else(demo_education == "grad", 4, demo_education_num)) %>%
  mutate(demo_ethnicity_africanOrBlack = if_else(demo_ethnicity == "africanOrBlack", 1, 0)) %>%
  mutate(demo_ethnicity_europeanOrWhite = if_else(demo_ethnicity == "europeanOrWhite", 1, 0)) %>%
  mutate(demo_hispanic_yes = if_else(demo_hispanic == "yes", 1, 0)) %>%
  mutate_at(vars("clinic_gluInRange", "clinic_gluBelow70", "clinic_gluBelow54", 
                 "clinic_gluAbove180", "clinic_gluAbove250", "clinic_gluCV"),
            funs(as.numeric(str_remove(., "%")))) %>%
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
    map(spread_draws, `(Intercept)`, b[group, term], sep = " u") %>%
    bind_rows(.id = "id") %>%
    filter(group == "(Intercept)") %>%
    mutate(condition_mean = `(Intercept)` + b) %>%
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
  
  # rename variables for interpretability
  colnames(TMB_plus) <- colnames(TMB_plus) %>% str_remove_all(., "demo_") %>% str_remove_all(., "clinic_")
  colnames(TMB_plus)[c(37, 38, 50:57)] <- 
    c("clinic_id", "SevereHypoEvents", "snoring_binary", "tired_binary", "gasping_binary", "bloodPressure_binary",
      "BMI_binary", "age_binary", "NeckCir_binary", "gender_male")
  
  vec <- c("80pct", "66pct", "50pct")
  for (j in 1:n_reps) {
    for (i in 1:length(vec)) {
      for_lasso <- 
        TMB_plus %>% filter(str_detect(id, vec[i])) %>% 
        filter(!is.na(condition_mean)) %>%
        mutate_all( ~ ifelse(is.na(.x), median(.x, na.rm = T), .)) %>%
        mutate_if(is.numeric, funs(z = (. - mean(.)) / sd(.)))
      x_df <- for_lasso %>% select(contains("_z"), -c(clinic_id_z, condition_mean_z))
      cv_model <- cv.glmnet(y = for_lasso$condition_mean, x = data.matrix(x_df), alpha = 1)
      coefs <- coef(cv_model, s = "lambda.min")
      
      if (i == 1) {
        lasso_results_inner <-
          data.frame(var = coefs@Dimnames[[1]][coefs@i + 1], 
                     coef = coefs@x,
                     inclusion = vec[i],
                     rep = j) 
        
        cvm <- 
          data.frame(lambda = cv_model$lambda,
                     mse = cv_model$cvm,
                     rsq = 1 - (cv_model$cvm/var(for_lasso$condition_mean)),
                     inclusion = vec[i],
                     rep = j) %>%
          mutate(lambda_min = if_else(lambda == cv_model$lambda.min, 1, 0),
                 lambda_1se = if_else(lambda == cv_model$lambda.1se, 1, 0))
      } else {
        tmp1 <-
          data.frame(var = coefs@Dimnames[[1]][coefs@i + 1], 
                     coef = coefs@x,
                     inclusion = vec[i],
                     rep = j) 
        lasso_results_inner <- bind_rows(lasso_results_inner, tmp1)
        
        tmp2 <- 
          data.frame(lambda = cv_model$lambda,
                     mse = cv_model$cvm,
                     rsq = 1 - (cv_model$cvm/var(for_lasso$condition_mean)),
                     inclusion = vec[i],
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
      print(paste0(vec[i], ", ", j))
    }
  }
  write.csv(lasso_results_outer, "Files/Output/Lasso_results_intercept.csv", row.names = FALSE)
  write.csv(cvm_outer, "Files/Output/cvm_results_intercept.csv", row.names = FALSE)
  write.csv((x_df %>% colnames(.) %>% str_remove_all(., "_z") %>% as.data.frame()),
            "Files/Output/Lasso_features_intercept.csv")
} else {
  lasso_results_outer <- read.csv("Files/Output/Lasso_results_intercept.csv")
  cvm_outer <- read.csv("Files/Output/cvm_results_intercept.csv")
}

lasso_results <- 
  lasso_results_outer %>%
  group_by(var, inclusion) %>%
  dplyr::summarise(coef_mean = mean(coef, na.rm = TRUE),
                   coef_sd = sd(coef, na.rm = TRUE),
                   coef_n = n(),
                   coef_lower = coef_mean - 1.96*coef_sd/sqrt(coef_n),
                   coef_upper = coef_mean + 1.96*coef_sd/sqrt(coef_n)) 

mean_RMSE <- 
  cvm_outer %>%
  filter(lambda_min == 1) %>%
  group_by(inclusion) %>%
  summarise(mean_mse = mean(mse),
            sd_mse = sd(mse),
            mean_rsq = mean(rsq),
            sd_rsq = sd(rsq)) %>%
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
lasso_tmp[2,] <-
  as.list(
    c("mean_R2", 
      paste0(as.character(round(mean_RMSE$mean_rsq[1], 2)), " (", 
             as.character(round(mean_RMSE$sd_rsq[1], 2)), ", NA)"), 
      paste0(as.character(round(mean_RMSE$mean_rsq[2], 2)), " (", 
             as.character(round(mean_RMSE$sd_rsq[2], 2)), ", NA)"), 
      paste0(as.character(round(mean_RMSE$mean_rsq[3], 2)), " (", 
             as.character(round(mean_RMSE$sd_rsq[3], 2)), ", NA)")))

lasso_flex <-
  bind_rows(lasso_tmp[1:2,], lasso_flex) 
options(knitr.kable.NA = '')
lasso_flex %>%
  kbl("html", digits = 3) %>%
  kable_classic(html_font = "Cambria", "basic", "center", full_width = F) %>%
  row_spec(1:2, background = "lightgray") %>%
  cat(., file = paste0(main_output_dir, "/lasso_reg_table_intercept.html"))

