# --------------------------
# Step 6: Post-hoc tests
# --------------------------

# ---------------------
# source functions & load libraries  
# ---------------------
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path ))
source("Functions.R")
load_packs()

# ---------------------
# Read in files
# ---------------------

main_output_dir <- "Files/Output" 
lasso_results <- read.csv(paste0(main_output_dir, "/Lasso_results.csv"))
TMB_tmp <- read.csv("Files/20221013_66pct_inclusion/TMB_tmp_20221014.csv")
models_c <- readRDS("Files/20221013_66pct_inclusion/testCentering_rawFalse/to_model_dsm_row2_cIV.rds")
dsm_speed <- readRDS("Files/20221013_66pct_inclusion/to_model_dsm_row2.rds")
models_comb <- list(dsm_speed = dsm_speed, models_cIV = models_c)

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
# Output model results
# ---------------------

prep_for_table <- 
  broom.mixed::tidy(models_c, conf.int = TRUE, 
                    effects = c("ran_pars", "fixed")) %>%
  mutate(pct_test_outcome = "66pct_dsm_RT") 

prep_for_table %>%
  select(-pct_test_outcome) %>%
  kbl("html", digits = 3) %>%
  kable_classic(html_font = "Cambria", "basic", "center", full_width = F) %>%
  pack_rows(index = table(prep_for_table$pct_test_outcome)) %>%
  cat(., file = paste0(main_output_dir, "/posthoc_model_statistics.html"))


# ---------------------
# Population-level effects
# ---------------------

summary(models_c, pars = "beta", probs = c(.025, .05, .95, .975))

# ---------------------
# Individual-level effects
# ---------------------

test_outcome <- names(models_comb)
for (i in 1:length(test_outcome)) {
  test_outcome_tmp <- test_outcome[i]
  if (str_detect(test_outcome_tmp, "cIV")) {
    wrangle_cat1 <-
      models_comb[str_detect(names(models_comb), test_outcome_tmp)] %>%
      map(spread_draws, `poly(MinLag_0_centered, 2, raw = FALSE)2`, b[group, term], sep = " u") %>%
      bind_rows(.id = "id") %>%
      rename(MinLag_0_WP_fixed = `poly(MinLag_0_centered, 2, raw = FALSE)2`) %>%
      pivot_wider(names_from = group, values_from = b) %>%
      mutate(term = str_remove(term, "ser_id:")) %>%
      mutate(condition_mean_quad = MinLag_0_WP_fixed + `poly(MinLag_0_centered, 2, raw = FALSE)2`) 
  } else {
    wrangle_cat1 <-
      models_comb[str_detect(names(models_comb), test_outcome_tmp)] %>%
      map(spread_draws, `poly(MinLag_0_WP, 2, raw = FALSE)2`, b[group, term], sep = " u") %>%
      bind_rows(.id = "id") %>%
      rename(MinLag_0_WP_fixed = `poly(MinLag_0_WP, 2, raw = FALSE)2`) %>%
      pivot_wider(names_from = group, values_from = b) %>%
      mutate(term = str_remove(term, "ser_id:")) %>%
      mutate(condition_mean_quad = MinLag_0_WP_fixed + `poly(MinLag_0_WP, 2, raw = FALSE)2`) 
  }
  wrangle_cat2 <- 
    wrangle_cat1 %>%
    group_by(id, term) %>%
    median_qi(condition_mean_quad, .width = c(.95, .9, .66)) %>%
    ungroup() %>%
    arrange(condition_mean_quad) %>%
    rownames_to_column() %>%
    mutate(rep = i)
  
  if (i == 1) {
    wrangle_cat <- wrangle_cat2
  } else {
    wrangle_cat <- bind_rows(wrangle_cat, wrangle_cat2)
  }
}

for_cor <-
  wrangle_cat %>% 
  filter(.width == .95) %>%
  select(rowname, term, rep) %>% 
  pivot_wider(names_from = rep, values_from = rowname) %>%
  mutate(across(.cols = -c(term), .fns = ~as.numeric(.)))
colnames(for_cor)[2:3] <- test_outcome
cor_orig_cIV <- cor.test(for_cor$dsm_speed, for_cor$models_cIV)

# examine H2: variation in individual estimates of gl acceleration
model <- names(models_comb)[2]
sd_y <- models_comb[[2]]$data$value %>% sd(., na.rm = TRUE)*.1
equivalence <- 
  equivalence_test(models_comb[[2]], ci = c(.9, .95), effects = "random", component = "all", range = c(0, sd_y*2)) %>%
  as.data.frame(.) %>% 
  filter(str_detect(Parameter, "igma") & 
           !str_detect(Parameter, "ntercept") & 
           str_detect(Parameter, "\\)2") & 
           !str_detect(Parameter, "\\)1")) %>% 
  select(-Parameter)

equivalence_df <- equivalence %>% select(CI, H0 = ROPE_Equivalence, `inside ROPE` = ROPE_Percentage, HDI_low, HDI_high)
equivalence_df %>%
  kbl("html", digits = 3) %>%
  kable_classic(html_font = "Cambria", "basic", "center", full_width = F) %>%
  cat(., file = paste0(main_output_dir, "/posthoc_H2.html"))

# -------------------
# Lasso regression
#--------------------

wrangle_cat_lasso <-
  wrangle_cat %>%
  filter(rep == 2) %>%
  dplyr::select(id, user_id = term, condition_mean = condition_mean_quad) 

TMB_plus <-
  TMB_clean %>% 
  group_by(user_id) %>%
  mutate(EMA_max = max(EMA_count, na.rm = TRUE)) %>%
  ungroup() %>%
  select(user_id, contains("demo"), contains("clinic")) %>%
  {
    bind_cols(
      select_if(., is.numeric),
      select_at(., "user_id")) 
  } %>%
  unique(.) %>%
  # right join so that we have one row per user, inclusion criterion 
  right_join(wrangle_cat_lasso) %>%
  select_if(~(sum(!is.na(.)) > .5*sum(is.na(.) | !is.na(.)))) 

# lasso results from primary analyses
cvm_outer <- read.csv("Files/Output/cvm_results.csv")
l_seq <- cvm_outer %>% filter(lambda_min == 1) %>% pull(lambda) %>% quantile(.)
l_seq_for_glm <- seq(l_seq[[4]], l_seq[[5]], by = .1)

n_reps <- 1000
for (j in 1:n_reps) {
  for_lasso <- 
    TMB_plus %>% 
    filter(!is.na(condition_mean)) %>%
    mutate_all( ~ ifelse(is.na(.x), median(.x, na.rm = T), .)) %>%
    mutate_if(is.numeric, funs(z = (. - mean(.)) / sd(.)))
  x_df <- for_lasso %>% select(contains("_z"), -c(clinic_id_z, condition_mean_z))
  cv_model <- cv.glmnet(y = for_lasso$condition_mean, 
                        x = data.matrix(x_df), 
                        alpha = 1, lambda = l_seq_for_glm)
  coefs <- coef(cv_model, s = "lambda.min")

  lasso_results_inner <-
    data.frame(var = coefs@Dimnames[[1]][coefs@i + 1], 
               coef = coefs@x,
               rep = j) 
  
  cvm <- 
    data.frame(lambda = cv_model$lambda,
               mse = cv_model$cvm,
               rsq = 1 - (cv_model$cvm/var(for_lasso$condition_mean)),
               rep = j) %>%
    mutate(lambda_min = if_else(lambda == cv_model$lambda.min, 1, 0),
           lambda_1se = if_else(lambda == cv_model$lambda.1se, 1, 0))
  
  if (j == 1) {
    lasso_results_outer <- lasso_results_inner
    cvm_outer <- cvm
  } else {
    lasso_results_outer <- bind_rows(lasso_results_outer, lasso_results_inner)
    cvm_outer <- bind_rows(cvm_outer, cvm)
  }
}

mean_RMSE <- 
  cvm_outer %>%
  filter(lambda_min == 1) %>%
  summarise(mean_mse = mean(mse),
            sd_mse = sd(mse),
            mean_rsq = mean(rsq),
            sd_rsq = sd(rsq)) %>%
  mutate(mean_rmse = sqrt(mean_mse),
         sd_rmse = sqrt(sd_mse))

lasso_flex_new <- 
  lasso_results_outer %>%
  group_by(var) %>%
  summarise(coef_mean = mean(coef, na.rm = TRUE),
            coef_sd = sd(coef, na.rm = TRUE),
            coef_n = n(),
            coef_lower = coef_mean - 1.96*coef_sd/sqrt(coef_n),
            coef_upper = coef_mean + 1.96*coef_sd/sqrt(coef_n)) %>%
  # Require that coefficient is retained in > 50% of reps
  filter((coef_n >= n_reps*.5) & var != ("(Intercept)")) %>%
  mutate(`Mean (SD, n)` = paste0(round(coef_mean, 2), 
                                 " (", round(coef_sd, 2), ", ",coef_n, ")"),
         inclusion = "Posthoc 66pct:\nMean (SD, n)") %>%
  select(Variable = var, inclusion, `Mean (SD, n)`, coef_mean) %>%
  pivot_wider(names_from = inclusion, values_from = `Mean (SD, n)`) %>%
  arrange(desc(coef_mean)) %>%
  select(-coef_mean)

lasso_tmp <- lasso_flex_new
lasso_tmp[1,] <-
  as.list(
    c("mean_RMSE", 
      paste0(as.character(round(mean_RMSE$mean_rmse[1], 2)), " (", 
             as.character(round(mean_RMSE$sd_rmse[1], 2)), ", NA)")))
lasso_tmp[2,] <-
  as.list(
    c("mean_R2", 
      paste0(as.character(round(mean_RMSE$mean_rsq[1], 2)), " (", 
             as.character(round(mean_RMSE$sd_rsq[1], 2)), ", NA)")))

lasso_flex_new <-
  bind_rows(lasso_tmp[1:2,], lasso_flex_new) 

lasso_flex_new %>%
  kbl("html", digits = 3) %>%
  kable_classic(html_font = "Cambria", "basic", "center", full_width = F) %>%
  row_spec(1:2, background = "lightgray") %>%
  cat(., file = paste0(main_output_dir, "/posthoc_lasso_reg_table.html"))

