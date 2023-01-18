# --------------------------
# Step 3: HB Modeling
# --------------------------

# From: R_Scripts/FinalDF_hypothesisDriven_20221014.R

# ---------------------
# source functions & load libraries  
# ---------------------
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path ))
source("Functions.R")
load_packs()

# ---------------------
# Set script parameters
# ---------------------
run_models <- FALSE
frac <- .5
main_output_dir <- "20221013_50pct_inclusion"

# ---------------------
# Read in data & prep for modeling
# ---------------------

mydf_clean <- read.csv(paste0("Files/", main_output_dir, "/mydf_clean_20221014.csv"))

to_model <-
  mydf_clean %>%
  select(user_id, test, outcome = key, value, MinLag_0) %>% 
  # calculate deviations from person-level means
  group_by(user_id) %>%
  mutate_if(is.numeric, 
            .funs = list(WP = ~ ((.x - mean(.x, na.rm = TRUE))/sd(.x, na.rm = TRUE)))) %>%
  ungroup() %>%
  # getting rid of duplicated rows
  distinct(.) %>%
  # nest prior to running models
  group_by(test, outcome) %>%
  nest() %>% ungroup()

# -----------------
# HBM to estimate participant-level slopes
# -----------------

setwd(paste0("Files/", main_output_dir))

#-----------------------
# DSM
#-----------------------

if (run_models) {
  to_model_dsm <- 
    to_model %>%
    unnest(data) %>%
    filter(test == "dsm" & !is.na(MinLag_0_WP)) %>%
    mutate(MinLag_0_WP_sq = MinLag_0_WP^2) %>%
    filter(
      (MinLag_0_WP_sq > mean(MinLag_0_WP_sq, na.rm = TRUE) - 3*sd(MinLag_0_WP_sq, na.rm = TRUE)) &
        (MinLag_0_WP_sq < mean(MinLag_0_WP_sq, na.rm = TRUE) + 3*sd(MinLag_0_WP_sq, na.rm = TRUE))
    ) %>%
    group_by(test, outcome) %>%
    nest() %>% ungroup() 
  write_rds(to_model_dsm, "to_model_dsm.rds")
  
  options(mc.cores = parallel::detectCores())
  to_model_dsm_row1 <-
    stan_lmer(value ~ poly(MinLag_0_WP, 2, raw = FALSE) +
                (poly(MinLag_0_WP, 2, raw = FALSE) | user_id),
              data = to_model_dsm$data[[1]],
              control=list(adapt_delta=0.999),
              iter = 20000, chains = 4)
  write_rds(to_model_dsm_row1, "to_model_dsm_row1.rds")
  
  to_model_dsm_row2 <-
    stan_lmer(value ~ poly(MinLag_0_WP, 2, raw = FALSE) +
                (poly(MinLag_0_WP, 2, raw = FALSE) | user_id),
              data = to_model_dsm$data[[2]],
              control=list(adapt_delta=0.999),
              iter = 15000, chains = 4)
  write_rds(to_model_dsm_row2, "to_model_dsm_row2.rds")
} else {
  to_model_dsm_row1 <- readRDS("to_model_dsm_row1.rds")
  to_model_dsm_row2 <- readRDS("to_model_dsm_row2.rds")
}

# Confirmed valid estimates for all parameters
# cf. https://mc-stan.org/rstan/reference/monitor.html
# Bulk & tail ESS > 100 * n_chains
launch_shinystan(to_model_dsm_row1, ppd = FALSE)
launch_shinystan(to_model_dsm_row2, ppd = FALSE)

row1 <- monitor(to_model_dsm_row1$stanfit)
row2 <- monitor(to_model_dsm_row2$stanfit)

row1 %>% as.data.frame() %>% filter((valid != 1) |  (Bulk_ESS == min(Bulk_ESS)) | (Tail_ESS == min(Tail_ESS)))
row2 %>% as.data.frame() %>% filter((valid != 1) |  (Bulk_ESS == min(Bulk_ESS)) | (Tail_ESS == min(Tail_ESS)))

#-----------------------
# GCPT
#-----------------------

if (run_models) {
  to_model_gcpt <-
    to_model %>%
    unnest(data) %>%
    filter(test == "gradcpt" & !is.na(MinLag_0_WP)) %>%
    mutate(MinLag_0_WP_sq = MinLag_0_WP^2) %>%
    filter(
      (MinLag_0_WP_sq > mean(MinLag_0_WP_sq, na.rm = TRUE) - 3*sd(MinLag_0_WP_sq, na.rm = TRUE)) &
        (MinLag_0_WP_sq < mean(MinLag_0_WP_sq, na.rm = TRUE) + 3*sd(MinLag_0_WP_sq, na.rm = TRUE))
    ) %>%
    group_by(test, outcome) %>%
    nest() %>% ungroup()
  write_rds(to_model_gcpt, "to_model_gcpt.rds")
  
  options(mc.cores = parallel::detectCores())
  to_model_gcpt_row1 <-
    stan_lmer(value ~ poly(MinLag_0_WP, 2, raw = FALSE) +
                (poly(MinLag_0_WP, 2, raw = FALSE) | user_id),
              data = to_model_gcpt$data[[1]],
              control=list(adapt_delta=0.999),
              iter = 15000, chains = 4)
  write_rds(to_model_gcpt_row1, "to_model_gcpt_row1.rds")
  
  to_model_gcpt_row2 <-
    stan_lmer(value ~ poly(MinLag_0_WP, 2, raw = FALSE) +
                (poly(MinLag_0_WP, 2, raw = FALSE) | user_id),
              data = to_model_gcpt$data[[2]],
              control=list(adapt_delta=0.999),
              iter = 8000, chains = 4)
  write_rds(to_model_gcpt_row2, "to_model_gcpt_row2.rds")
} else {
  to_model_gcpt_row1 <- readRDS("to_model_gcpt_row1.rds")
  to_model_gcpt_row2 <- readRDS("to_model_gcpt_row2.rds")
}

launch_shinystan(to_model_gcpt_row1, ppd = FALSE)
launch_shinystan(to_model_gcpt_row2, ppd = FALSE)

row1 <- monitor(to_model_gcpt_row1$stanfit)
row2 <- monitor(to_model_gcpt_row2$stanfit)

row1 %>% as.data.frame() %>% filter((valid != 1) |  (Bulk_ESS == min(Bulk_ESS)) | (Tail_ESS == min(Tail_ESS)))
row2 %>% as.data.frame() %>% filter((valid != 1) |  (Bulk_ESS == min(Bulk_ESS)) | (Tail_ESS == min(Tail_ESS)))

#-----------------------
# DSM, 2022-11-10 post-hoc: with centered (but not scaled) IVs
# Rerunning DSM medRTc, in keeping with manuscript's focus
#-----------------------
 if (frac == .66) {
   if (run_models) {
     to_model_dsm_OC <- 
       to_model %>%
       unnest(data) %>%
       filter(test == "dsm" & !is.na(MinLag_0_WP)) %>%
       mutate(MinLag_0_WP_sq = MinLag_0_WP^2) %>%
       filter(
         (MinLag_0_WP_sq > mean(MinLag_0_WP_sq, na.rm = TRUE) - 3*sd(MinLag_0_WP_sq, na.rm = TRUE)) &
           (MinLag_0_WP_sq < mean(MinLag_0_WP_sq, na.rm = TRUE) + 3*sd(MinLag_0_WP_sq, na.rm = TRUE))
       ) %>%
       group_by(test, outcome, user_id) %>%
       mutate(MinLag_0_centered = MinLag_0 - mean(MinLag_0, na.rm = TRUE)) %>%
       ungroup() %>%
       group_by(test, outcome) %>%
       nest() %>% ungroup() 
     write_rds(to_model_dsm_OC, "testCentering_rawFalse/to_model_dsm_OC.rds")
     
     options(mc.cores = parallel::detectCores())
     to_model_dsm_row2_centered <-
       stan_lmer(value ~ poly(MinLag_0_centered, 2, raw = FALSE) +
                   (poly(MinLag_0_centered, 2, raw = FALSE) | user_id),
                 data = to_model_dsm_OC$data[[2]],
                 control=list(adapt_delta=0.99),
                 iter = 20000, chains = 4)
     write_rds(to_model_dsm_row2_centered, "testCentering_rawFalse/to_model_dsm_row2_cIV.rds")
   } else {
     to_model_dsm_row2_cIV <- readRDS("testCentering_rawFalse/to_model_dsm_row2_cIV.rds")
   }
   
   # Bulk & tail ESS > 100 * n_chains
   launch_shinystan(to_model_dsm_row2_cIV, ppd = FALSE)
   row1_cIV <- monitor(to_model_dsm_row2_cIV$stanfit)
   row1_cIV %>% as.data.frame() %>% filter((valid != 1) |  (Bulk_ESS == min(Bulk_ESS)) | (Tail_ESS == min(Tail_ESS)))
 }


