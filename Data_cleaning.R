# --------------------------
# Step 1: Data cleaning
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
frac <- .5
main_output_dir <- "20221013_50pct_inclusion"

# ---------------------
# Read in & clean data
# ---------------------

# Gl variables & EMA-level data
mydf <- read.csv("Files/overlapWindows_gl_and_EMA_2022_08_05.csv")

# identify subjects who were withdrawn from analysis
glu_withdrawn <- c(
  #never started at clinic or reassigned user_id number
  "GLUCOG042016", "GLUCOG045011", "GLUCOG045027",
  # technical reasons (before reaching 50%)
  "GLUCOG005004",
  # involuntary (before reaching 50%)
  "GLUCOG045002", "GLUCOG045045", "GLUCOG005007", "GLUCOG005005", 
  "GLUCOG045029", "GLUCOG045013","GLUCOG005076", "GLUCOG042003",
  # voluntary (before reaching 50%)
  "GLUCOG035034", "GLUCOG035054")

# Clean EMA using flags defined in glucog_pipeline_functions.R
TMB_dat <-
  read.csv("Files/ema_level_data.csv") %>%
  filter(!(user_id %in% glu_withdrawn)) %>%
  filter(!isOnboarding) %>%
  mutate(time = ymd_hms(test_start_time_clinic)) %>%
  # participant has repeated timestamps (two devices)
  filter(user_id != "GLUCOG005024") %>%
  # fix EMA count issue
  mutate(EMA_count_orig = EMA_count) %>%
  group_by(user_id) %>%
  arrange(user_id, time) %>%
  # fix EMA count coding error (EMA count = number of completed EMAs)
  mutate(EMA_count = row_number(),
         time = ymd_hms(time)) %>%
  ungroup() %>%
  {. ->> TMB_tmp} %>%
  select(user_id, sitting_id, EMA_count_orig,
         dsm_accuracy_flag, dsm_ncorrect_flag, gradcpt_omission_flag, mot_framerate_flag,
         interruptions_anyInterruptions,
         gradcpt_commission_errors, gradcpt_sdRTc, gradcpt_cvRTc,
         survey_stress_score, survey_depression_emotion_score, survey_anxiety_emotion_score, survey_context_score,
         survey_glucose_estimate_numeric, survey_sleep_quality, survey_sleep_time) %>%
  gather(key = key, value = value, contains("flag")) %>%
  # clean cog data as described below (refer to functions in glucog_pipeline_functions.R)
  # exclude DS accuracy < 50% & n correct < 6
  # exclude GCPT ommision errors > 50%
  # exclude MOT in power save (frame time > 30 ms)
  filter(value == 0 & !is.na(value)) %>%
  separate(key, into = c("test", "omit"), sep = "_", extra = "merge") %>%
  select(-omit, -value) %>%
  mutate(new = paste0(user_id, "_", EMA_count_orig, "_", test)) %>%
  mutate(interruptions_anyInterruptions = if_else(interruptions_anyInterruptions == "yes", 1, 0),
         survey_sleep_quality_tmp = if_else(survey_sleep_quality == "Very poorly", 0, 999),
         survey_sleep_quality_tmp = if_else(survey_sleep_quality == "Somewhat poorly", 1, survey_sleep_quality_tmp),
         survey_sleep_quality_tmp = if_else(survey_sleep_quality == "Somewhat well", 2, survey_sleep_quality_tmp),
         survey_sleep_quality_tmp = if_else(survey_sleep_quality == "Very well", 3, survey_sleep_quality_tmp),
         survey_sleep_quality = survey_sleep_quality_tmp) %>%
  select(-survey_sleep_quality_tmp)

# Write EMA data to file
write.csv(TMB_dat,
          paste0("Files/", main_output_dir, "/TMB_dat_20221014.csv"),
          row.names = FALSE)
write.csv(TMB_tmp,
          paste0("Files/", main_output_dir, "/TMB_tmp_20221014.csv"),
          row.names = FALSE)

# Merge & clean EMA + CGM data
mydf_clean <-
  mydf %>%
  mutate(EMA_count_orig = EMA_count) %>%
  group_by(user_id) %>%
  arrange(user_id, time) %>%
  # fix EMA count coding error (EMA count = number of completed EMAs)
  mutate(EMA_count = row_number(),
         time = ymd_hms(time)) %>%
  ungroup() %>%
  left_join(TMB_dat %>% unique(.)) %>%
  filter(!is.na(new)) %>%
  # require touch observations
  filter(touch == "true" & finished == 1) %>%
  # convert wide to long
  gather(key = key, value = value,
         c((starts_with("gradcpt") | starts_with("dsm") | starts_with("mot")) & !contains("position"))) %>%
  separate(key, into = c("test", "key"), sep = "_", extra = "merge") %>%
  mutate(battery_position = if_else(test == "dsm", dsm_battery_position, as.integer(999)),
         battery_position = if_else(test == "mot", mot_battery_position, battery_position),
         battery_position = if_else(test == "gradcpt", gradcpt_battery_position, battery_position)) %>%
  select(-dsm_battery_position, -mot_battery_position,
         -gradcpt_battery_position, -new, -EMA_count_orig) %>%
  # Call unique to eliminate rows that were duplicated during left_join with "new"
  unique(.) %>%
  # Compute total number of EMAs completed after cleaning
  group_by(user_id, EMA_count) %>%
  nest() %>% ungroup() %>%
  arrange(user_id, EMA_count) %>%
  group_by(user_id) %>%
  mutate(EMA_count_clean = row_number()) %>%
  ungroup() %>%
  unnest(data) %>%
  # Require `frac` completion rate AFTER data cleaning
  group_by(user_id) %>%
  nest() %>% ungroup() %>%
  mutate(max_clean = map_dbl(.x = data, ~(.x %>% filter(EMA_count_clean == max(EMA_count_clean)) %>%
                                            pull(EMA_count_clean) %>% unique(.))),
         max = map_dbl(.x = data, ~(.x %>% filter(EMA_count == max(EMA_count)) %>% pull(EMA_count) %>% unique(.))),
         EMA_count_binary = if_else((max >= (frac*45)) & (max <= 45), 1, 0),
         EMA_count_clean_binary = if_else((max_clean >= (frac*45)) & (max <= 45), 1, 0)) %>%
  filter(EMA_count_clean_binary == 1) %>%
  select(-contains("max"), -contains("binary")) %>%
  unnest(data) %>%
  # Fill in sleep time & sleep quality data
  mutate(day = date(time)) %>%
  group_by(user_id, day) %>%
  nest() %>%
  ungroup() %>%
  mutate(sleep_quality_temp = map(.x = data, ~.x %>% pull(survey_sleep_quality) %>% na.omit(.) %>% unique(.)),
         data = map2(.x = data, .y = sleep_quality_temp,
                     ~.x %>% mutate(survey_sleep_quality = if_else(length(.y) == 0, NA_real_, .y[1]))),
         sleep_time_temp = map(.x = data, ~.x %>% pull(survey_sleep_time) %>% na.omit(.) %>% unique(.)),
         data = map2(.x = data, .y = sleep_time_temp,
                     ~.x %>% mutate(survey_sleep_time = if_else(length(.y) == 0, NA_real_, .y[1])))) %>%
  select(-ends_with("temp"), -day) %>%
  unnest(cols = c(data))

# Write to file
write.csv(mydf_clean,
          paste0("Files/", main_output_dir, "/mydf_clean_20221014.csv"),
          row.names = FALSE)

# ---------------------
# Do participants included in analysis differ from the full sample?
# ---------------------

mydf_clean <- read.csv(paste0("Files/", main_output_dir, "/mydf_clean_20221014.csv"))

full_sample <-
  read.csv("Files/ema_level_data.csv") %>%
  filter(!(user_id %in% glu_withdrawn)) %>%
  filter(!isOnboarding) %>%
  mutate(time = ymd_hms(test_start_time_clinic),
         included = if_else(user_id %in% mydf_clean$user_id, 1, 0),
         sample = "Full") 

for_analysis <- 
  full_sample %>% 
  filter(included == 1) %>%
  mutate(sample = "Included") %>%
  bind_rows(full_sample) %>%
  select(user_id, sample, demo_age, demo_gender, demo_education, 
         demo_race = demo_ethnicity, demo_ethnicity = demo_hispanic) %>%
  # simplify coding for chisq analysis - otherwise, too few levels for comparison
  mutate(demo_gender_male = if_else(demo_gender == "male", 1, 0),
         demo_gender_female = if_else(demo_gender == "female", 1, 0),
         demo_education_num = if_else(demo_education == "high", 0, NA_real_),
         demo_education_num = if_else(demo_education == "some_college", 1, demo_education_num),
         demo_education_num = if_else(demo_education == "technical", 1, demo_education_num),
         demo_education_num = if_else(demo_education == "college", 2, demo_education_num),
         demo_education_num = if_else(demo_education == "masters", 3, demo_education_num),
         demo_education_num = if_else(demo_education == "grad", 4, demo_education_num),
         demo_ethnicity_africanOrBlack = if_else(demo_race == "africanOrBlack", 1, 0),
         demo_ethnicity_europeanOrWhite = if_else(demo_race == "europeanOrWhite", 1, 0),
         demo_hispanic_yes = if_else(demo_ethnicity == "yes", 1, 0)) %>%
  select(-demo_education, -demo_ethnicity, -demo_race, -demo_gender) %>%
  distinct(.) %>%
  gather(key = key, value = value, -user_id, -sample) %>%
  group_by(key) %>% nest() %>%
  ungroup()

# numeric variables
for_analysis %>%
  filter(key == "demo_age" | key == "demo_education_num") %>%
  mutate(data = map(.x = data, ~.x %>% 
                      mutate(value = as.numeric(value),
                             sample = as.factor(sample))),
         lm = map(.x = data, ~lm(value ~ sample, data = .x)),
         tidied = map(lm, tidy)) %>%
  unnest(tidied) %>%
  filter(term != "(Intercept)")
  
# factor variables
for_analysis %>%
  filter(key != "demo_age" & key != "demo_education_num") %>%
  mutate(data = map(.x = data, 
                    ~.x %>% 
                      mutate(value = as.factor(value),
                             sample = as.factor(sample))),
         chisq = map(.x = data, ~chisq_test(.x, value ~ sample))) %>%
  unnest(chisq)

# ---------------------
# When were CGM data from second device scrubbed? 
# ---------------------
# script to generate CSV:
# ~/Dropbox (Partners HealthCare)/BachTech_Mats/GluCogAim1/R_Scripts/DerivingCGMMetrics_overlappingWindows.R
mydf_clean <- read.csv(paste0("Files/", main_output_dir, "/mydf_clean_20221014.csv"))
CGM_scrubbed <- read.csv("Files/CGM_scrubbed.csv") %>%
  group_by(user_id) %>%
  filter(time == max(time, na.rm = TRUE)) %>%
  rename(CGM_scrubbed = time)

mydf_clean %>%
  mutate(time = as.Date(time)) %>%
  select(user_id, time, EMA_count) %>%
  filter(EMA_count == 1) %>%
  distinct(.) %>%
  left_join(CGM_scrubbed) %>%
  mutate(dif = difftime(CGM_scrubbed, time, units = 'days')) %>%
  filter(dif > 2) %>%
  summarise(mean = mean(dif))





