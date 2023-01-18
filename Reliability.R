# --------------------------
# Step 2: Reliability
# --------------------------

# ---------------------
# source functions & load libraries  
# ---------------------

#### ---- load packages and functions ---- #### 
if (!require("lme4")) install.packages("lme4"); require("lme4")
if (!require("lmerTest")) install.packages("lmerTest"); require("lmerTest")
if (!require("ggforce")) install.packages("ggforce"); require("ggforce")
if (!require("performance")) install.packages("performance"); require("performance")
if (!require("dplyr")) install.packages("dplyr"); require("dplyr")
if (!require("psych")) install.packages("psych"); require("psych")
if (!require("tidyverse")) install.packages("tidyverse"); require("tidyverse")
if (!require("patchwork")) install.packages("patchwork"); require("patchwork")
if (!require("flextable")) install.packages("flextable"); require("flextable")
if (!require("scales")) install.packages("scales"); require("scales")
if (!require("kableExtra")) install.packages("kableExtra"); require("kableExtra")
#detach plyr if loaded
if(any(grepl("package:plyr", search()))) detach("package:plyr")
source('https://www.dropbox.com/s/4mppmzaff55ra01/gradCPT_HelperFunctions.R?raw=1')
source('https://www.dropbox.com/s/bmhdm3d5xa7j7oa/json_splitter.R?raw=1')

current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path ))
frac <- .5
main_output_dir <- "20221013_50pct_inclusion"

# ---------------------
# load cognitive data
# ---------------------
#### ---- load cognitive data ---- ####
EXCLUDE_ONBOARDING = TRUE
REMOVE_FLAGGED = TRUE #do we want to remove EMAs flagged for non-compliance?

samples = c('glucog')
ema_data = list()
subject_data = list()
min_emas = list()
total_emas = list()

#glucog participants in progress
glucog_in_progress <- c()

##glucog
#row-per-ema data
ema_data$glucog <- read.csv("Files/ema_level_data.csv", stringsAsFactors = F)
subject_data$glucog <- read.csv("Files/subject_level_data.csv", stringsAsFactors = F)

#exclude onboarding?
if (EXCLUDE_ONBOARDING){
  ema_data$glucog = ema_data$glucog[ema_data$glucog$EMA_number > 0,]
}

# only include participants included in primary analyses
mydf_clean_20221014 <- read.csv(paste0("Files/", main_output_dir, "/mydf_clean_20221014.csv"))
id_list <- mydf_clean_20221014$user_id %>% unique(.)
ema_data$glucog = ema_data$glucog[ema_data$glucog$user_id %in% id_list,]
subject_data$glucog = subject_data$glucog[subject_data$glucog$user_id %in% id_list,]

enrolled_glucog = nrow(subject_data$glucog)

#rename "DSC" outcomes in subject_data "DSM" instead
dsc_rename = which(grepl('dsc_', names(subject_data$glucog)))
names(subject_data$glucog)[dsc_rename] = gsub('dsc_', replacement = 'dsm_', names(subject_data$glucog)[dsc_rename])

#row-per-subject data
#min emas completed for inclusion
min_emas$glucog = frac*45 #minimum amount of emas completed out of 45 total to be in analyses done here
total_emas$glucog = 45

##restrict subject-level data to participants who did minimum emas
for (s in 1:length(samples)){
  this_sample = samples[s]
  this_sub_data = subject_data[[this_sample]]
  subject_data[[this_sample]] = this_sub_data[this_sub_data$ema_batteries_completed >= min_emas[this_sample] & !is.na(this_sub_data$ema_batteries_completed),]
  
  this_ema_data = ema_data[[this_sample]]
  ema_data[[this_sample]] = this_ema_data[this_ema_data$user_id %in% subject_data[[this_sample]]$user_id,]
}

# -----------------
# Create half-test data (for WP reliability)
# -----------------

half_test_data = list()
for (s in 1:length(samples)){
  half_test_data[[samples[s]]] = list()
}

#### ---- MOT ---- ####
for (s in 1:length(samples)){
  this_sample = samples[s]
  this_ema_data = ema_data[[this_sample]]
  
  #get trial-level ema data, test trials only
  if (REMOVE_FLAGGED){
    mot_only = this_ema_data[!this_ema_data$mot_any_flag, c('session_id', 'mot_data')]
  } else {
    mot_only = this_ema_data[, c('session_id', 'mot_data')]
  }
  mot_only = mot_only[!is.na(mot_only$mot_data),]
  mot_long0 = deJSON.data(mot_only, 'mot_data')
  mot_long = mot_long0[mot_long0$type == 'test',]
  mot_long$hits = as.numeric(as.character(mot_long$hits))
  
  #add trial number column
  mot_long = mot_long %>% 
    group_by(session_id) %>% 
    mutate(trial_number = row_number())
  
  #assign trials 1,4,6 to half 1, and trials 2,3,5 to half 2 (did this to approximately match % correct on each half)
  mot_long$half = NA
  mot_long$half[mot_long$trial_number %in% c(1,4,6)] = 1
  mot_long$half[mot_long$trial_number %in% c(2,3,5)] = 2
  
  #get mean score by half
  mot_half_score  = mot_long %>%
    group_by(session_id, half) %>%
    summarize(half_score = sum(hits)) %>%
    select(session_id, half, half_score)
  half_test_data[[samples[s]]]$mot = merge(mot_half_score, this_ema_data[!is.na(this_ema_data$mot_data),], by='session_id')
}

#### ---- DSM ---- ####
for (s in 1:length(samples)){
  this_sample = samples[s]
  this_ema_data = ema_data[[this_sample]]
  
  #get trial-level ema data, test trials only
  
  if (REMOVE_FLAGGED){
    dsm_only = this_ema_data[!this_ema_data$dsm_any_flag, c('session_id', 'dsm_data')]
  } else {
    dsm_only = this_ema_data[, c('session_id', 'dsm_data')]
  }
  
  dsm_only = dsm_only[!is.na(dsm_only$dsm_data) & nchar(dsm_only$dsm_data) > 2,]
  dsm_long0 = deJSON.data(dsm_only, 'dsm_data')
  dsm_long = dsm_long0[dsm_long0$type == 'test',]
  dsm_long$correct = as.numeric(as.character(dsm_long$correct))
  dsm_long$rt = as.numeric(as.character(dsm_long$rt))
  
  ## medianRTc
  #limit to correct trials only
  dsm_long_correct = dsm_long[dsm_long$correct == 1,]
  #add trial number column
  dsm_long_correct <- dsm_long_correct %>% 
    group_by(session_id) %>% 
    mutate(trial_number = row_number())
  #split into even and odd trials (based on correct trials only)
  dsm_long_correct$half = NA
  dsm_long_correct$half[dsm_long_correct$trial_number %% 2 == 1] = 1
  dsm_long_correct$half[dsm_long_correct$trial_number %% 2 == 0] = 2
  
  #get mean score by half
  dsm_half_score_correct  = dsm_long_correct %>%
    group_by(session_id, half) %>%
    summarize(half_score_mRTC = median(rt)) %>%
    select(session_id, half, half_score_mRTC)
  
  # make column for trial start time and delay time
  dsm_half_by_time <- dsm_long
  dsm_half_by_time$delay_time = NA #time between pressing down on key and next trial starting
  
  #change these variables to numeric
  dsm_half_by_time$rt = as.numeric(dsm_half_by_time$rt)
  dsm_half_by_time$timestamp = as.numeric(dsm_half_by_time$timestamp)
  
  #first test trial starts at 0
  dsm_half_by_time <-
    dsm_half_by_time %>%
    group_by(session_id) %>%
    mutate(trial_start_time = if_else(timestamp == min(timestamp), 0, NA_real_)) %>%
    nest()
  
  #probably a more efficient way to do this that could be done in dplyr per ema
  for (j in 1:nrow(dsm_half_by_time)) {
    df <- dsm_half_by_time$data[[j]]
    for (t in 2:nrow(df)){
      df$delay_time[t-1] = (df$timestamp[t] - df$timestamp[t-1]) - df$rt[t]
      df$trial_start_time[t] =  df$trial_start_time[t-1] + df$rt[t-1] + df$delay_time[t-1]
    }
    # Splitting into quadrants & testing 1, 4 vs. 2, 3 to balance practice effects
    df$half_time <- ifelse((df$trial_start_time < max(df$trial_start_time)*.25) |
                             ((df$trial_start_time >= max(df$trial_start_time)*.5) &
                                (df$trial_start_time < max(df$trial_start_time)*.75)), 1, 2)
    dsm_half_by_time$data[[j]] <- df
  }
  
  dsm_half_by_time <-
    dsm_half_by_time %>%
    unnest(cols = c(data)) %>%
    ungroup() %>%
    group_by(session_id, half_time) %>%
    summarise(half_score_correct = sum(correct)) %>%
    select(session_id, half_time, half_score_correct) %>%
    rename(half = half_time)
  
  dsm_half_score_correct <- left_join(dsm_half_score_correct, dsm_half_by_time)
  
  half_test_data[[samples[s]]]$dsm = merge(dsm_half_score_correct, this_ema_data[!is.na(dsm_only$dsm_data) & nchar(dsm_only$dsm_data) > 2,], by='session_id')
}

#### ---- gradCPT ---- ####
for (s in 1:length(samples)){
  this_sample = samples[s]
  this_ema_data = ema_data[[this_sample]]
  
  #60 go, 15 no go per test...we will split by go/nogo even/odd (not overall trial number)
  #get trial-level ema data, test trials only
  
  if (REMOVE_FLAGGED){
    gradcpt_only = this_ema_data[!this_ema_data$gradcpt_any_flag, c('session_id', 'gradcpt_data')]
  } else {
    gradcpt_only = this_ema_data[, c('session_id', 'gradcpt_data')]
  }
  
  gradcpt_only = gradcpt_only[!is.na(gradcpt_only$gradcpt_data) & nchar(gradcpt_only$gradcpt_data) > 2,]
  gradcpt_long0 = deJSON.data(gradcpt_only, 'gradcpt_data')
  gradcpt_long = gradcpt_long0[gradcpt_long0$type == 'test',]
  gradcpt_long$correct = as.numeric(as.character(gradcpt_long$correct))
  gradcpt_long$accuracy = as.numeric(as.character(gradcpt_long$accuracy))
  gradcpt_long$accuracy[gradcpt_long$accuracy == -1] = 0
  gradcpt_long$rt = as.numeric(as.character(gradcpt_long$rt))
  gradcpt_long$trial = as.numeric(as.character(gradcpt_long$trial))
  names(gradcpt_long)[names(gradcpt_long) == 'correct'] = 'go'
  
  #add trial number column by go/nogo
  gradcpt_long <- gradcpt_long %>% 
    group_by(session_id, go) %>% 
    mutate(trial_number = row_number())
  
  #split into even and odd trials
  gradcpt_long$half = NA
  gradcpt_long$half[gradcpt_long$trial_number %% 2 == 1] = 1
  gradcpt_long$half[gradcpt_long$trial_number %% 2 == 0] = 2
  
  #get info by half
  gradcpt_half_info = gradcpt_long %>%
    group_by(session_id, half, go) %>%
    summarize(correct = sum(accuracy), count = n()) %>%
    ungroup() %>%
    complete(session_id, half, go,
             fill = list(correct = 0, freq = 0)) %>%
    select(session_id, half, go, correct, count)
  
  #adjust 100% and 0% for dprime calculation
  gradcpt_half_info_adjusted = gradcpt_half_info
  #100% hit rate
  gradcpt_half_info_adjusted$correct[gradcpt_half_info_adjusted$go == 1 & gradcpt_half_info_adjusted$correct == gradcpt_half_info_adjusted$count] = 
    gradcpt_half_info_adjusted$count[gradcpt_half_info_adjusted$go == 1 & gradcpt_half_info_adjusted$correct == gradcpt_half_info_adjusted$count] - .5
  #0% hit rate
  gradcpt_half_info_adjusted$correct[gradcpt_half_info_adjusted$go == 1 & gradcpt_half_info_adjusted$correct == 0] = .5
  #100% false alarm rate
  gradcpt_half_info_adjusted$correct[gradcpt_half_info_adjusted$go == 0 & gradcpt_half_info_adjusted$correct == 0] = .5
  #0% false alarm rate
  gradcpt_half_info_adjusted$correct[gradcpt_half_info_adjusted$go == 0 & gradcpt_half_info_adjusted$correct == gradcpt_half_info_adjusted$count] = 
    gradcpt_half_info_adjusted$count[gradcpt_half_info_adjusted$go == 0 & gradcpt_half_info_adjusted$correct == gradcpt_half_info_adjusted$count] - .5
  gradcpt_half_info_adjusted$perc = NA
  gradcpt_half_info_adjusted$perc[gradcpt_half_info_adjusted$go==1] = gradcpt_half_info_adjusted$correct[gradcpt_half_info_adjusted$go==1]/
    gradcpt_half_info_adjusted$count[gradcpt_half_info_adjusted$go==1]
  gradcpt_half_info_adjusted$perc[gradcpt_half_info_adjusted$go==0] = 
    (gradcpt_half_info_adjusted$count[gradcpt_half_info_adjusted$go==0] - gradcpt_half_info_adjusted$correct[gradcpt_half_info_adjusted$go==0])/
    gradcpt_half_info_adjusted$count[gradcpt_half_info_adjusted$go==0]
  gradcpt_odd_go = gradcpt_half_info_adjusted[gradcpt_half_info_adjusted$go == 1 & gradcpt_half_info_adjusted$half == 1, c('session_id', 'perc')]
  gradcpt_odd_nogo = gradcpt_half_info_adjusted[gradcpt_half_info_adjusted$go == 0 & gradcpt_half_info_adjusted$half == 1, c('session_id', 'perc')]
  gradcpt_odd = merge(gradcpt_odd_go, gradcpt_odd_nogo, by = 'session_id')
  gradcpt_odd$half = 1
  gradcpt_even_go = gradcpt_half_info_adjusted[gradcpt_half_info_adjusted$go == 1 & gradcpt_half_info_adjusted$half == 2, c('session_id', 'perc')]
  gradcpt_even_nogo = gradcpt_half_info_adjusted[gradcpt_half_info_adjusted$go == 0 & gradcpt_half_info_adjusted$half == 2, c('session_id', 'perc')]
  gradcpt_even = merge(gradcpt_even_go, gradcpt_even_nogo, by = 'session_id')
  gradcpt_even$half = 2
  gradcpt_half_score = rbind(gradcpt_odd, gradcpt_even)
  names(gradcpt_half_score) = c('session_id', 'HR', 'FAR', 'half')
  gradcpt_half_level = merge(gradcpt_half_score, this_ema_data[!is.na(gradcpt_only$gradcpt_data) & nchar(gradcpt_only$gradcpt_data) > 2,], by='session_id')
  
  # dprime calculation by half
  d_and_c = dprime(gradcpt_half_level$HR, gradcpt_half_level$FAR)
  gradcpt_half_level = cbind(gradcpt_half_level, d_and_c)
  names(gradcpt_half_level)[names(gradcpt_half_level) == 'bias'] = 'crit'
  names(gradcpt_half_level)[names(gradcpt_half_level) == 'd'] = 'dprime'
  
  ## medianRT calculation by half
  medRT <-
    gradcpt_long %>%
    select(session_id, half, go, accuracy, rt) %>%
    # median RT on correct trials
    # need to filter for trials where they make a response (i.e., GO trials)
    filter(accuracy == 1 & go == 1) %>%
    group_by(session_id, half) %>%
    summarize(gradhalf_score_mRTC = median(rt)) %>%
    select(session_id, half, gradhalf_score_mRTC) %>%
    ungroup()
  gradcpt_half_level <- left_join(gradcpt_half_level, medRT)
  
  # combine together
  half_test_data[[samples[s]]]$gradcpt = gradcpt_half_level
}

# -----------------------
# Calculate M, SD, & BP/WP reliability
# -----------------------
ema_summary_data = data.frame(matrix(nrow = 5, ncol = 6))
names(ema_summary_data) = c('Outcome', 'Mean (SD)', 'Range', 'WPR', 'BPR', 'N (T)')

test_names = c("mot","dsm", "dsm", "gradcpt", "gradcpt")
ema_outcome_name = c('mot_correct', 'dsm_medianRTc', 'dsm_num_correct', 'gradcpt_medianRTc', 'gradcpt_dprime')
half_level_score = c('half_score', 'half_score_mRTC', 'half_score_correct', 'gradhalf_score_mRTC', 'dprime')
outcome_names = c('MOT Accuracy', 'DSM medRT', 'DSM num correct', 'GradCPT medRT', 'GradCPT d-prime')
if (REMOVE_FLAGGED){
  ema_mean_names = c('mot_correct_ema_mean_postQC', 'dsm_medianRTc_ema_mean_postQC', 
                     'dsm_num_correct_ema_mean_postQC', 'gradcpt_medianRTc_ema_mean_postQC',
                     'gradcpt_dprime_ema_mean_postQC')
  ema_flag_names = c('mot_any_flag', 'dsm_any_flag', 
                     'dsm_any_flag', 'gradcpt_any_flag', 'gradcpt_any_flag')
} else {
  ema_mean_names = c('mot_correct_ema_mean', 'dsm_medianRTc_ema_mean', 
                     'dsm_num_correct_ema_mean', 'gradcpt_medianRTc_ema_mean', 'gradcpt_dprime_ema_mean')
  ema_flag_names = c('mot_any_flag', 'dsm_any_flag', 
                     'dsm_any_flag', 'gradcpt_any_flag', 'gradcpt_any_flag')
}
#scaling and rounding outcomes for table
scale_factor = c(100, 1, 2, 1, 1)
round_factor = c(1, 0, 1, 0, 2)

for (t in 1:length(test_names)){
  this_test = test_names[t]
  this_half_score = half_level_score[t]
  ema_summary_data$Outcome[t] = outcome_names[t]
  for (s in 1:length(samples)){
    this_sample = samples[s]
    if (!(this_test == 'crt' & this_sample == 'glucog')){
      
      x = scale_factor[t]
      r = round_factor[t]
      
      #ema mean, sd, min, max, N
      # mean
      this_data = x*subject_data[[this_sample]][[ema_mean_names[t]]]
      this_df = subject_data[[this_sample]]
      sdata4comparison = this_df[c('user_id', 'demo_age', 'demo_gender', ema_mean_names[t])]
      
      # mean, sd, min, max (summarizing participant-avg. stats)
      m_mean = round(mean(this_data, na.rm = T),r)
      sd_mean = round(sd(this_data, na.rm = T),r)
      min_mean = round(min(this_data, na.rm = T),r)
      max_mean = round(max(this_data, na.rm = T),r)
      N <- length(this_data)
      obs <- 
        ema_data$glucog %>% 
        select(user_id, all_of(ema_flag_names[t])) %>%
        filter_at(2, all_vars(. == 0)) %>% 
        group_by(user_id) %>% 
        summarise(count = n()) %>% 
        ungroup() %>% 
        summarise(mean = mean(count)) %>% 
        pull(mean) %>% round(., 1)
      
      msd_mean = paste0(m_mean, ' (', sd_mean, ')')
      c1 = which(str_detect(names(ema_summary_data), "Mean"))
      ema_summary_data[t, c1] = msd_mean
      
      range_mean = paste0(min_mean, ' - ', max_mean)
      c2 = which(str_detect(names(ema_summary_data), "Range"))
      ema_summary_data[t, c2] = range_mean
      
      n_t = paste0(N, ' (', obs, ')')
      c3 = which(str_detect(names(ema_summary_data), "N"))
      ema_summary_data[t, c3] = n_t
      
      #reliability
      mlr_output =  mlr(half_test_data[[this_sample]][[this_test]], 
                        grp = 'user_id', Time = 'EMA_number', items = 'half', 
                        values = this_half_score, lmer = TRUE, lme = FALSE, long = TRUE)
      c5 = str_detect(names(ema_summary_data), "BPR")
      ema_summary_data[t, c5] = round(mlr_output$RkRn,2)
      c6 = str_detect(names(ema_summary_data), "WPR")
      ema_summary_data[t, c6] = round(mlr_output$Rcn,2)
    }
  }
}

# -------------------
# Save to file
# -------------------

ema_summary_data %>%
  kbl("html") %>%
  kable_classic(html_font = "Cambria", "basic", "center", full_width = F) %>%
  cat(., file = "Files/Output/reliability_50pct.html")






