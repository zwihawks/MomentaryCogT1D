# --------------------------
# Functions
# --------------------------

load_packs <- function() {
  library(corrplot)
  library(corrr)
  library(kableExtra)
  library(magick)
  library(infer)
  library(magrittr)
  library(purrr)
  library(forcats)
  library(tidyr)
  library(modelr)
  library(ggdist)
  library(tidybayes)
  library(ggplot2)
  library(cowplot)
  library(rstan)
  library(rstanarm)
  library(RColorBrewer)
  library(corrr)
  library(corrplot)
  library(ggdag)
  library(gghighlight)
  library(DescTools)
  library(ggrepel)
  library(tidyverse)
  library(lubridate)
  library(ggnewscale)
  library(ggpubr)
  library(glmnet)
  library(flextable)
  library(matrixcalc)
  library(igraph)
  library(infotheo)
  library(psych)
  library(egg)
  library(matrixStats)
  library(lme4)
  library(broom)
  library(sjPlot)
  library(MuMIn)
  library(tidymodels)
  library(parallel)
  library(fs)
  library(dplyr)
}

wrangle_subject_file <- function(df, select_SB = TRUE) {
  subject_data <- 
    df %>% 
    # remove first row (all NA)
    filter(user_id != "") %>%
    # recode severe hypo & DKA events as ordinal
    # For 12-mo. hx, limited categories to 0 (coded as 0), 1 (coded as 1), and more than 1 (coded as 2))
    dplyr::mutate(clinic__SHNumEverB = if_else(str_detect(SHNumEverB, ">"), "6", SHNumEverB)) %>%
    dplyr::mutate(clinic__SHNumEverB = as.numeric(if_else(str_detect(clinic__SHNumEverB, "-"), "5", clinic__SHNumEverB))) %>%
    dplyr::mutate(clinic__SHLast12MonthsB = if_else(SHLast12MonthsB == "", "0", SHLast12MonthsB)) %>%
    dplyr::mutate(clinic__SHLast12MonthsB = as.numeric(if_else(str_detect(clinic__SHLast12MonthsB, "-") |
                                                                 str_detect(clinic__SHLast12MonthsB, "2") |
                                                                 str_detect(clinic__SHLast12MonthsB, "3"), 
                                                               "2", clinic__SHLast12MonthsB))) %>%
    dplyr::mutate(clinic__SHSeizComaLast12MonthsB = SHSeizComaLast12MonthsB) %>%
    dplyr::mutate(clinic__DKANumEverB = if_else(str_detect(DKANumEverB, ">"), "6", DKANumEverB)) %>%
    dplyr::mutate(clinic__DKANumEverB = as.numeric(if_else(str_detect(clinic__DKANumEverB, "-"), "5", clinic__DKANumEverB))) %>%
    dplyr::mutate(clinic__DKALast12MonthsB = if_else(DKALast12MonthsB == "", "0", DKALast12MonthsB)) %>%
    dplyr::mutate(clinic__DKALast12MonthsB = as.numeric(if_else(str_detect(clinic__DKALast12MonthsB, ">") |
                                                                  str_detect(clinic__DKALast12MonthsB, "-"), 
                                                                "2", clinic__DKALast12MonthsB))) %>%
    dplyr::mutate(clinic__cgmusestat = if_else(CGMUseStat == "Current", 1, 0)) %>%
    dplyr::mutate(clinic__microvascular_binary = microvascular_binary) %>%
    dplyr::mutate(clinic__microvascular_count = microvascular_count) %>%
    # recode clinic variables
    # convert to cm
    dplyr::mutate(clinic__WaistCir = if_else(WaistCirUnits == "in", WaistCir*2.54, WaistCir)) %>%
    dplyr::mutate(clinic__NeckCir = if_else(NeckCirUnits == "in", NeckCir*2.54, NeckCir)) %>%
    dplyr::mutate(clinic__Height = if_else(HeightUnits == "in", Height*2.54, Height)) %>%
    # convert to kg
    dplyr::mutate(clinic__Weight = if_else(WeightUnits == "lbs", Weight*0.4536, Weight)) %>%
    # compute BMI
    dplyr::mutate(clinic__BMI = (clinic__Weight/clinic__Height^2)*10000) %>%
    # STOP-BANG
    dplyr::mutate(clinic__STOP_snoring_test = if_else(sleepComfort_snoring == "yes", 1, 0)) %>%
    dplyr::mutate(clinic__STOP_tired_test = if_else(sleepComfort_tired == "yes", 1, 0)) %>%
    dplyr::mutate(clinic__STOP_observed_test = if_else(sleepComfort_observed == "yes", 1, 0)) %>%
    dplyr::mutate(clinic__STOP_pressure_test = if_else(sleepComfort_pressure == "yes", 1, 0)) %>%
    dplyr::mutate(clinic__BANG_BMI_test = if_else(clinic__BMI > 35, 1, 0)) %>%
    dplyr::mutate(clinic__BANG_age_test = if_else(demo_age > 50, 1, 0)) %>%
    dplyr::mutate(clinic__BANG_NeckCir_test = if_else(clinic__NeckCir > 40, 1, 0)) %>%
    dplyr::mutate(clinic__BANG_gender_male = if_else(Gender == "M", 1, 0)) %>%
    # compute stop-bang score
    rowwise() %>%
    dplyr::mutate(clinic__STOP_sum = sum(c_across(contains("_STOP_")), na.rm = TRUE)) %>%
    dplyr::mutate(clinic__BANG_sum = sum(c_across(contains("_BANG_")), na.rm = TRUE)) %>%
    dplyr::mutate(clinic__STOP_BANG_sum = sum(c(clinic__STOP_sum, clinic__BANG_sum), na.rm = TRUE)) %>%
    ungroup() %>%
    dplyr::mutate(clinic__STOP_BANG_risk = if_else(clinic__STOP_BANG_sum <= 2, "low", "intermediate")) %>%
    dplyr::mutate(clinic__STOP_BANG_risk = if_else(clinic__STOP_BANG_sum >= 5, "high", clinic__STOP_BANG_risk)) %>%
    dplyr::mutate(clinic__STOP_BANG_risk = if_else(clinic__STOP_sum >= 2 & 
                                                     clinic__BANG_gender_male == 1 &
                                                     !is.na(clinic__BANG_gender_male), "high", clinic__STOP_BANG_risk)) %>%
    dplyr::mutate(clinic__STOP_BANG_risk = if_else(clinic__STOP_sum >= 2 & 
                                                     clinic__BANG_BMI_test == 1 &
                                                     !is.na(clinic__BANG_BMI_test), "high", clinic__STOP_BANG_risk)) %>%
    dplyr::mutate(clinic__STOP_BANG_risk = if_else(clinic__STOP_sum >= 2 & 
                                                     clinic__BANG_NeckCir_test == 1 & 
                                                     !is.na(clinic__BANG_NeckCir_test), "high", clinic__STOP_BANG_risk)
    ) %>% 
    select(user_id, starts_with("clinic__")) 
  colnames(subject_data) <- colnames(subject_data) %>% str_replace_all(., "__", "_")
  
  if (select_SB) {
    subject_data <- 
      subject_data %>% select(user_id, contains("STOP"), contains("BANG"))
  }
  return(subject_data)
}

post_intervals <- function (df) {
  plot <- ggplot(df, aes(x = mean, y = fixef_num, group = inclusion)) +
    geom_vline(xintercept = 0, lty = 2) +
    geom_segment(aes(x = `2.5%`, xend = `97.5%`, 
                     y = fixef_num + jitter, yend = fixef_num + jitter, 
                     color = `95% CI`), size = 1) +
    scale_color_manual(values = c("deepskyblue4", "coral4")) +
    new_scale_color() +
    geom_segment(aes(x = `5%`, xend = `95%`, 
                     y = fixef_num + jitter, yend = fixef_num + jitter,
                     color = `90% CI`), size = 2) +
    scale_color_manual(values = c("deepskyblue3", "coral3")) +
    geom_point(aes(fill = inclusion, y = fixef_num + jitter), shape = 21, size = 3) +
    geom_text(aes(label = sig, y = fixef_num + jitter - .002), size = 3, fontface = "bold") +
    scale_fill_manual(values = c("gray47", "gray73", "white")) +
    theme_classic() +
    labs(fill = "EMA inclusion\ncut-off",
         y = "Model term (linear vs. quadratic)", x = "Magnitude of group effect") +
    guides(fill = guide_legend(order = 1)) 
}

wrangle_func <- function(models, test_outcome) {
  wrangled <-
    models[str_detect(names(models), test_outcome)] %>%
    map(spread_draws, `poly(MinLag_0_WP, 2, raw = FALSE)2`, b[group, term], sep = " u") %>%
    bind_rows(.id = "id") %>%
    rename(MinLag_0_WP_fixed = `poly(MinLag_0_WP, 2, raw = FALSE)2`) %>%
    pivot_wider(names_from = group, values_from = b) %>%
    mutate(term = str_remove(term, "ser_id:")) %>%
    mutate(condition_mean_quad = MinLag_0_WP_fixed + `poly(MinLag_0_WP, 2, raw = FALSE)2`) %>%
    group_by(id, term) %>%
    median_qi(condition_mean_quad, .width = c(.95, .9, .66)) %>%
    ungroup() %>%
    separate(id, into = c("omit0", "inclusion", "test"), sep = "/") %>% 
    separate(inclusion, into = c("omit1", "inclusion", "omit2"), sep = "_") %>% 
    separate(test, into = c("omit3", "omit4", "test", "outcome"), sep = "_") %>% 
    select(-contains("omit")) %>%
    mutate(outcome = if_else(outcome == "row1.rds", "accuracy", "speed"),
           test_outcome = paste0(test, "_", outcome),
           color = "No (66-95%)",
           color = if_else(((.lower < 0 & .upper < 0) | (.lower > 0 & .upper > 0)) & .width == .95, "Yes (95%)", color),
           color = if_else(((.lower < 0 & .upper < 0) | (.lower > 0 & .upper > 0)) & .width == .9, "Yes (90%)", color),
           color = if_else(((.lower < 0 & .upper < 0) | (.lower > 0 & .upper > 0)) & .width == .66, "Yes (66%)", color)) %>%
    group_by(test_outcome) %>%
    nest() %>%
    ungroup()
  
  return(wrangled)
}

cat_plot_func <- function(df) {
  df_new <- 
    df %>%
    mutate(plots = map(.x = data, .y = test_outcome,
                     ~ggplot(data = .x, aes(x = reorder(term, condition_mean_quad), y = condition_mean_quad, 
                                            ymin = .lower, ymax = .upper)) +
                       geom_pointinterval(size = 2, aes(color = color)) +
                       geom_vline(xintercept = 0, lty = 2) +
                       theme_bw() +
                       labs(x = "Participants (ordered by individual estimate of glucose acceleration)", 
                            y = "Individual estimate of glucose acceleration", 
                            color = "CI excludes 0", fill = "Random intercept",
                            title = toupper(.y)) +
                       theme(axis.text.x=element_blank(),  
                             axis.ticks.x=element_blank()) +
                       scale_color_manual(values = c("gray73", "deepskyblue3", "coral3", "darkred")) +
                       facet_wrap(~inclusion, nrow = 3) +
                       guides(color = 'none')
  ))
  
  return(df_new)
}

summary_stats_func <- function(df) {
  df_new <- 
    df$data[[1]] %>% 
    select(condition_mean = condition_mean_quad, inclusion) %>%
    group_by(inclusion) %>% 
    summarise(mean = mean(condition_mean), 
              sd = sd(condition_mean), 
              min = min(condition_mean), 
              max = max(condition_mean)) %>%
    mutate(outcome = df$test_outcome[1])
  
  return(df_new)
}

percent_strong_func <- function(df) {
  df_new <-
    df$data[[1]] %>% 
    mutate(color = if_else(str_detect(color, "Yes"), 1, 0)) %>%
    group_by(term) %>%
    summarise(color = sum(color)) %>%
    ungroup() %>%
    mutate(color = if_else(color >= 1, 1, 0)) %>%
    mutate(outcome = df$test_outcome[1])
  
  return(df_new)
}

quad_plots <- function(model_input, 
                       EMA_cutoffs = c("66pct"), 
                       TMB_list = list(NULL),
                       lasso = FALSE,
                       lasso_results = NULL,
                       n_reps = NA,
                       make_plots = TRUE) {
  
  for (i in EMA_cutoffs) {
    stanmodel <- model_input[str_detect(names(model_input), i)]
    user_id <- stanmodel[[1]]$data$user_id %>% unique(.) 
    epreds <- expand.grid(MinLag_0_WP = seq(-4, 4, by = .05),
                          user_id = user_id) %>%
      as_tibble(.) %>%
      add_epred_draws(stanmodel[[1]], ndraws = 1000)
    
    user_acceleration <- 
      stanmodel[[1]] %>% 
      spread_draws(`(Intercept)`, `poly(MinLag_0_WP, 2, raw = FALSE)2`, b[group, term], sep = " u") %>% 
      dplyr::rename(Intercept = `(Intercept)`, quad_fixed = `poly(MinLag_0_WP, 2, raw = FALSE)2`) %>%
      filter(group == "poly(MinLag_0_WP, 2, raw = FALSE)2" | group == "(Intercept)") %>%
      pivot_wider(names_from = group, values_from = b) %>%
      dplyr::mutate(user_acceleration = quad_fixed + `poly(MinLag_0_WP, 2, raw = FALSE)2`) %>%
      dplyr::mutate(user_intercept = Intercept + `(Intercept)`) %>% 
      group_by(term) %>% median_qi(user_intercept, user_acceleration) %>% 
      separate(term, into = c("omit", "user_id"), sep = ":") %>% 
      select(-omit)
    
    if (i == "66pct") {
      epred_summary <- 
        epreds %>% 
        group_by(MinLag_0_WP, user_id) %>% 
        dplyr::summarise(mean = mean(.epred)) %>%
        ungroup() %>% 
        left_join(user_acceleration %>% ungroup()) %>%
        mutate(mean_centered = mean - user_intercept) %>%
        mutate(inclusion = i)
    } else {
      tmp <- 
        epreds %>% 
        group_by(MinLag_0_WP, user_id) %>% 
        summarise(mean = mean(.epred)) %>%
        ungroup() %>% 
        left_join(user_acceleration %>% ungroup()) %>%
        mutate(mean_centered = mean - user_intercept) %>%
        mutate(inclusion = i)
      epred_summary <- bind_rows(epred_summary, tmp)
    }
  }
  
  if (!make_plots) {
    print("Returning epred_summary")
    return(epred_summary)
  }
  
  individual_plots <- list()
  legend <- list()
  metric_list <- c("")
  if (lasso) {
    metric_list <- 
      c("", 
        lasso_results %>% 
          filter((coef_n > n_reps) & var != ("(Intercept)")) %>%
          dplyr::select(var, inclusion, coef = coef_mean) %>%
          pivot_wider(names_from = inclusion, values_from = coef) %>% 
          ungroup() %>%
          filter(complete.cases(.)) %>% 
          dplyr::mutate(var = str_remove(var, "_z")) %>% pull(var))
  } 
  
  if (length(TMB_list) != length(metric_list)) {
    print("ERROR")
  }
  
  for (j in 1:length(TMB_list)) {
    TMB <- TMB_list[[j]]
    metric <- metric_list[j]
    
    text <- ""
    epred_summary_tmp <- epred_summary 
    if (!is.null(TMB)) {
      metric_vec <- TMB %>% pull(metric)
      TMB <- TMB %>% mutate(metric = metric_vec)
      epred_summary_tmp <- epred_summary %>%
        left_join(TMB %>% select(user_id, metric))
      text <- metric
      print(text)
    }
    
    individual_plots[[j]] <-
      ggplot(data = epred_summary_tmp, aes(x = MinLag_0_WP, y = mean_centered)) +
      {if(is.null(TMB)) geom_line(aes(group = user_id), color = "gray55", alpha = .4) } +
      {if(!is.null(TMB)) geom_line(aes(group = user_id, color = metric), alpha = .4) } +
      stat_smooth(color = "black", method = "lm", formula = y ~ x + I(x^2), aes(color = ntile)) +
      labs(x = "Glucose (centered, scaled within-person)", 
           y = paste0("DSM RT\n(centered around person-level intercept)"), 
           color = text %>% str_remove_all(., "clinic_") %>% str_remove_all(., "demo_")) +
      theme_bw() +
      guides(fill = 'none') +
      {if(!is.null(TMB)) scale_color_gradient2(low = 'red', mid = "salmon", high = "blue",
                                              guide = guide_colorbar(direction = "horizontal", 
                                                                     title.position = "top"),
                                              midpoint = median(epred_summary_tmp$metric, na.rm = TRUE)) } +
      {if(!is.null(TMB)) ggtitle(text %>% str_remove_all(., "clinic_") %>% str_remove_all(., "demo_")) } +
      {if(metric == "clinic_gluBelow70") scale_color_gradient2(low = 'red', mid = "salmon", high = "blue",
                                                               guide = guide_colorbar(direction = "horizontal", 
                                                                                      title.position = "top"),
                                                               limits = c(0, 10), oob = squish,
                                                               midpoint = median(epred_summary_tmp$metric, na.rm = TRUE)) } +
      {if(metric == "clinic_gluCV") scale_color_gradient2(low = 'red', mid = "salmon", high = "blue",
                                                          guide = guide_colorbar(direction = "horizontal", 
                                                                                 title.position = "top"),
                                                          limits = c(20, 50), oob = squish,
                                                          midpoint = median(epred_summary_tmp$metric, na.rm = TRUE)) }
    
    legend[[j]] <- get_legend(individual_plots[[j]])
    individual_plots[[j]] <- individual_plots[[j]] + guides(color = 'none')
  }
  
  return(list(legend, individual_plots))
}

plot_min_max <- function(epred_summary, optimal = c("min", "max"), text) {
  if (optimal == "min") {
    optimal_func <- function(x) {
      min(x, na.rm = TRUE)
    } 
  } else {
    optimal_func <- function(x) {
      max(x, na.rm = TRUE)
    }
  }
  
  minima_maxima <- 
    epred_summary %>% group_by(user_id) %>% 
    filter(mean == optimal_func(mean)) %>% ungroup()
  
  minima_maxima <- minima_maxima %>% 
    mutate(percent_deviation = 100*((mean - user_intercept)/user_intercept),
           deviation = mean - user_intercept,
           ntile_x = ntile(MinLag_0_WP, 100), 
           ntile_y = ntile(percent_deviation, 100)) 
  
  summary_stats_tmp <-
    minima_maxima %>% 
    select(contains("deviation"), MinLag_0_WP, user_id, inclusion) %>%
    group_by(inclusion) %>%
    summarise_if(is.numeric, mean) 
  
  plotting_minima_maxima <-
    ggplot(minima_maxima, aes(x = MinLag_0_WP, 
                              y = percent_deviation)) +
    theme_bw() +
    stat_density_2d_filled(contour_var = 'ndensity', alpha = .8) +
    geom_hline(yintercept = 0, size = 1) +
    geom_vline(xintercept = 0, size = 1) +
    geom_point(alpha = .8, size = .2) +
    facet_wrap(~inclusion) +
    labs(x = "Glucose (centered, scaled within-person)", 
         title = "", 
         y = paste0(text, "\n% deviation from WP average"),
         fill = "Density") +
    scale_fill_manual(values = c("white", brewer.pal(n = 9, "Reds"))) +
    theme(legend.position = "bottom") 
  
  return(list(plotting_minima_maxima, minima_maxima))
}

plot_graph <- function (df, decile, color_key, out_dir, 
                        with_names = TRUE, user_id) {
  c_scale <- colorRamp(c('pink', "red", "darkred"))
  outcome <- "MinLag_0_WP"
  outcome_BP <- "MinLag_0"
  
  user_avgs_BP <- df %>% 
    select((!contains("_WP") & contains("survey")), value, 
           all_of(outcome_BP), -all_of(outcome)) %>%
    summarise_all(mean, na.rm = TRUE) %>%
    gather(key = key, value = num)
  
  n_reps <- 1000
  for (boot in 0:n_reps) {
    user_avgs_WP <- df %>%
      select((contains("_WP") & contains("survey")), value_WP, all_of(outcome)) 
    
    if (boot > 0) {
      user_avgs_WP <- user_avgs_WP %>% slice_sample(prop = 1, replace = TRUE)
    } 
    
    my_mat <- matrix(0, ncol(user_avgs_WP), ncol(user_avgs_WP))
    for (g in 1:ncol(user_avgs_WP)) {
      for (h in 1:ncol(user_avgs_WP)) {
        my_mat[g, h] <- mutinformation(infotheo::discretize(user_avgs_WP[, g]), 
                                       infotheo::discretize(user_avgs_WP[, h]))
      }
    }
    rownames(my_mat) <- colnames(user_avgs_WP)
    colnames(my_mat) <- colnames(user_avgs_WP)
    my_mat[lower.tri(my_mat, diag = TRUE)] <- NA
    
    comput_cors_tmp <- 
      my_mat %>%
      as.data.frame(.) %>% 
      rownames_to_column() %>% 
      gather(key = key, value = variable, -rowname) %>% 
      mutate(variable = if_else(is.na(variable), 0, variable)) %>%
      mutate(retain = ntile(variable, decile), 
             retain = 
               if_else((rowname == "value_WP" & key == outcome & variable > 0), as.integer(decile), retain)) %>%
      mutate(variable = if_else(retain != decile, 0, variable)) %>%
      select(-retain) %>% 
      pivot_wider(names_from = key, values_from = variable) 
    
    if (boot == 0) {
      compute_cors_long <- comput_cors_tmp %>% gather(key = key, value = value, -rowname)
    } else {
      comput_cors_tmp_long <- comput_cors_tmp %>% gather(key = key, value = value, -rowname)
      colnames(comput_cors_tmp_long)[3] <- paste0("value_", boot)
      compute_cors_long <-
        compute_cors_long %>% left_join(comput_cors_tmp_long)
    }
    print(boot)
  }
  
  bootstrapped_paths <- 
    compute_cors_long %>% 
    mutate_at(vars(contains("value")), ~if_else(. > 0, 1, 0)) %>%
    mutate(sum = rowSums(select(., contains("value")))) %>%
    select(rowname, key, sum) %>%
    arrange(desc(sum)) %>%
    filter(sum > n_reps/2)
  
  mean_pathsA <- 
    compute_cors_long %>% 
    gather(key = rep, value = value, -rowname, -key) %>%
    mutate(value = if_else(value == 0, NA_real_, value)) %>%
    group_by(rowname, key) %>% 
    summarise(mean = mean(value, na.rm = TRUE)) %>%
    right_join(bootstrapped_paths) %>%
    select(1:3) 
  mean_pathsB <- data.frame(rowname = mean_pathsA$key,
                            key = mean_pathsA$rowname, 
                            mean = mean_pathsA$mean) %>%
    bind_rows(mean_pathsA) %>%
    pivot_wider(names_from = key, values_from = mean) %>%
    {. ->> tmp} %>%
    select(rowname, tmp$rowname)
  
  my_mat2 <- mean_pathsB %>% select(-rowname) %>% as.matrix(.)
  rownames(my_mat2) <- mean_pathsB$rowname
  
  if (sum(rownames(my_mat2) %in% colnames(my_mat2)) != length(rownames(my_mat2))) {
    print("Warning: rownames != colnames in matrix conversion")
  }
  
  g <- graph_from_adjacency_matrix(my_mat2, mode = "undirected", 
                                   weighted=TRUE, diag = FALSE)
  g <- delete.vertices(g, which(igraph::degree(g)==0))
  g <- delete.edges(g, which(is.na(E(g)$weight)))
  lay <- layout_with_fr(g)
  vert_labels <-
    data.frame(key = V(g)$name) %>% 
    mutate(key_orig = key,
           key = str_remove_all(key, "_WP"),
           key = str_remove_all(key, "_sq")) %>% 
    left_join(user_avgs_BP) %>%
    mutate(key = str_remove(key, "survey_"),
           key = str_remove(key, "_WP"),
           key = str_remove(key, "_score"),
           key = str_remove(key, "_numeric"),
           key = if_else(key == "value", "DSM_medRT", key),
           num = as.character(round(num, 2)),
           label = paste0(key, "\n", num)) %>%
    left_join(color_key %>% select(key_orig = var, col)) 
  for (m in 1:nrow(vert_labels)) {
    V(g)[vert_labels$key_orig[m]]$color <- vert_labels$col[m]
    if (with_names) {
      V(g)[vert_labels$key_orig[m]]$name <- vert_labels$label[m]
    } else {
      # V(g)[vert_labels$key_orig[m]]$name <- vert_labels$num[m]
      V(g)[vert_labels$key_orig[m]]$name <- ""
    }
  }
  
  index_gl <-
    vert_labels %>% 
    mutate(rownum = row_number()) %>%
    filter(key == "DSM_medRT" | key == "MinLag_0") %>%
    pull(rownum)
  index_other <-
    vert_labels %>% 
    mutate(rownum = row_number()) %>%
    filter(key != "DSM_medRT" & key != "MinLag_0") %>%
    pull(rownum)
  vert_labels$path <- NA
  for (m in index_other) {
    path_length_1 <-
      distances(g, 
                v = V(g)[index_gl[1]],
                to = V(g)[m], 
                algorithm = "unweighted")[[1]]
    path_length_2 <-
      distances(g, 
                v = V(g)[index_gl[2]],
                to = V(g)[m], 
                algorithm = "unweighted")[[1]]
    vert_labels$path[m] <- min(c(path_length_1, path_length_2))
  }
  
  png(filename = paste0(out_dir, i, "_withNames", with_names, ".png"), 
      width = 4, height = 5, units = "in", res = 300)
  plot(g, 
       layout=lay, 
       edge.width = (E(g)$weight + 1)*5,
       edge.color = apply(c_scale(scales::rescale(E(g)$weight, to = c(0, 1))), 1, 
                          function(x) rgb(x[1]/255,x[2]/255,x[3]/255)),
       vertex.label.family = "Helvetica",
       vertex.label.color = "black",
       vertex.size = 20,
       edge.label.color = "darkred",
       main = user_id)
  dev.off()
  return(vert_labels)
}


