# get data ----------
get_result_data <- function(path = NULL) {
  
  if (is.null(path)) {
    path <- "/Users/maxsalvatore/Dropbox (University of Michigan)/epiq/EPIQ-Data-09242020/EPIQ_SurveyQuestionResults.xlsx"
    message("Using default path: '/Users/maxsalvatore/Dropbox (University of Michigan)/epiq/EPIQ-Data-09242020/EPIQ_SurveyQuestionResults.xlsx'")
  }
  
  readxl::read_xlsx(path) %>%
    janitor::clean_names()
}

# get info -----------
get_life_meaning_info <- function() {
  
  read_tsv(here::here("data", "life_meaning_info.txt"), col_types = cols())
  
}

# clean data ----------
clean_lm_data <- function(dat) {
  
  dat %>%
    mutate(
      values = dplyr::recode(answers, !!!lm_labs)
    ) %>%
    dplyr::select(
      id       = participant_identifier,
      question = result_identifier,
      answers,
      values
    ) %>%
    mutate(
      question  = as.factor(question)
    ) %>%
    left_join(
      lm_info, by = "question"
    )
  
}

clean_covariates <- function(dat, questions = NULL) {
  
  if (is.null(questions)) {
    message("selecting default questions")
    questions <- c(
      "house_type",
      "income",
      "employment_status",
      "education",
      "race_eth",
      "sex",
      "gender",
      "cancer_ever"
    )
  }
  
  tmp_age <- dat %>%
    group_by(id) %>%
    summarize(
      age = round(min(start_dsb)/365)
    )
  
  dat %>%
    dplyr::select(id, res_id, answer) %>%
    mutate(res_id  = as.factor(res_id)) %>%
    dplyr::filter(
      res_id %in% questions
    ) %>%
    pivot_wider(
      names_from  = res_id,
      values_from = answer,
      id_cols     = id
    ) %>%
    mutate(
      cancer_ever = case_when(
        tolower(cancer_ever) == "true" ~ "Yes",
        tolower(cancer_ever) == "false" ~ "No",
        T ~ cancer_ever
      )
    ) %>%
    left_join(tmp_age, by = "id")
  
}

get_subdim_scores <- function(dat) {
  
  dat %>%
    group_by(id, dimension, sub_dimension) %>%
    summarize(
      mean = mean(values, na.rm = TRUE),
    ) %>%
    ungroup() %>%
    mutate(
      tmp_name = paste0(dimension, "_", sub_dimension)
    ) %>%
    dplyr::select(id, tmp_name, mean) %>%
    pivot_wider(
      names_from  = tmp_name,
      values_from = mean,
      id_cols     = id
    ) %>%
    dplyr::select(-id)
  
}

# quick descriptive summary ---------
quick_descriptive_summary <- function(dat) {
  
  dat %>%
    group_by(question) %>%
    summarize(
      n    = sum(!is.na(values)),
      n_miss = sum(is.na(values)),
      mean = mean(values, na.rm = T),
      sd   = sd(values, na.rm = T),
      min  = min(values, na.rm = T),
      max  = max(values, na.rm = T),
      .groups = "drop_last"
    )
  
}

# get itemwise corplot -----------
get_itemwise_corplot <- function(dat) {
  
  dat %>%
    dplyr::select(id, question, values) %>%
    pivot_wider(
      names_from = question,
      values_from = values,
      id_cols = id
    ) %>%
    dplyr::select(-id) %>%
    cor(use = "complete.obs") %>%
    ggcorrplot(
      hc.order    = TRUE,
      outline.col = "white",
      type        = "lower",
      ggtheme     = ggplot2::theme_minimal,
      lab         = TRUE
    ) +
    labs(
      title = "CMM itemwise correlation"
    ) +
    lm_theme
  
}

# get sub-dimension-wise corplot ----------
get_subdimwise_corplot <- function(dat) {
  
  dat %>%
    group_by(id, dimension, sub_dimension) %>%
    summarize(
      mean = mean(values, na.rm = TRUE),
    ) %>%
    ungroup() %>%
    mutate(
      tmp_name = paste0(dimension, "_", sub_dimension)
    ) %>%
    dplyr::select(id, tmp_name, mean) %>%
    pivot_wider(
      names_from  = tmp_name,
      values_from = mean,
      id_cols     = id
    ) %>%
    dplyr::select(-id) %>%
    cor(use = "complete.obs") %>%
    ggcorrplot(
      hc.order    = TRUE,
      outline.col = "white",
      type        = "lower",
      ggtheme     = ggplot2::theme_minimal,
      lab         = TRUE
    ) +
    labs(
      title = "CMM subdimension-wise correlation",
      caption = "**Note:** Pearson correlation by simple average of within-subject subdimension items."
    ) +
    lm_theme
  
}

# get cmm score -----------
get_cmm_score <- function(dat) {
  
  dat %>%
    group_by(id) %>%
    summarize(
      cmm_score = sum(values, na.rm = TRUE),
      .groups   = "drop_last"
    ) %>%
    ungroup()
  
}

# cmm score distribution plot ----------
get_cmm_distplot <- function(dat, bins = 30) {
  
  dat %>%
    ggplot(aes(x = cmm_score)) +
    geom_histogram(bins = bins) +
    labs(
      title   = "Distribution of total CMM score",
      x       = "CMM Score",
      y       = "Count",
      caption = glue("**Note:** n = {dim(dat)[1]}; range from {min(dat$cmm_score, na.rm = TRUE)} to {max(dat$cmm_score, na.rm= TRUE)} out of a possible 147")
    ) +
    theme_minimal() +
    lm_theme
  
}

# cmm score distribution by variables -----------
get_cmm_distby <- function(dat, var, desc) {
  
  n_tot <- dim(dat %>% filter(!is.na(cmm_score)))[1]
  tmp   <- dat %>% filter(!is.na(cmm_score) & !is.na({{ var }}))
  n_miss <- n_tot - dim(tmp)[1]
  
  tmp_aov <- aov(tmp %>% pull(cmm_score) ~ tmp %>% pull({{ var }}))
  anova_p <- summary(tmp_aov)[[1]][["Pr(>F)"]][1]
  
  tmp %>%
    ggplot(aes(x = {{ var }}, y = cmm_score)) +
    geom_boxplot() +
    stat_summary(fun.data = give_n, geom = "text") +
    coord_flip() +
    labs(
      title = glue("Distribution of CMM by {desc}"),
      x     = desc,
      y     = "CMM Score",
      caption = glue("**ANOVA p-value:** {round(anova_p, 4)}<br>",
                     "**Note:** n = {dim(tmp)[1]} ({n_miss} missing {desc} data)")
    ) +
    theme_minimal() +
    lm_theme
  
}

# give n -----------
  # from: https://stackoverflow.com/questions/15660829/how-to-add-a-number-of-observations-per-group-and-use-group-mean-in-ggplot2-boxp
give_n <- function(x, mult = 1.05, man = TRUE, man_y = 150){
  
  if (man == FALSE) {
    return(c(y = median(x) * mult, label = length(x))) 
  }
  
  if (man == TRUE) {
    return(c(y = man_y, label = length(x))) 
  }
  
}

# lm_theme ----------
lm_theme <- ggplot2::theme(
    text          = element_text(family = "Lato"),
    plot.title    = element_text(face = "bold", size = 18),
    plot.subtitle = element_text(size = 14),
    axis.title    = element_text(size = 12, face = "italic"),
    axis.text     = element_text(size = 12),
    plot.caption  = element_markdown(hjust = 0, size = 11),
    panel.grid    = element_blank(),
    axis.line     = element_line(size = 0.25, color = "black"),
    panel.background = element_blank()
  )

# life meaning values -----------
lm_labs <- c(
  "Strongly disagree"         = 1,
  "Disagree"                  = 2,
  "Slightly disagree"         = 3,
  "Neither agree or disagree" = 4,
  "Slightly agree"            = 5,
  "Agree"                     = 6,
  "Strongly agree"            = 7
)

# scatterplot ----------
get_scatterplot <- function(dat, x_var, y_var, ..., ann_x = 6.5, ann_y = 1.5) {
  
  tmp_cor <- cor.test(dat %>% pull({{ x_var }}), dat %>% pull({{ y_var }}))
  
  dat %>%
    drop_na({{ x_var }}, {{ y_var }}) %>%
    ggplot(aes(x = {{ x_var }}, y = {{ y_var }})) +
    geom_hline(yintercept = dat %>% pull({{ y_var }}) %>% median(na.rm = TRUE), color = "darkgray") +
    geom_vline(xintercept = dat %>% pull({{ x_var }}) %>% median(na.rm = TRUE), color = "darkgray") +
    geom_smooth(method = "loess", se = FALSE, formula = "y ~ x", span = 0.5) + 
    geom_jitter() +
    annotate(geom  = "text",
             label = glue("Correlation = {round(tmp_cor$estimate, 2)}\nP = {formatC(tmp_cor$p.value, format = 'e', digits = 2)}"),
             x     = ann_x,
             y     = ann_y,
             hjust = 0.5
    )+
    labs(...,
         caption = glue("**Note:** Simple average within each subdomain used. Guidelines represent axis median.")
    ) +
    theme_minimal() +
    lm_theme
  
}