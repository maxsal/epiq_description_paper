# script used to prepare:
#     ~/projects/epiq/linkage_test/data/magic_data_20210223.rds

library("glue")

root <- "~/projects/epiq/linkage_test/"

source(glue("{root}libraries.R"))
source(glue("{root}functions.R"))

# long format ----------
survey <- read_tsv("~/projects/epiq/data/magic_survey_data_20210223.txt",
                   col_types = cols()) %>%
  janitor::clean_names() %>%
  dplyr::rename(
    id          = de_id_patient_id,
    module      = survey_name,
    version     = survey_version,
    step        = step_identifier,
    res_id      = result_identifier,
    step_type   = survey_step_type,
    task_status = survey_task_status,
    start_dsb   = survey_start_date_days_since_birth,
    end_dsb     = survey_end_date_days_since_birth,
    question    = survey_question,
    answer      = survey_answer
  )

# pull duplicate IDs for GENDER res_id
#   appears 8 participants responded to this question twice,
#   Male and M or Female and F
#   dropping responses that are M or F
gender_dup_ids <- survey %>%
  dplyr::select(id, res_id, answer) %>%
  distinct() %>%
  group_by(id, res_id) %>%
  mutate(count = n()) %>%
  filter(count > 1 & res_id == "GENDER") %>%
  ungroup() %>%
  pull(id) %>%
  unique()

# upon previously exploration, one participant appears to
#   have answered the occupational exposure module twice
#   with generally inconsistent responses
#   dropping occupational exposure module responses for
#   this participant
other_dup_ids <- survey %>%
  dplyr::select(id, res_id, answer) %>%
  distinct() %>%
  group_by(id, res_id) %>%
  mutate(count = n()) %>%
  filter(count > 1 & res_id != "GENDER") %>%
  ungroup() %>%
  pull(id) %>%
  unique()

survey <- survey %>%
  # drop duplicate GENDER responses for 8 individuals (all concordant)
  dplyr::filter(!(res_id == "GENDER" & (id %in% gender_dup_ids) & (answer %in% c("M", "F")))) %>%
  # drop occupational exposure module responses for individual who answered twice (inconsistent answers)
  dplyr::filter(!(id %in% other_dup_ids & module == "(EPI-Q Optional) Occupational Exposure"))

# save cleaned file
readr::write_rds(survey, glue("{root}data/magic_data_20210223.rds"), compress = "gz")

# wide format ------------
cov <- survey %>%
  clean_covariates()

survey_wide <- survey %>%
  dplyr::select(id, res_id, answer) %>%
  distinct() %>%
  tidyr::pivot_wider(
    names_from  = res_id,
    values_from = answer
  ) %>%
  left_join(cov %>% dplyr::select(id, age), by = "id") %>%
  janitor::clean_names() %>%
  dplyr::filter(id %in% cov$id)

write_tsv(survey_wide, glue("{root}data/tmp_wide.txt"))
tmp <- read_tsv(glue("{root}data/tmp_wide.txt"), col_types = cols())

readr::write_rds(tmp, glue("{root}data/magic_data_wide_20210223.rds"), compress = "gz")
file.remove(glue("{root}data/tmp_wide.txt"))

