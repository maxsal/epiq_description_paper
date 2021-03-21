library("glue")

root <- "~/projects/epiq/linkage_test/"

source(glue("{root}libraries.R"))
source(glue("{root}functions.R"))

# load epiq data -----------
survey <- readr::read_rds(glue("{root}data/magic_data_20210223.rds"))
wide   <- readr::read_rds(glue("{root}data/magic_data_wide_20210223.rds"))

cov <- survey %>%
  clean_covariates()
length(unique(cov$id))

# load magic data ----------
magic_pim <- fread("/net/junglebook/magic_data/MGI_Phenome_20200323/MGI_20200323_PEDMASTER.txt")
magic_pim <- setnames(magic_pim, "IID", "id")

magic_cov <- fread("/net/junglebook/magic_data/genotypes_20200316/MGI_sample_info_20200316.txt")
setnames(magic_cov, "IID", "id")
setnames(magic_cov, "INFERRED_ANCESTRY", "ancestry_geno")
magic_cov <- magic_cov[, sex_geno := ifelse(INFERRED_SEX == 1, "Male", ifelse(INFERRED_SEX == 2, "Female", NA))][
  , .(id, sex_geno, ancestry_geno)
]

super <- as.data.table(merge(wide, magic_pim, all.x = TRUE))
super <- as_tibble(merge(super, magic_cov, all.x = TRUE))

table(super$sex, super$sex_geno, useNA = "ifany")
table(super$gender, super$sex_geno, useNA = "ifany")
table(super$race_eth, super$ancestry_geno, useNA = "ifany")

as.data.table(super)[, .(id, sex, sex_geno, gender)
                     ][, agree := as.numeric(sex == sex_geno)
                       ][agree != 1,]

write_rds(super, "~/projects/epiq/linkage_test/data/super.rds", compress = "gz")

dist_plot <- function(data, x = age, y = NULL, color = NULL, type = "hist") {
  
  if (type == "hist" & is.null(y)) {
    
    tmp <- data %>%
      ggplot(aes(x = {{ x }}, color = {{ color }})) +
      geom_histogram() +
      labs(
        x = "Age",
        y = "Count",
        # title = glue("{x} by {color}"),
        caption = glue("\uA9 MAGIC {format(Sys.Date(), '%e %B %Y')}")
      ) +
      lm_theme
    
    ggsave(
      tmp,
      filename = glue("{root}fig/{x}_x_{color}_hist.pdf"),
      width  = 7,
      height = 5,
      device = cairo_pdf
    )
    
    return(TRUE)
    
  }
  
  if (type == "density" & is.null(y)) {
    
    tmp <- data %>%
      ggplot(aes(x = {{ x }}, color = {{ color }})) +
      geom_density() +
      labs(
        x = "Age",
        y = "Count",
        title = glue("{x} by {color}"),
        caption = glue("\uA9 MAGIC {format(Sys.Date(), '%e %B %Y')}")
      ) +
      lm_theme
    
    ggsave(
      tmp,
      filename = glue("{root}fig/{x}_x_{color}_density.pdf"),
      width  = 7,
      height = 5,
      device = cairo_pdf
    )
    
    return(TRUE)
    
  }
  
  return(FALSE)
  
}

super %>% dist_plot(x = "age", color = "cancer_ever", type = "density")

# load comorbidity codes -----------
com <- read_tsv("/net/junglebook/michiganmedicine/larsf/data/comorbidities/ComoScore_Phecodes.txt",
                col_types = cols())

taco <- super %>%
  select(c("id", "cancer_ever", com %>% filter(AnyCancer == 1) %>% pull(phecode) %>% paste0("X", .))) %>%
  replace(is.na(.), 0) %>%
  mutate(cancer_pheno = rowSums(across(com %>% filter(AnyCancer == 1) %>% pull(phecode) %>% paste0("X", .)))) %>%
  mutate(cancer_ehr   = as.numeric(cancer_pheno > 0)) %>%
  select(id, cancer_ever, cancer_pheno, cancer_ehr)
    
table(taco$cancer_ever, taco$cancer_ehr, useNA = "ifany")

check_ids <- taco %>%
  dplyr::filter(cancer_ever == FALSE & cancer_ehr == 1) %>%
  pull(id)
     

tmp <- super %>%
  filter(id %in% check_ids) %>%
  dplyr::select(com %>% filter(AnyCancer == 1) %>% pull(phecode) %>% paste0("X", .)) %>%
  replace(is.na(.), 0) %>%
  sapply(., function(x) sum(x, na.rm = T)) %>% t() %>% t() %>% as.data.frame()

tmp$phecode <- rownames(tmp)

tmp %>%
  as_tibble() %>%
  dplyr::select(
    phecode,
    n = V1
  ) %>%
  arrange(desc(n))

tmp2 <- super %>%
  dplyr::select("id", com %>% filter(AnyCancer == 1) %>% pull(phecode) %>% paste0("X", .)) %>%
  pivot_longer(
    names_to = "phecode",
    values_to = "value",
    -id
  ) %>%
  dplyr::mutate(
    weird = case_when(
      id %in% check_ids ~ "strange",
      T ~ "fine"
    )
  ) %>%
  group_by(weird, phecode) %>%
  summarize(n = sum(value, na.rm = T)) %>%
  pivot_wider(
    names_from = "weird",
    values_from = "n"
  ) %>%
  arrange(desc(fine)) %>%
  left_join(com %>% 
              mutate(phecode = paste0("X", phecode)
                     ) %>%
              dplyr::select(phecode, description),
            by = "phecode") %>%
  dplyr::select(phecode, description, fine, strange)

write_tsv(tmp2, glue("{root}tmp_cancer_comp.txt"))

super %>%
  replace(is.na(.), 0) %>%
  mutate(sum = rowSums(across(com %>% filter(Comorbidity == "cancer" %>% pull(xphecode))))) %>%
  select(id, cancer_ever, sum)
super <- super %>%
  dplyr::mutate(
    cancer_pheno = 
  )

cancer_codes <- com %>% filter(Comorbidity == "cancer") %>% dplyr::select(xphecode, description)
write_tsv(cancer_codes, glue("{root}tmp_cancer.txt"))

tmp2 <- super %>%
  select(c(com %>% filter(Comorbidity == "cancer") %>% pull(xphecode))) %>%
  replace(is.na(.), 0) %>%
  sapply(., function(x) sum(x, na.rm = T)) %>%
  t() %>% t() %>% as.data.frame()
tmp2$xphecode <- rownames(tmp2)
tmp3 <- tmp2 %>% as_tibble() %>% dplyr::select(xphecode, n = V1) %>% arrange(desc(n))

cancer_codes <- cancer_codes %>% left_join(tmp3, by = "xphecode")
write_tsv(cancer_codes, glue("{root}tmp_cancer.txt"))


# cancer_ever by age -----------
p <- super %>%
  drop_na(cancer_ever) %>%
  ggplot(aes(x = age, color = cancer_ever)) +
  geom_density() +
  labs(
    x = "Age", 
    y = "Density",
    title = "Age by self-reported cancer status (ever)",
    caption = glue("23 February 2021<br>**\uA9 EPI-Q**")
  ) +
  lm_theme
  ggsave(p, filename = glue("{root}fig/age_x_cancer_ever_density.pdf"), width = 9, height = 6, device = cairo_pdf)
  
  p <- super %>%
    drop_na(cancer_ever) %>%
    ggplot(aes(x = age, fill = cancer_ever)) +
    geom_histogram(alpha=0.5, position="identity") +
    labs(
      x = "Age", 
      y = "Density",
      title = "Age by self-reported cancer status (ever)",
      caption = glue("23 February 2021<br>**\uA9 EPI-Q**")
    ) +
    lm_theme
  ggsave(p, filename = glue("{root}fig/age_x_cancer_ever_hist.pdf"), width = 9, height = 6, device = cairo_pdf)
