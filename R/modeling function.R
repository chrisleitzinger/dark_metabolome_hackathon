library(tidyverse)
library(tidymodels)
library(themis)
library(doParallel)

omics_data <- read_delim(paste0(here::here(), "/lusq_reprocess_2021-12-08_iron_log2_merged.txt")) %>% 
  janitor::clean_names()
metabolites_data <- read_delim(paste0(here::here(), "/hmdb_keep_v3_python_blessed.txt")) %>% 
  janitor::clean_names()

omics_data <- omics_data %>% 
  mutate(across((starts_with("x06")), .fns = ~ replace_na(., min(., na.rm = TRUE))))

mldata <- omics_data %>% 
  # Filter identified metabolites
  filter(non_heavy_identified_flag == 1 & !is.na(hmdb) & !str_detect(hmdb, "\\|")) %>% 
  # Join with known metabolites
  inner_join(., metabolites_data %>% 
               select(accession, taxonomy_super_class),
             by = c("hmdb" = "accession")) %>% 
  filter(!is.na(taxonomy_super_class))

test_data <- omics_data %>% 
  # Filter identified metabolites
  filter(non_heavy_identified_flag != 1) %>% 
  select(hmdb, row_id, 
         row_m_z, row_retention_time,
         starts_with("x06"))

train_data <- mldata %>% 
  select(hmdb, row_m_z, row_retention_time, 
         starts_with("x06"), taxonomy_super_class) %>% 
  group_by(taxonomy_super_class) %>% 
  mutate(number_sample_class = n()) %>% 
  ungroup() %>% 
  filter(str_detect(taxonomy_super_class, "Lipids|acids")) %>% 
  select(-number_sample_class) %>% 
  mutate_if(is.character, factor)


model_metabolite <- function(train_data, test_data) {


  mldata_recipe <-
    # 1.model formula
    recipe(taxonomy_super_class ~ ., data = train_data)  %>% 
    # 2.keep these variables but not use them as either outcomes or predictors
    update_role(hmdb, new_role = "ID") %>% 
    # Use nearest neighbor to create new synthetic observation almost similar
    step_smote(taxonomy_super_class)
  set.seed(123)
  mldata_folds <- vfold_cv(train_data, strata = taxonomy_super_class)
  
  set.seed(345)
  
  ranger_spec <- rand_forest(
    mtry = tune(), 
    min_n = tune(),
    trees = tune()) %>% 
    set_mode("classification") %>% 
    set_engine("ranger")
  
  # Set up a workflow container
  ranger_workflow <- 
    workflow() %>% 
    # Add preprocessor
    add_recipe(mldata_recipe) %>%
    add_model(ranger_spec)
  
  # Tuning
  registerDoParallel()
  set.seed(345)
  ranger_tune <-
    tune_grid(ranger_workflow, 
              resamples = mldata_folds, 
              grid = 20)
  set.seed(678)
  final_rf <- ranger_workflow %>% 
    finalize_workflow(select_best(ranger_tune, "accuracy"))
  
  set.seed(101)
  
  last_fit <- final_rf %>% 
    fit(train_data)
  
  test_fit <- last_fit %>% 
    predict(test_data)
  
  predicted_metabolites <- test_data %>% 
    bind_cols(., test_fit) %>% 
    rename(taxonomy_super_class = .pred_class) %>% 
    mutate(pred = "prediction") %>% 
    bind_rows(., train_data %>% 
                mutate(pred = "measured"))
  
  predicted_metabolites %>% 
    ggplot(aes(x=row_m_z, row_retention_time, color= pred))+
    geom_point()+
    theme_minimal()+
    facet_wrap(. ~ taxonomy_super_class, ncol = 2)

}

model_metabolite(train_data = train_data, test_data = test_data)
