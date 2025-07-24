# TITLE: Supplementary statistical analyses
# AUTHOR: Hannah Keefe
# DATE OF LAST EDIT: March 19, 2025
# BACKGROUND:
## for each of the 8 GLMS in the SAR SDM paper, I tested two different transformations of range size (log and sqrt)
## all of the transformations met model assumptions, but we had to choose one, so I chose the one for each model that best fit the assumptions
## comparisons of the two transformations for each model can be found in the Quarto printouts for each model
## the code for those Quarto docs can be found: sarsdm/analyses/build_mX.qmd and sarsdm/analyses/build_mX_summary.qmd
## the statistical analysis using the best fit transformation can be found under 03_statistical_analysis.R, and is run as part of the 04_fig_tables.qmd to create tables and a forest plot for the paper
## this code is almost identical to the 03_statistical_analysis.R code, but it subs the range size transformation to that with a worse fit so that we can provide tables to be included in the supplementary information
# GOALS: 
## 1) run GLMs for SAR SDM literature review paper using the alternate range size transformations
## 2) save model results in tables - these will be included in the supplementary materials of the paper
# NOTES: 
## more details on the models can be found in the quarto model summaries

# Load libraries, set path #====================================================

#load libraries
library(dplyr) #for data management
library(tidyverse) #for data management
library(GGally) #for pairs plots
library(car)
library(performance) #for VIF, testing model assumptions, evaluating model
library(broom) #for looking at model outputs
library(knitr)

# setwd("~/projects/sarsdm")
path = here::here()

# Load and clean database #=====================================================

# load database
sarsdm <- read.csv(paste0(path,"/data/SAR-SDM_database_clean.csv"))

# clean database
sarsdm <- sarsdm %>%
  # remove species without range data, and Extirpated species
  filter (range_data == "1", sara_status_adj != "Extirpated") %>%
  # remove mosses because no sdms so can't model
  filter(taxonomic_group != "Mosses") %>%
  # convert binary response variables into factors
  mutate(year = min_year_of_listing,
         sdm = factor(sdm),
         cc_sdm = factor(cc_sdm),
         sdm_ca = factor(sdm_ca),
         cc_sdm_ca = factor(cc_sdm_ca))%>%
  # set taxonomic group and sara status as factors, and apply reference group
  mutate(taxonomic_group = factor(taxonomic_group,
                                  levels = c("Amphibians", "Arthropods",
                                             "Birds","Lichens","Mammals (terrestrial)",
                                             "Molluscs","Reptiles",
                                             "Vascular Plants"),
                                  labels = c("Amphibians", "Arthropods",
                                             "Birds","Lichens","Mammals",
                                             "Molluscs","Reptiles",
                                             "VascularPlants")),
         taxonomic_group = relevel(taxonomic_group, ref="Arthropods"),
         sara_status_adj = factor(sara_status_adj,
                                  levels = c("Special Concern", "Threatened",
                                             "Endangered")),
         sara_status_adj = relevel(sara_status_adj, ref = "Special Concern")) %>%
  # set cc as a threat to factor
  mutate(cc_threat_ca_adj = factor(cc_threat_ca_adj)) %>%
  # apply transformations to range size
  mutate(log_range_am = log(range_size_America),
         sqrt_range_am = sqrt(range_size_America)) %>%
  select( 
    # species information
    species, range_data,
    # response variables
    sdm, sdm_ca, cc_sdm, cc_sdm_ca, n_sdm, n_cc_sdm, n_sdm_ca, n_cc_sdm_ca, 
    # covariates
    taxonomic_group, sara_status_adj, log_range_am, sqrt_range_am, year, cc_threat_ca_adj) %>%
  unique()

# Model 1: probability of inclusion in a North American SDM #===================

## Notes: Binomial distribution, logit link, excludes Amphibians, log transformation of range size

m1_supp_db <- sarsdm %>%
  filter (taxonomic_group != "Amphibians") %>%
  mutate(taxonomic_group = droplevels(taxonomic_group)) %>%
  select(sdm, taxonomic_group, sara_status_adj, log_range_am, year)

## Pairs plot ##================================================================

m1_supp_db %>% ggpairs()

## Standardize data ##==========================================================

m1_supp_db_std <- model.matrix(sdm ~ taxonomic_group + log_range_am + sara_status_adj + year,
                          data = m1_supp_db) %>% #creates a dataframe with dummy binary variables
  as.data.frame() %>%
  mutate_all(scale) %>% #centres a variable around it's mean (i.e., centers it)
  add_column(sdm = m1_supp_db$sdm) %>%
  rename()

# create version of database with cleaned up column names
m1_supp_db_std_clean <- m1_supp_db_std
colnames(m1_supp_db_std_clean) <- c("Intercept","Birds","Lichens","Mammals","Molluscs","Reptiles","VascularPlants","Rangesize","Threatened","Endangered","Year","sdm")

m1_supp_db <- m1_supp_db %>%
  mutate(log_range_std = log_range_am - mean(log_range_am),
         year_std = year - mean(year))

## Build model ##===============================================================

# for plotting model estimates
m1_supp_std <- glm (sdm ~ 
                 taxonomic_groupBirds +
                 taxonomic_groupLichens +
                 taxonomic_groupMammals +
                 taxonomic_groupMolluscs +
                 taxonomic_groupReptiles +
                 taxonomic_groupVascularPlants +
                 sara_status_adjThreatened + 
                 sara_status_adjEndangered +
                 log_range_am +
                 year,
               family = binomial,
               data = m1_supp_db_std)

# for saving model summary
m1_supp_std_clean <- glm (sdm ~ Birds + Lichens + Mammals +Molluscs + Reptiles +
                       VascularPlants + Threatened + Endangered + Rangesize + Year,
                     family = binomial,
                     data = m1_supp_db_std_clean)

# will use this for plotting model assumptions
m1_supp_ptl <- glm (sdm ~ taxonomic_group + sara_status_adj + log_range_std + year_std,
               family = binomial,
               data = m1_supp_db)

summary(m1_supp_std)

# write csv file with model results
m1_supp_results <- tidy(m1_supp_std_clean) %>%
  mutate(p.value = round(p.value, digits=3),
         estimate = round(estimate, digits =3),
         std.error = round(std.error, digits =3),
         statistic = round(statistic, digits =3)) %>%
  rename(Variables = term,
         "Standardized effects" = estimate,
         "Standard Error" = std.error,
         Statistic = statistic,
         "p-value" = p.value) %>%
  mutate(model = c("m1 Supplementary: North American SDM inclusion"))
write.csv(m1_supp_results, paste0(path,"/outputs/model_results/model1_alt_rangesize.csv"))

## Check assumptions ##=========================================================

# check collinearity, outliers, binned residuals
check_model(m1_supp_ptl)

# check linearity with logit for continuous covariates
op <- par(mfrow = c(1, 2))
car::crPlot(m1_supp_ptl, variable = "log_range_std", id = TRUE)
car::crPlot(m1_supp_ptl, variable = "year_std", id = TRUE)

## Clean environment ##=========================================================

rm(m1_supp_db, m1_supp_db_std, m1_supp_db_std_clean, m1_supp_results, m1_supp_ptl, m1_supp_std, m1_supp_std_clean)

# Model 2: probability of inclusion in North American CC-SDM #==================

## Notes: Binomial distribution, logit link, sqrt transformation of range size

m2_supp_db <- sarsdm %>%
  select (cc_sdm, taxonomic_group, sqrt_range_am, sara_status_adj, year, cc_threat_ca_adj)

## Pairs plot ##================================================================

m2_supp_db %>% ggpairs()

## Standardize data ##==========================================================

m2_supp_db_std <- model.matrix(cc_sdm ~ taxonomic_group + sqrt_range_am + sara_status_adj + year + cc_threat_ca_adj,
                          data = m2_supp_db) %>% #creates a dataframe with dummy binary variables
  as.data.frame() %>%
  mutate_all(scale) %>% #centres a variable around it's mean (i.e., centers it)
  add_column(cc_sdm = m2_supp_db$cc_sdm) %>%
  rename()

m2_supp_db_std_clean <- m2_supp_db_std
colnames(m2_supp_db_std_clean) <- c("Intercept","Amphibians","Birds","Lichens","Mammals","Molluscs","Reptiles","VascularPlants","Rangesize","Threatened","Endangered","Year","CCThreat","cc_sdm")

m2_supp_db <- m2_supp_db %>%
  mutate(sqrt_range_std = sqrt_range_am - mean(sqrt_range_am),
         year_std = year - mean(year))

## Build model ##===============================================================

# for plotting covariates
m2_supp_std <- glm (cc_sdm ~ 
                 taxonomic_groupAmphibians +
                 taxonomic_groupBirds +
                 taxonomic_groupLichens +
                 taxonomic_groupMammals +
                 taxonomic_groupMolluscs +
                 taxonomic_groupReptiles +
                 taxonomic_groupVascularPlants +
                 sara_status_adjThreatened + 
                 sara_status_adjEndangered +
                 sqrt_range_am +
                 year +
                 cc_threat_ca_adj1,
               family = binomial,
               data = m2_supp_db_std)

# for printing summary
m2_supp_std_clean <- glm (cc_sdm ~ Amphibians + Birds + Lichens + Mammals + Molluscs + 
                            Reptiles + VascularPlants + Threatened + Endangered +
                       Rangesize + Year + CCThreat,
                     family = binomial,
                     data = m2_supp_db_std_clean)


# will use this for plotting model assumptions
m2_supp_ptl <- glm (cc_sdm ~ taxonomic_group + sara_status_adj + sqrt_range_std + year_std + cc_threat_ca_adj,
               family = binomial,
               data = m2_supp_db)

summary(m2_supp_std)

# write csv file with model results
m2_supp_results <- tidy(m2_supp_std_clean) %>%
  mutate(p.value = round(p.value, digits=3),
         estimate = round(estimate, digits =3),
         std.error = round(std.error, digits =3),
         statistic = round(statistic, digits =3)) %>%
  rename(Variables = term,
         "Standardized effects" = estimate,
         "Standard Error" = std.error,
         Statistic = statistic,
         "p-value" = p.value) %>%
  mutate(model = c("m2 Supplementary: North American SDM-CC inclusion"))
write.csv(m2_supp_results, paste0(path,"/outputs/model_results/model2_alt_rangesize.csv"))

## Check assumptions ##=========================================================

# check collinearity, outliers, binned residuals
check_model(m2_supp_ptl)

# check linearity with logit for continuous covariates
op <- par(mfrow = c(1, 2))
car::crPlot(m2_supp_ptl, variable = "sqrt_range_std", id = TRUE)
car::crPlot(m2_supp_ptl, variable = "year_std", id = TRUE)

## Clean environment ##=========================================================

rm(m2_supp_db, m2_supp_db_std, m2_supp_db_std_clean, m2_supp_results, m2_supp_ptl, m2_supp_std, m2_supp_std_clean)

# Model 3: probability of inclusion in Canada-inclusive SDM #===================

## Notes: Binomial distribution, logit link, excludes Amphibians, log root transformation of range size

m3_supp_db <- sarsdm %>%
  filter (taxonomic_group != "Amphibians") %>%
  mutate(taxonomic_group = droplevels(taxonomic_group)) %>%
  select (sdm_ca, taxonomic_group, sara_status_adj, year, log_range_am)

## Pairs plot ##================================================================

m3_supp_db %>% ggpairs()

## Standardize data ##==========================================================

m3_supp_db_std <- model.matrix(sdm_ca ~ taxonomic_group + log_range_am + sara_status_adj + year,
                          data = m3_supp_db) %>% #creates a dataframe with dummy binary variables
  as.data.frame() %>%
  mutate_all(scale) %>% #centres a variable around it's mean (i.e., centers it)
  add_column(sdm_ca = m3_supp_db$sdm_ca) %>%
  rename()

m3_supp_db_std_clean <- m3_supp_db_std
colnames(m3_supp_db_std_clean) <- c("Intercept","Birds","Lichens","Mammals","Molluscs","Reptiles","VascularPlants","Rangesize","Threatened","Endangered","Year","sdm_ca")

m3_supp_db <- m3_supp_db %>%
  mutate(log_range_std = log_range_am - mean(log_range_am),
         year_std = year - mean(year))

## Build model ##===============================================================

# for plotting model estimates
m3_supp_std <- glm (sdm_ca ~ 
                 taxonomic_groupBirds +
                 taxonomic_groupLichens +
                 taxonomic_groupMammals +
                 taxonomic_groupMolluscs +
                 taxonomic_groupReptiles +
                 taxonomic_groupVascularPlants +
                 sara_status_adjThreatened + 
                 sara_status_adjEndangered +
                 log_range_am +
                 year,
               family = binomial,
               data = m3_supp_db_std)

# for saving model summary
m3_supp_std_clean <- glm (sdm_ca ~ Birds +Lichens + Mammals + Molluscs +Reptiles +
                       VascularPlants +Threatened +  Endangered +Rangesize +Year,
                     family = binomial,
                     data = m3_supp_db_std_clean)

# for plotting model assumptions, etc.
m3_supp_ptl <- glm (sdm_ca ~ taxonomic_group + sara_status_adj + log_range_std + year_std,
               family = binomial,
               data = m3_supp_db)

summary(m3_supp_std)

# write csv file with model results
m3_supp_results <- tidy(m3_supp_std_clean) %>%
  mutate(p.value = round(p.value, digits=3),
         estimate = round(estimate, digits =3),
         std.error = round(std.error, digits =3),
         statistic = round(statistic, digits =3)) %>%
  rename(Variables = term,
         "Standardized effects" = estimate,
         "Standard Error" = std.error,
         Statistic = statistic,
         "p-value" = p.value) %>%
  mutate(model = c("m3 Supplementary: Canada-inclusive SDM inclusion"))
write.csv(m3_supp_results, paste0(path,"/outputs/model_results/model3_alt_rangesize.csv"))

## Check assumptions ##=========================================================

# check collinearity, outliers, binned residuals
check_model(m3_supp_ptl)

# check linearity with logit for continuous covariates
op <- par(mfrow = c(1, 2))
car::crPlot(m3_supp_ptl, variable = "log_range_std", id = TRUE)
car::crPlot(m3_supp_ptl, variable = "year_std", id = TRUE)

## Clean environment ##=========================================================

rm(m3_supp_db, m3_supp_db_std, m3_supp_db_std_clean, m3_supp_results, m3_supp_ptl, m3_supp_std, m3_supp_std_clean)

# Model 4: probability of inclusion in Canada-inclusive CC-SDM #================

## Notes: Binomial distribution, logit link, log root transformation of range size

m4_supp_db <- sarsdm %>%
  select(cc_sdm_ca, taxonomic_group, sara_status_adj, log_range_am, year, cc_threat_ca_adj)

## Pairs plot ##================================================================

m4_supp_db %>% ggpairs()

## Standardize data ##==========================================================

m4_supp_db <- m4_supp_db %>%
  mutate(log_range_std = (log_range_am - mean(log_range_am))/sd(log_range_am),
         year_std = (year - mean(year))/sd(year))


m4_supp_db_std <- model.matrix(cc_sdm_ca ~ taxonomic_group + log_range_am + sara_status_adj + year + cc_threat_ca_adj,
                          data = m4_supp_db) %>% #creates a dataframe with dummy binary variables
  as.data.frame() %>%
  mutate_all(scale) %>% #centres a variable around it's mean (i.e., centers it)
  add_column(cc_sdm_ca = m4_supp_db$cc_sdm_ca) %>%
  rename()

m4_supp_db_std_clean <- m4_supp_db_std
colnames(m4_supp_db_std_clean) <- c("Intercept","Amphibians","Birds","Lichens","Mammals","Molluscs","Reptiles","VascularPlants","Rangesize","Threatened","Endangered","Year","CCThreat","cc_sdm_ca")

## Build model ##===============================================================

# for plotting model estimates
m4_supp_std <- glm (cc_sdm_ca ~ 
                 taxonomic_groupAmphibians +
                 taxonomic_groupBirds +
                 taxonomic_groupLichens +
                 taxonomic_groupMammals +
                 taxonomic_groupMolluscs +
                 taxonomic_groupReptiles +
                 taxonomic_groupVascularPlants +
                 sara_status_adjThreatened + 
                 sara_status_adjEndangered +
                 log_range_am +
                 cc_threat_ca_adj1 +
                 year,
               family = binomial,
               data = m4_supp_db_std)

# for saving model estimates
m4_supp_std_clean <- glm (cc_sdm_ca ~ Amphibians + Birds +Lichens +Mammals + Molluscs +
                       Reptiles +VascularPlants +Threatened + Endangered +
                       Rangesize +CCThreat +Year,
                     family = binomial,
                     data = m4_supp_db_std_clean)

# model where some of the covariates are standardized - will use this for plotting
m4_supp_ptl <- glm (cc_sdm_ca ~ taxonomic_group + sara_status_adj + log_range_std + year_std + cc_threat_ca_adj,
               family = binomial,
               data = m4_supp_db)

summary(m4_supp_std)

# write csv file with model results
m4_supp_results <- tidy(m4_supp_std_clean) %>%
  mutate(p.value = round(p.value, digits=3),
         estimate = round(estimate, digits =3),
         std.error = round(std.error, digits =3),
         statistic = round(statistic, digits =3)) %>%
  rename(Variables = term,
         "Standardized effects" = estimate,
         "Standard Error" = std.error,
         Statistic = statistic,
         "p-value" = p.value) %>%
  mutate(model = c("m4 Supplementary: Canada-inclusive SDM-CC inclusion"))
write.csv(m4_supp_results, paste0(path,"/outputs/model_results/model4_alt_rangesize.csv"))

## Check assumptions ##=========================================================

# check collinearity, outliers, binned residuals
check_model(m4_supp_ptl)

# check linearity with logit for continuous covariates
op <- par(mfrow = c(1, 2))
car::crPlot(m4_supp_ptl, variable = "log_range_std", id = TRUE)
car::crPlot(m4_supp_ptl, variable = "year_std", id = TRUE)

## Clean environment ##=========================================================

rm(m4_supp_db, m4_supp_db_std, m4_supp_db_std_clean, m4_supp_results, m4_supp_ptl, m4_supp_std, m4_supp_std_clean)

# Model 5: Canada-inclusive SDM research effort #===============================

# Notes: Poisson distribution, log link, log transformed range size, excludes Birds

m5_supp_db <- sarsdm %>%
  filter (taxonomic_group != "Birds") %>%
  mutate(taxonomic_group = droplevels(taxonomic_group)) %>%
  select(n_sdm_ca, taxonomic_group, sara_status_adj, log_range_am, year)

## Pairs plot ##================================================================

m5_supp_db %>% ggpairs()

## Standardize data ##==========================================================

m5_supp_db <- m5_supp_db %>%
  mutate(log_range_std = (log_range_am - mean(log_range_am))/sd(log_range_am),
         year_std = (year - mean(year))/sd(year))

m5_supp_db_std <- model.matrix(n_sdm_ca ~ taxonomic_group + log_range_am + sara_status_adj + year,
                          data = m5_supp_db) %>% #creates a dataframe with dummy binary variables
  as.data.frame() %>%
  mutate_all(scale) %>% #centres a variable around it's mean (i.e., centers it)
  add_column(n_sdm_ca = m5_supp_db$n_sdm_ca) %>%
  rename()

m5_supp_db_std_clean <- m5_supp_db_std
colnames(m5_supp_db_std_clean) <- c("Intercept","Amphibians","Lichens","Mammals","Molluscs","Reptiles","VascularPlants","Rangesize","Threatened","Endangered","Year","n_sdm_ca")

## Build model ##===============================================================

# for plotting model estimates
m5_supp_std <- glm (n_sdm_ca ~ 
                 taxonomic_groupAmphibians +
                 taxonomic_groupLichens +
                 taxonomic_groupMammals +
                 taxonomic_groupMolluscs +
                 taxonomic_groupReptiles +
                 taxonomic_groupVascularPlants +
                 sara_status_adjThreatened + 
                 sara_status_adjEndangered +
                 log_range_am +
                 year,
               family = poisson,
               data = m5_supp_db_std)

# for saving model estimates
m5_supp_std_clean <- glm (n_sdm_ca ~  Amphibians +Lichens + Mammals +Molluscs +
                       Reptiles +VascularPlants +Threatened + Endangered + Rangesize +
                       Year,
                     family = poisson,
                     data = m5_supp_db_std_clean)

# model where some of the covariates are standardized - will use this for plotting
m5_supp_ptl <- glm (n_sdm_ca ~
                 taxonomic_group +
                 sara_status_adj +
                 log_range_std +
                 year_std,
               family = poisson,
               data = m5_supp_db)

summary(m5_supp_std)

# write csv file with model results
m5_supp_results <- tidy(m5_supp_std_clean) %>%
  mutate(p.value = round(p.value, digits=3),
         estimate = round(estimate, digits =3),
         std.error = round(std.error, digits =3),
         statistic = round(statistic, digits =3)) %>%
  rename(Variables = term,
         "Standardized effects" = estimate,
         "Standard Error" = std.error,
         Statistic = statistic,
         "p-value" = p.value) %>%
  mutate(model = c("m5 Supplementary: Canada-inclusive SDM research effort"))
write.csv(m5_supp_results, paste0(path,"/outputs/model_results/model5_alt_rangesize.csv"))

## Check assumptions ##=========================================================

# posterior predictive check, collinearity, outliers
check_model(m5_supp_ptl)

# check for overdispersion
check_overdispersion(m5_supp_ptl)

# residuals vs fitted
op <- par(mfrow=c(1,1),mar=c(5,4,1,2)) #graphing window with white space
plot(m5_supp_ptl, which = 1) 

# linearity with log for continuous covariates
op <- par(mfrow=c(1,2)) #graphing window with white space between panels
car::crPlot(m5_supp_ptl, variable = "log_range_std", id = TRUE)
car::crPlot(m5_supp_ptl, variable = "year_std", id = TRUE)

## Clean environment ##=========================================================

rm(m5_supp_db, m5_supp_db_std, m5_supp_db_std_clean, m5_supp_results, m5_supp_ptl, m5_supp_std, m5_supp_std_clean)

# Model 6: Canada-inclusive CC-SDM research effort #============================

# Notes: Poisson distribution, log link, log transformed range size, excludes Birds

m6_supp_db <- sarsdm %>%
  filter(taxonomic_group != "Birds") %>%
  mutate(taxonomic_group = droplevels(taxonomic_group)) %>%
  select (n_cc_sdm_ca, taxonomic_group, sara_status_adj, log_range_am, year, cc_threat_ca_adj)

## Pairs plot ##================================================================

m6_supp_db %>% ggpairs()

## Standardize data ##==========================================================

m6_supp_db <- m6_supp_db %>%
  mutate(log_range_std = (log_range_am - mean(log_range_am))/sd(log_range_am),
         year_std = (year - mean(year))/sd(year))

m6_supp_db_std <- model.matrix(n_cc_sdm_ca ~ taxonomic_group + log_range_am + sara_status_adj + year + cc_threat_ca_adj,
                          data = m6_supp_db) %>% #creates a dataframe with dummy binary variables
  as.data.frame() %>%
  mutate_all(scale) %>% #centres a variable around it's mean (i.e., centers it)
  add_column(n_cc_sdm_ca = m6_supp_db$n_cc_sdm_ca) %>%
  rename()

m6_supp_db_std_clean <- m6_supp_db_std
colnames(m6_supp_db_std_clean) <- c("Intercept","Amphibians","Lichens","Mammals","Molluscs","Reptiles","VascularPlants","Rangesize","Threatened","Endangered","Year","CCThreat","n_cc_sdm_ca")

## Build model ##===============================================================

# for plotting model estimates
m6_supp_std <- glm (n_cc_sdm_ca ~ 
                 taxonomic_groupAmphibians +
                 taxonomic_groupLichens +
                 taxonomic_groupMammals +
                 taxonomic_groupMolluscs +
                 taxonomic_groupReptiles +
                 taxonomic_groupVascularPlants +
                 sara_status_adjThreatened + 
                 sara_status_adjEndangered +
                 log_range_am +
                 year +
                 cc_threat_ca_adj1,
               family = poisson,
               data = m6_supp_db_std)

# for saving model estimates
m6_supp_std_clean <- glm (n_cc_sdm_ca ~ Amphibians +Lichens +Mammals + Molluscs +
                        Reptiles +VascularPlants +Threatened + Endangered +
                       Rangesize +Year +CCThreat,
                     family = poisson,
                     data = m6_supp_db_std_clean)

# for plotting model assumptions, etc
m6_supp_ptl <- glm (n_cc_sdm_ca ~
                 taxonomic_group +
                 sara_status_adj +
                 log_range_std +
                 year_std +
                 cc_threat_ca_adj,
               family = poisson,
               data = m6_supp_db)

summary(m6_supp_std)

# write csv file with model results
m6_supp_results <- tidy(m6_supp_std_clean) %>%
  mutate(p.value = round(p.value, digits=3),
         estimate = round(estimate, digits =3),
         std.error = round(std.error, digits =3),
         statistic = round(statistic, digits =3)) %>%
  rename(Variables = term,
         "Standardized effects" = estimate,
         "Standard Error" = std.error,
         Statistic = statistic,
         "p-value" = p.value) %>%
  mutate(model = c("m6 Supplementary: Canada-inclusive SDM-CC research effort"))
write.csv(m6_supp_results, paste0(path,"/outputs/model_results/model6_alt_rangesize.csv"))

## Check assumptions ##=========================================================

# posterior predictive check, collinearity, outliers
check_model(m6_supp_ptl)

# check for overdispersion
check_overdispersion(m6_supp_ptl)

# residuals vs fitted
op <- par(mfrow=c(1,1),mar=c(5,4,1,2)) #graphing window with white space
plot(m6_supp_ptl, which = 1) 

# linearity with log for continuous covariates
op <- par(mfrow=c(1,2)) #graphing window with white space between panels
car::crPlot(m6_supp_ptl, variable = "log_range_std", id = TRUE)
car::crPlot(m6_supp_ptl, variable = "year_std", id = TRUE)

## Clean environment ##=========================================================

rm(m6_supp_db, m6_supp_db_std, m6_supp_db_std_clean, m6_supp_results, m6_supp_ptl, m6_supp_std, m6_supp_std_clean)

# Model 7: North American SDM research effort #=================================

# Notes: Quasipoisson distribution, log link, log transformed range size, excludes Birds
# Poisson distribution was overdispersed
# NegBin distribution was underdispersed

m7_supp_db <- sarsdm %>% 
  filter(taxonomic_group != "Birds") %>%
  mutate(taxonomic_group = droplevels(taxonomic_group)) %>%
  select (n_sdm, taxonomic_group, sara_status_adj, log_range_am, year)

## Pairs plot ##================================================================

m7_supp_db %>% ggpairs()

## Standardize data ##==========================================================

m7_supp_db <- m7_supp_db %>%
  mutate(log_range_std = (log_range_am - mean(log_range_am))/sd(log_range_am),
         year_std = (year - mean(year))/sd(year))

m7_supp_db_std <- model.matrix(n_sdm ~ taxonomic_group + log_range_am + sara_status_adj + year,
                          data = m7_supp_db) %>% #creates a dataframe with dummy binary variables
  as.data.frame() %>%
  mutate_all(scale) %>% #centres a variable around it's mean (i.e., centers it)
  add_column(n_sdm = m7_supp_db$n_sdm) %>%
  rename()


m7_supp_db_std_clean <- m7_supp_db_std
colnames(m7_supp_db_std_clean) <- c("Intercept","Amphibians","Lichens","Mammals","Molluscs","Reptiles","VascularPlants","Rangesize","Threatened","Endangered","Year","n_sdm")

## Build model ##===============================================================

# for plotting model estimates
m7_supp_std <- glm (n_sdm ~ 
                 taxonomic_groupAmphibians +
                 taxonomic_groupLichens +
                 taxonomic_groupMammals +
                 taxonomic_groupMolluscs +
                 taxonomic_groupReptiles +
                 taxonomic_groupVascularPlants +
                 sara_status_adjThreatened + 
                 sara_status_adjEndangered +
                 log_range_am +
                 year,
               family = quasipoisson,
               data = m7_supp_db_std)

# for saving model summary
m7_supp_std_clean <- glm (n_sdm ~Amphibians +Lichens + Mammals +Molluscs +
                       Reptiles +VascularPlants +Threatened + Endangered +Rangesize +
                       Year,
                     family = quasipoisson,
                     data = m7_supp_db_std_clean)

# for plotting model assumptions, etc
m7_supp_ptl <- glm (n_sdm ~
                 taxonomic_group +
                 sara_status_adj +
                 log_range_std +
                 year_std,
               family = quasipoisson,
               data = m7_supp_db)

summary(m7_supp_std)

# write csv file with model results
m7_supp_results <- tidy(m7_supp_std_clean) %>%
  mutate(p.value = round(p.value, digits=3),
         estimate = round(estimate, digits =3),
         std.error = round(std.error, digits =3),
         statistic = round(statistic, digits =3)) %>%
  rename(Variables = term,
         "Standardized effects" = estimate,
         "Standard Error" = std.error,
         Statistic = statistic,
         "p-value" = p.value) %>%
  mutate(model = c("m7 Supplementary: North American SDM research effort"))
write.csv(m7_supp_results, paste0(path,"/outputs/model_results/model7_alt_rangesize.csv"))

## Check assumptions ##=========================================================

# collinearity, outliers, homogeneity of variance, misspecified dispersion / zero-inflation
check_model(m7_supp_ptl)

# Dispersion ratio - Quasipoisson can handle some overdispersion (ratio between 1-2 is good)
dispersion <- sum(resid(m7_supp_ptl, type = "deviance")^2) / df.residual(m7_supp_ptl)
dispersion

# Residuals vs fitted
op <- par(mfrow=c(1,1),mar=c(5,4,1,2)) #graphing window with white space
plot(m7_supp_ptl, which = 1) 

# Linearity with log
op <- par(mfrow=c(1,2)) #graphing window with white space between panels
car::crPlot(m7_supp_ptl, variable = "log_range_std", id = TRUE)
car::crPlot(m7_supp_ptl, variable = "year_std", id = TRUE)

## Clean environment ##=========================================================

rm(m7_supp_db, m7_supp_db_std, m7_supp_db_std_clean, m7_supp_results, m7_supp_ptl, m7_supp_std, m7_supp_std_clean)

# Model 8: North American CC-SDM research effort #==============================

# Notes: Quasipoisson distribution, log link, log transformed range size, excludes Birds
# Poisson distribution was overdispersed
# NegBin distribution was underdispersed

m8_supp_db <- sarsdm %>% 
  filter(taxonomic_group != "Birds") %>%
  mutate(taxonomic_group = droplevels(taxonomic_group)) %>%
  select(n_cc_sdm, taxonomic_group, sara_status_adj, log_range_am, year, cc_threat_ca_adj)

## Pairs plot ##================================================================

m8_supp_db %>% ggpairs()

## Standardize data ##==========================================================

m8_supp_db <- m8_supp_db %>%
  mutate(log_range_std = (log_range_am - mean(log_range_am))/sd(log_range_am),
         year_std = (year - mean(year))/sd(year))

m8_supp_db_std <- model.matrix(n_cc_sdm ~ taxonomic_group + log_range_am + sara_status_adj + year + cc_threat_ca_adj,
                          data = m8_supp_db) %>% #creates a dataframe with dummy binary variables
  as.data.frame() %>%
  mutate_all(scale) %>% #centres a variable around it's mean (i.e., centers it)
  add_column(n_cc_sdm = m8_supp_db$n_cc_sdm) %>%
  rename()

m8_supp_db_std_clean <- m8_supp_db_std
colnames(m8_supp_db_std_clean) <- c("Intercept","Amphibians","Lichens","Mammals","Molluscs","Reptiles","VascularPlants","Rangesize","Threatened","Endangered","Year","CCThreat","n_cc_sdm")

## Build model ##===============================================================

# for plotting model results
m8_supp_std <- glm (n_cc_sdm ~ 
                 taxonomic_groupAmphibians +
                 taxonomic_groupLichens +
                 taxonomic_groupMammals +
                 taxonomic_groupMolluscs +
                 taxonomic_groupReptiles +
                 taxonomic_groupVascularPlants +
                 sara_status_adjThreatened + 
                 sara_status_adjEndangered +
                 log_range_am +
                 year +
                 cc_threat_ca_adj1,
               family = quasipoisson,
               data = m8_supp_db_std)

# for saving model results
m8_supp_std_clean <- glm (n_cc_sdm ~ Amphibians + Lichens +Mammals +Molluscs +
                       Reptiles +VascularPlants +Threatened + Endangered +
                       Rangesize +Year + CCThreat,
                     family = quasipoisson,
                     data = m8_supp_db_std_clean)

# for plotting model assumptions, etc.
m8_supp_ptl <- glm (n_cc_sdm ~
                 taxonomic_group +
                 sara_status_adj +
                 log_range_std +
                 year_std +
                 cc_threat_ca_adj,
               family = quasipoisson,
               data = m8_supp_db)

summary(m8_supp_std)

# write csv file with model results
m8_supp_results <- tidy(m8_supp_std_clean) %>%
  mutate(p.value = round(p.value, digits=3),
         estimate = round(estimate, digits =3),
         std.error = round(std.error, digits =3),
         statistic = round(statistic, digits =3)) %>%
  rename(Variables = term,
         "Standardized effects" = estimate,
         "Standard Error" = std.error,
         Statistic = statistic,
         "p-value" = p.value) %>%
  mutate(model = c("m8 Supplementary: North American SDM-CC research effort"))
write.csv(m8_supp_results, paste0(path,"/outputs/model_results/model8_alt_rangesize.csv"))

## Check assumptions ##=========================================================

# collinearity, outliers, homogeneity of variance, misspecified dispersion / zero-inflation
check_model(m8_supp_ptl)

# Dispersion ratio - Quasipoisson can handle some overdispersion (ratio between 1-2 is good)
dispersion <- sum(resid(m8_supp_ptl, type = "deviance")^2) / df.residual(m8_supp_ptl)
dispersion

# Residuals vs fitted
op <- par(mfrow=c(1,1),mar=c(5,4,1,2)) #graphing window with white space
plot(m8_supp_ptl, which = 1) 

# Linearity with log
op <- par(mfrow=c(1,2)) #graphing window with white space between panels
car::crPlot(m8_supp_ptl, variable = "log_range_std", id = TRUE)
car::crPlot(m8_supp_ptl, variable = "year_std", id = TRUE)

## Clean environment ##=========================================================

rm(m8_supp_db, m8_supp_db_std, m8_supp_db_std_clean, m8_supp_results, m8_supp_ptl, m8_supp_std, m8_supp_std_clean)

