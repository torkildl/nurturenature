library(huxtable)
library(here)
#library(tidyverse)
library(ggplot2)
library(tidyr)
library(stringr)
library(broom)
library(purrr)
library(data.table)
library(broom)
library(lme4)
library(broom.mixed)
library(tibble)
library(multcomp)
library(dplyr)

###
rm(list=ls())
source(here("../my-useful-R-functions.R"))

###
### Step 3. Conventional and parent-sib models per subject/grade
###
complete_trios <- fread(here("noncognitive-cognitive/final/data/complete-trios.csv")) %>% as.data.frame
sibtrios <- fread(here("noncognitive-cognitive/final/data/sibtrios.csv")) %>% as.data.frame

###
### Defintions
### 
### Outcomes: Three subjects and three measurement points, but no English test in 9th grade.
outcomes <- c("NPREG05","NPREG08","NPREG09","NPLES05","NPLES08","NPLES09","NPENG05","NPENG08")
names(outcomes) <- c("Math5th","Math8th","Math9th","Reading5th","Reading8th","Reading9th","English5th","English8th")

### Control variables definition
controls <- c("child_female",
              "as.factor(child_cohort)",str_c("child_PC",1:10))

### Estimation function
estimate_ols_model <- function(outcome,model_def) {
  m <- lm(data=analysis_df, formula=as.formula(str_c(outcome,"~", model_def)))
  return(m)
}

###
### ORDINARY OLS MODELS, SEPARATE SUBJECTS
###
analysis_df <- complete_trios %>% as.data.frame %>% 
  # Scale the PGS with N(0,1)
  filter(!is.na(child_EA4)) %>% 
  filter(!is.na(mother_EA4)) %>% 
  filter(!is.na(father_EA4)) %>% 
  mutate_at(vars(one_of(outcomes)), as.numeric)
analysis_df$avgtest <- rowMeans(dplyr::select(analysis_df,all_of(outcomes)),na.rm = TRUE)

# Cognitive/noncognitive
child_genes <- c("child_Cog","child_NonCog")
genetic_nurture <- c("father_Cog","mother_Cog",
                     "father_NonCog","mother_NonCog")
cog_model_defs <- list(
  child_model = str_c(paste(child_genes, collapse=" + ")," + ", 
                      paste(controls, collapse = " + ")),
  nurture_model = str_c(paste(child_genes, collapse=" + "), " + ", 
                        paste(genetic_nurture, collapse = " + ")," +",
                        paste(controls, collapse = " + "))
)
### Estimate regression models, cog/noncog
cog_modnames <- str_c(rep(names(outcomes), times=length(cog_model_defs)), ":", rep(names(cog_model_defs), each=length(outcomes)))
cog_models <- map(cog_model_defs, function(x) map(outcomes, estimate_ols_model, x)) %>% 
  flatten %>% 
  setNames(nm = str_c("Std:Cog/NonCog:",cog_modnames))

# EA4
child_genes <- c("child_EA4")
genetic_nurture <- c("father_EA4","mother_EA4")

EA4_model_defs <- list(
  child_model = str_c(paste(child_genes, collapse=" + ")," + ", 
                      paste(controls, collapse = " + ")),
  nurture_model = str_c(paste(child_genes, collapse=" + "), " + ", 
                        paste(genetic_nurture, collapse = " + ")," +",
                        paste(controls, collapse = " + "))
)
### Estimate regression models, EA4
EA4_modnames <- str_c(rep(names(outcomes), times=length(EA4_model_defs)), ":", rep(names(EA4_model_defs), each=length(outcomes)))
EA4_models <- map(EA4_model_defs, function(x) map(outcomes, estimate_ols_model, x)) %>% 
  flatten %>% 
  setNames(nm = str_c("Std:EA4:",EA4_modnames))


# WFEA
child_genes <- c("child_WFEA")
genetic_nurture <- c("father_WFEA","mother_WFEA")

WFEA_model_defs <- list(
  child_model = str_c(paste(child_genes, collapse=" + ")," + ", 
                      paste(controls, collapse = " + ")),
  nurture_model = str_c(paste(child_genes, collapse=" + "), " + ", 
                        paste(genetic_nurture, collapse = " + ")," +",
                        paste(controls, collapse = " + "))
)
### Estimate regression models, WFEA
WFEA_modnames <- str_c(rep(names(outcomes), times=length(WFEA_model_defs)), ":", rep(names(WFEA_model_defs), each=length(outcomes)))
WFEA_models <- map(WFEA_model_defs, function(x) map(outcomes, estimate_ols_model, x)) %>% 
  flatten %>% 
  setNames(nm = str_c("Std:WFEA:",WFEA_modnames))

ols_separate_models <- c(cog_models, EA4_models, WFEA_models)

###
### OLS separate-outcome analysis with parent-sibs
###
analysis_df <- sibtrios %>% as.data.frame %>% 
  # Scale the PGS with N(0,1)
  mutate_at(vars(contains("BATCH")), as.factor) %>% 
  mutate_at(vars(one_of(outcomes)), as.numeric) %>% 
  filter(!is.na(child_EA4)) %>% 
  filter(!is.na(sibparent_dev_EA4)) %>% 
  filter(!is.na(sibparent_mean_EA4)) %>% 
  filter(!is.na(otherparent_EA4))
analysis_df$avgtest <- rowMeans(select(analysis_df,all_of(outcomes)),na.rm = TRUE)


### Control variables definitions
controls <- c("child_female", "as.factor(child_cohort)",  str_c("child_PC",1:10))
parent_controls <- c()

### Helper funtction to write out the model
writemodel <- function(...) { 
  vars <- as.character(as.list(match.call(expand.dots = TRUE))) %>% .[2:length(.)]
  nvars <- map(vars, get) %>% 
    map(function(x) str_c(paste(x, collapse=" + ")))
  amodel <- str_c(paste(nvars, collapse = " + ")) 
  return(amodel)
}

### Separate models per outcome, sib-parent analysis

# EA4
pgses <- c("EA4")
child_pgs <- str_c("child_",pgses)
child_model <- writemodel(child_pgs, controls)
nurture_pgs <- crossing(c("mother","father"), pgses) %>% unite(col = "v") %>% pull
nurture_model <- writemodel(child_pgs, controls, nurture_pgs)
sibling_pgs <- crossing(c("sibparent_mean","sibparent_dev","otherparent"), pgses) %>% unite(col = "v") %>% pull
sibling_model <- writemodel(child_pgs, controls, sibling_pgs)
uncle_pgs <- crossing(c("sibparent","othersibling","otherparent"), pgses) %>% unite(col = "v") %>% pull
uncle_model <- writemodel(child_pgs, controls, uncle_pgs)
sib_EA4_model_defs <- list("child_model" = child_model,
                           "nurture_model" = nurture_model,
                           "sibling_model" = sibling_model,
                           "uncle_model" = uncle_model)
sib_EA4_modnames <- str_c(rep(names(outcomes), times=length(sib_EA4_model_defs)), ":", rep(names(sib_EA4_model_defs), each=length(outcomes)))
sib_EA4_models <- map(sib_EA4_model_defs, function(x) map(outcomes, estimate_ols_model, x)) %>% 
  flatten %>% 
  setNames(nm = str_c("Sib:EA4:",sib_EA4_modnames))

# Cog/NonCog
pgses <- c("Cog","NonCog")
child_pgs <- str_c("child_",pgses)
child_model <- writemodel(child_pgs, controls)
nurture_pgs <- crossing(c("mother","father"), pgses) %>% unite(col = "v") %>% pull
nurture_model <- writemodel(child_pgs, controls, nurture_pgs)
sibling_pgs <- crossing(c("sibparent_mean","sibparent_dev","otherparent"), pgses) %>% unite(col = "v") %>% pull
sibling_model <- writemodel(child_pgs, controls, sibling_pgs)
uncle_pgs <- crossing(c("sibparent","othersibling","otherparent"), pgses) %>% unite(col = "v") %>% pull
uncle_model <- writemodel(child_pgs, controls, uncle_pgs)
sib_cog_model_defs <- list("child_model" = child_model,
                           "nurture_model" = nurture_model,
                   "sibling_model" = sibling_model,
                   "uncle_model" = uncle_model)
sib_cog_modnames <- str_c(rep(names(outcomes), times=length(sib_cog_model_defs)), ":", rep(names(sib_cog_model_defs), each=length(outcomes)))
sib_cog_models <- map(sib_cog_model_defs, function(x) map(outcomes, estimate_ols_model, x)) %>% 
  flatten %>% 
  setNames(nm = str_c("Sib:Cog/NonCog:",sib_cog_modnames))

# WFEA
pgses <- c("WFEA")
child_pgs <- str_c("child_",pgses)
child_model <- writemodel(child_pgs, controls)
nurture_pgs <- crossing(c("mother","father"), pgses) %>% unite(col = "v") %>% pull
nurture_model <- writemodel(child_pgs, controls, nurture_pgs)
sibling_pgs <- crossing(c("sibparent_mean","sibparent_dev","otherparent"), pgses) %>% unite(col = "v") %>% pull
sibling_model <- writemodel(child_pgs, controls, sibling_pgs)
uncle_pgs <- crossing(c("sibparent","othersibling","otherparent"), pgses) %>% unite(col = "v") %>% pull
uncle_model <- writemodel(child_pgs, controls, uncle_pgs)
sib_WFEA_model_defs <- list("child_model" = child_model,
                            "nurture_model" = nurture_model,
                           "sibling_model" = sibling_model,
                           "uncle_model" = uncle_model)
sib_WFEA_modnames <- str_c(rep(names(outcomes), times=length(sib_WFEA_model_defs)), ":", rep(names(sib_WFEA_model_defs), each=length(outcomes)))
sib_WFEA_models <- map(sib_WFEA_model_defs, function(x) map(outcomes, estimate_ols_model, x)) %>% 
  flatten %>% 
  setNames(nm = str_c("Sib:WFEA:",sib_WFEA_modnames))

sib_ols_separate_models <- c(sib_EA4_models,sib_cog_models,sib_WFEA_models)

all_ols_models <- c(ols_separate_models,sib_ols_separate_models)


###
### Random effects models
###
re_controls <- c(controls,"testname")
outcomes_nonames <- outcomes %>% setNames(nm=outcomes)
stdlong <- complete_trios %>% 
  select(any_of(outcomes_nonames), child_cohort, child_female, child_w19_0634_lnr,
         child_EA4, mother_EA4, father_EA4,
         child_WFEA, mother_WFEA, father_WFEA,
         child_Cog, mother_Cog, father_Cog,
         child_NonCog, mother_NonCog, father_NonCog,
         any_of(controls)) %>% 
  gather(testname, testscore, any_of(outcomes)) %>% 
  dplyr::select(-starts_with("NP")) %>% 
  filter(!is.na(testscore)) %>%
  mutate(testscore = as.numeric(testscore))
stdavg <- stdlong %>% 
  group_by(child_w19_0634_lnr) %>% 
  mutate(avgtest = mean(testscore)) %>% slice(1) %>% 
  ungroup

siblong <- sibtrios %>% 
  gather(testname, testscore, any_of(outcomes)) %>% 
  dplyr::select(-starts_with("NP")) %>% 
  filter(!is.na(testscore)) %>%
  mutate(testscore = as.numeric(testscore))
sibavg <- siblong %>% 
  group_by(child_w19_0634_lnr) %>% 
  mutate(avgtest = mean(testscore)) %>% slice(1) %>% 
  ungroup


###
### Complete trio sample RE models
###

# Child models
gen_std_child_model <- function(PGI) {
  pgis <- unlist(map(PGI, ~str_c(c("child_"),.x)))
  aformula <- str_c("testscore ~ ",
                    str_c(pgis,collapse = " + ")," + ",
                    str_c(re_controls,collapse=" + "),
                    " + (1 | child_w19_0634_lnr) ")
  return(aformula)
}
std_ea4_child_model <- lmer(data=stdlong, formula = gen_std_child_model("EA4"))
std_cog_child_model <- lmer(data=stdlong, formula = gen_std_child_model(c("Cog","NonCog")))
std_wfea_child_model <- lmer(data=stdlong, formula = gen_std_child_model("WFEA"))

# Nurture models
gen_std_nurture_model <- function(PGI) {
  pgis <- unlist(map(PGI, ~str_c(c("child_","mother_","father_"),.x)))
  aformula <- str_c("testscore ~ ",
                    str_c(pgis,collapse = " + ")," + ",
                    str_c(re_controls,collapse=" + "),
                    " + (1 | child_w19_0634_lnr) ")
  return(aformula)
}
std_ea4_nurture_model <- lmer(data=stdlong, formula = gen_std_nurture_model("EA4"))
std_cog_nurture_model <- lmer(data=stdlong, formula = gen_std_nurture_model(c("Cog","NonCog")))
std_wfea_nurture_model <- lmer(data=stdlong, formula = gen_std_nurture_model("WFEA"))

std_re_child_models <- c(
  "Std:EA4:Combined:child_model" = std_ea4_child_model,
  "Std:Cog:Combined:child_model" = std_cog_child_model,
  "Std:WFEA:Combined:child_model" = std_wfea_child_model
)
std_re_nurture_models <- c(
  "Std:EA4:Combined:nurture_model" = std_ea4_nurture_model,
  "Std:Cog:Combined:nurture_model" = std_cog_nurture_model,
  "Std:WFEA:Combined:nurture_model" = std_wfea_nurture_model
)
std_re_models <- c(std_re_child_models,
                   std_re_nurture_models)

###
### Parent-sib RE models
###

# Child models (1)
gen_sib_child_model <- function(PGI) {
  pgis <- unlist(map(PGI, ~str_c(c("child_"),.x)))
  aformula <- str_c("testscore ~ ",
                    str_c(pgis,collapse = " + ")," + ",
                    str_c(re_controls,collapse=" + "),
                    " + (1 | child_w19_0634_lnr) ")
  return(aformula)
}
sib_ea4_child_model <- lmer(data=siblong, formula = gen_sib_child_model("EA4"))
sib_cog_child_model <- lmer(data=siblong, formula = gen_sib_child_model(c("Cog","NonCog")))
sib_wfea_child_model <- lmer(data=siblong, formula = gen_sib_child_model("WFEA"))

# Nurture models (2)
gen_sib_nurture_model <- function(PGI) {
  pgis <- unlist(map(PGI, ~str_c(c("child_","mother_","father_"),.x)))
  aformula <- str_c("testscore ~ ",
                    str_c(pgis,collapse = " + ")," + ",
                    str_c(re_controls,collapse=" + "),
                    " + (1 | child_w19_0634_lnr) ")
  return(aformula)
}
sib_ea4_nurture_model <- lmer(data=siblong, formula = gen_sib_nurture_model("EA4"))
sib_cog_nurture_model <- lmer(data=siblong, formula = gen_sib_nurture_model(c("Cog","NonCog")))
sib_wfea_nurture_model <- lmer(data=siblong, formula = gen_sib_nurture_model("WFEA"))

# Sibling models (3)
gen_sib_sibling_model <- function(PGI) {
  pgis <- unlist(map(PGI, ~str_c(c("child_","sibparent_mean_","sibparent_dev_","otherparent_"),.x)))
  aformula <- str_c("testscore ~ ",
                  str_c(pgis,collapse = " + ")," + ",
                  str_c(re_controls,collapse=" + "),
                  " + (1 | child_w19_0634_lnr) ")
  return(aformula)
}
sib_ea4_sibling_model <- lmer(data=siblong, formula = gen_sib_sibling_model("EA4"))
sib_cog_sibling_model <- lmer(data=siblong, formula = gen_sib_sibling_model(c("Cog","NonCog")))
sib_wfea_sibling_model <- lmer(data=siblong, formula = gen_sib_sibling_model("WFEA"))

# Uncle/aunt models
gen_sib_uncle_model <- function(PGI) {
  pgis <- unlist(map(PGI, ~str_c(c("child_","sibparent_","otherparent_","othersibling_"),.x)))
  aformula <- str_c("testscore ~ ",
                    str_c(pgis,collapse = " + ")," + ",
                    str_c(re_controls,collapse=" + "),
                    " + (1 | child_w19_0634_lnr) ")
  return(aformula)
}
sib_ea4_uncle_model <- lmer(data=siblong, formula = gen_sib_uncle_model("EA4"))
sib_cog_uncle_model <- lmer(data=siblong, formula = gen_sib_uncle_model(c("Cog","NonCog")))
sib_wfea_uncle_model <- lmer(data=siblong, formula = gen_sib_uncle_model("WFEA"))


sib_re_child_models <- c(
  "Sib:EA4:Combined:child_model" = sib_ea4_child_model,
  "Sib:Cog:Combined:child_model" = sib_cog_child_model,
  "Sib:WFEA:Combined:child_model" = sib_wfea_child_model
)
sib_re_nurture_models <- c(
  "Sib:EA4:Combined:nurture_model" = sib_ea4_nurture_model,
  "Sib:Cog:Combined:nurture_model" = sib_cog_nurture_model,
  "Sib:WFEA:Combined:nurture_model" = sib_wfea_nurture_model
)
sib_re_sibling_models <- c(
  "Sib:EA4:Combined:sibling_model" = sib_ea4_sibling_model,
  "Sib:Cog:Combined:sibling_model" = sib_cog_sibling_model,
  "Sib:WFEA:Combined:sibling_model" = sib_wfea_sibling_model
)
sib_re_uncle_models <- c(
  "Sib:EA4:Combined:uncle_model" = sib_ea4_uncle_model,
  "Sib:Cog:Combined:uncle_model" = sib_cog_uncle_model,
  "Sib:WFEA:Combined:uncle_model" = sib_wfea_uncle_model
)

sib_re_models <- c(sib_re_child_models,
                   sib_re_nurture_models, 
                   sib_re_sibling_models,
                   sib_re_uncle_models)
all_re_models <- c(std_re_models, sib_re_models)


ols_results <- map(all_ols_models, function(x) broom::tidy(x,conf.int=T)) %>% 
  bind_rows(.id="model") %>% 
  separate(col=model, into=c("design","PGIdef","outcome","model_def"), sep=":", remove=T) %>% 
  mutate(method="OLS") 

re_results <- map(all_re_models, function(x) broom.mixed::tidy(x,conf.int=T)) %>% 
  bind_rows(.id="model") %>% 
  separate(col=model, into=c("design","PGIdef","outcome","model_def"), sep=":", remove=T) %>% 
  mutate(method="RE") 

results <- bind_rows(ols_results, re_results)
fwrite(results, here("noncognitive-cognitive/final/output/regression-results-output.txt"))

all_models <- c(all_ols_models, all_re_models)
names(all_models)


###
### Tables and graphs for main text and supplement
###

### Sample descriptive statistics
selected_vars <- c(
  "Child EA4 Score" = "child_EA4",
  "Child Cog Score"= "child_Cog",
  "Child NonCog Score"= "child_NonCog",
  "Child WFEA Score"= "child_WFEA",
  "Mother EA4 Score" = "mother_EA4",
  "Mother Cog Score"= "mother_Cog",
  "Mother NonCog Score"= "mother_NonCog",
  "Mother WFEA Score"= "mother_WFEA",
  "Father EA4 Score" = "father_EA4",
  "Father Cog Score"= "father_Cog",
  "Father NonCog Score"= "father_NonCog",
  "Father WFEA Score"= "father_WFEA",
  "Math 5th" = "NPREG05",
  "Math 8th" = "NPREG08",
  "Math 9th" = "NPREG09",
  "Reading 5th" = "NPLES05",
  "Reading 8th" = "NPLES08",
  "Reading 9th" = "NPLES09",
  "English 5th" = "NPENG05",
  "English 8th" = "NPENG08",
  "Child Year of Birth" = "child_cohort",
  "Mother Year of Birth" = "mother_cohort",
  "Father Year of Birth" = "father_cohort",
  "Child female" = "child_female"
)

getstats <- function(x) {
  x %>% 
    filter(!is.na(child_EA4)) %>% 
    filter(!is.na(mother_EA4)) %>% 
    filter(!is.na(father_EA4)) %>% 
    select(one_of(selected_vars)) %>% 
    psych::describe(x = .) %>% 
    as.data.frame %>% 
    rownames_to_column("rowname") %>% 
    select(rowname, N = n, Mean = mean, SD = sd, Min = min, Max = max) %>% 
    left_join(data.frame(rowname = selected_vars, Variable=names(selected_vars))) %>% 
    select(Variable, N, Mean, SD, Min, Max) %>% 
    mutate_at(vars(-Variable), ~round(.x,digits=3))
}
tab1a_df <- complete_trios %>% getstats %>% mutate(Sample = "Trios")
tab1b_df <- sibtrios %>% getstats %>% mutate(Sample = "Siblings")
table1 <- bind_rows(tab1a_df,tab1b_df) %>% select(Sample, Variable, N, Mean, SD, Min, Max)
xlsx::write.xlsx2(table1, file = here("noncognitive-cognitive/final/output/table1-descriptives.xlsx"), col.names = TRUE, sheetName = "Table 1")



### Figure labels
theylab <- "Association with children's academic achievement (Beta)"
thexlab <- ""

### Figure 1
###
### Plot results for RE models
###
myterms <- c("child_EA4","child_Cog","child_NonCog","child_WFEA",
             "mother_EA4","mother_Cog","mother_NonCog","mother_WFEA",
             "father_EA4","father_Cog","father_NonCog","father_WFEA",
             "sibparent_dev_EA4","sibparent_mean_EA4", "otherparent_EA4",
             "sibparent_dev_Cog","sibparent_mean_Cog", "otherparent_Cog",
             "sibparent_dev_NonCog","sibparent_mean_NonCog", "otherparent_NonCog",
             "sibparent_dev_WFEA","sibparent_mean_WFEA", "otherparent_WFEA",
             "sibparent_EA4", "sibparent_WFEA", "sibparent_Cog", "sibparent_NonCog",
             "othersibling_EA4","othersibling_Cog","othersibling_NonCog","othersibling_WFEA"
)

childterm <- "Child's\n PGI"
motherterm <- "Mother's\n PGI"
fatherterm <- "Father's\n PGI"
devterm <- "Parent's\n Deviation\n from\n Parent-Sibling\n Pair Mean PGI"
meanterm <- "Mean of\n Parent's\n and\n Sibling's\n PGIs"
otherterm <- "Other\n Parent's\n PGI"
parentterm <- "Parent's\n PGI"
uncleterm <- "Parent's\n Sibling's\n PGI"
names(myterms) <- c(rep(childterm,4),
                    rep(motherterm,4),
                    rep(fatherterm,4),
                    c(devterm, meanterm, otherterm),
                    c(devterm, meanterm, otherterm),
                    c(devterm, meanterm, otherterm),
                    c(devterm, meanterm, otherterm),
                    rep(parentterm,4),
                    rep(uncleterm,4)
)
termsdf <- data.frame(term = myterms, label = names(myterms)) %>% 
  mutate(label = factor(label, levels = c(childterm, motherterm, fatherterm, meanterm, devterm, parentterm, otherterm, uncleterm)))

broadplot <- results %>% 
  filter(outcome=="Combined") %>%
  filter(design=="Sib") %>% 
  filter(model_def %in% c("child_model","nurture_model","sibling_model","uncle_model")) %>% 
  mutate(model_def = case_when(model_def=="child_model" ~ "Model 1",
                               model_def=="nurture_model" ~ "Model 2",
                               model_def=="sibling_model" ~ "Model 3",
                               model_def=="uncle_model" ~ "Model 4")) %>% 
  filter(term %in% myterms) %>% 
  left_join(termsdf, "term") %>%
  mutate(PGI = case_when(
    PGIdef=="Cog" ~ "Cognitive/Non-cognitive\n (Cog/NonCog) PGI",
    PGIdef=="EA4" ~ "Population Educational\n attainment 4 (EA4) PGI",
    PGIdef=="WFEA" ~"Within-family \nEducational attainment\n (WFEA) PGI"
    
  )) %>% 
  mutate(PGI = factor(PGI, levels = c("Population Educational\n attainment 4 (EA4) PGI",
                                      "Cognitive/Non-cognitive\n (Cog/NonCog) PGI",
                                      "Within-family \nEducational attainment\n (WFEA) PGI"))) %>% 
  mutate(shapefix = case_when(PGIdef=="EA4" ~ "SAME",
                              PGIdef=="WFEA" ~ "SAME",
                              str_detect(term,"NonCog") ~ "Non-cognitive",
                              str_detect(term,"Cog") ~ "Cognitive")) %>% 
  ggplot(mapping = aes(x=label,y=estimate,ymin=conf.low,ymax=conf.high, color=shapefix, shape=shapefix)) +
  geom_pointrange(position = position_dodge2(width=0.20)) +
  geom_abline(slope = 0, intercept=0) + 
  geom_abline(slope = 0, intercept=0.08, linetype=2) + 
  ylab(theylab) + xlab(thexlab) + 
  scale_color_manual(values=c("Non-cognitive"="orange","Cognitive"="blue","EA4"="black","SAME"="black")) +
  facet_grid(PGI~model_def, shrink = TRUE,scales="free_x",space = "free_x", margins = ) +
  theme_minimal() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme(strip.text.x = element_text(size = 11)) + 
  theme(strip.text.y = element_text(size = 11)) + 
  theme(axis.text=element_text(size=9)) +
  theme(axis.title=element_text(size=11)) + 
  theme(legend.position = "none") + theme(panel.background = element_rect(fill="grey90", color=NA))

ggsave(filename = here("noncognitive-cognitive/final/output/figure1.pdf"), broadplot,
       width = 297,height = 210,units = "mm",  dpi = 320,device = "pdf")
ggsave(filename = here("noncognitive-cognitive/final/output/figure1.png"), broadplot,
       width = 297,height = 210,units = "mm",  dpi = 320,device = "png")


bracketplot <- results %>% 
  filter(outcome=="Combined") %>%
  filter(design=="Sib") %>% 
  filter(model_def %in% c("child_model","nurture_model","sibling_model","uncle_model")) %>% 
  mutate(model_def = case_when(model_def=="child_model" ~ "Model 1",
                               model_def=="nurture_model" ~ "Model 2",
                               model_def=="sibling_model" ~ "Model 3",
                               model_def=="uncle_model" ~ "Model 4")) %>% 
  filter(term %in% myterms) %>% 
  left_join(termsdf, "term") %>%
  mutate(PGI = case_when(
    PGIdef=="Cog" ~ "Cognitive/Non-cognitive\n (Cog/NonCog) PGI",
    PGIdef=="EA4" ~ "Population Educational\n attainment 4 (EA4) PGI",
    PGIdef=="WFEA" ~"Within-family \nEducational attainment\n (WFEA) PGI"
    
  )) %>% 
  mutate(PGI = factor(PGI, levels = c("Population Educational\n attainment 4 (EA4) PGI",
                                      "Cognitive/Non-cognitive\n (Cog/NonCog) PGI",
                                      "Within-family \nEducational attainment\n (WFEA) PGI"))) %>% 
  mutate(shapefix = case_when(PGIdef=="EA4" ~ "SAME",
                              PGIdef=="WFEA" ~ "SAME",
                              str_detect(term,"NonCog") ~ "Non-cognitive",
                              str_detect(term,"Cog") ~ "Cognitive")) %>% 
  ggplot(mapping = aes(x=label,y=estimate,ymin=conf.low,ymax=conf.high, color=shapefix, shape=shapefix)) +
  geom_pointrange(position = position_dodge2(width=0.20)) +
  geom_abline(slope = 0, intercept=0) + 
  geom_abline(slope = 0, intercept=0.08, linetype=2) + 
  ylab(theylab) + xlab(thexlab) + 
  scale_color_manual(values=c("Non-cognitive"="orange","Cognitive"="blue","EA4"="black","SAME"="black")) +
  facet_grid(PGI~model_def, shrink = TRUE,scales="free_x",space = "free_x", margins = ) +
  theme_minimal() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme(strip.text.x = element_text(size = 11)) + 
  theme(strip.text.y = element_text(size = 11)) + 
  theme(axis.text=element_text(size=9)) +
  theme(axis.title=element_text(size=11)) + 
  theme(legend.position = "none") + theme(panel.background = element_rect(fill="grey90", color=NA))

bracketplot_draw <- ggdraw(plot=bracketplot) +
  # EA4
  draw_line(x = c(0.455, 0.555), y=c(0.3+0.6, 0.3+0.6), color="black", size=0.4) +
  draw_line(x = c(0.455, 0.455), y=c(0.3+0.6, 0.285+0.6), color="black", size=0.4) +
  draw_line(x = c(0.555, 0.555), y=c(0.3+0.6, 0.285+0.6), color="black", size=0.4) +
  draw_label("*", x = 0.505, y=0.31+0.6, size = 9) +
  # Cog
  draw_line(x = c(0.455, 0.555), y=c(0.3+0.3, 0.3+0.3), color="black", size=0.4) +
  draw_line(x = c(0.455, 0.455), y=c(0.3+0.3, 0.285+0.3), color="black", size=0.4) +
  draw_line(x = c(0.555, 0.555), y=c(0.3+0.3, 0.285+0.3), color="black", size=0.4) +
  draw_label("**", x = 0.505, y=0.31+0.3, size = 9) +
  # WFEA test
  draw_line(x = c(0.455, 0.555), y=c(0.3, 0.3), color="black", size=0.4) +
  draw_line(x = c(0.455, 0.455), y=c(0.3, 0.285), color="black", size=0.4) +
  draw_line(x = c(0.555, 0.555), y=c(0.3, 0.285), color="black", size=0.4) +
  draw_label("*", x = 0.505, y=0.31, size = 9)
  
ggsave(filename = here("noncognitive-cognitive/final/output/figure1B.pdf"), bracketplot_draw,
       width = 297,height = 210,units = "mm",  dpi = 320,device = "pdf")
ggsave(filename = here("noncognitive-cognitive/final/output/figure1B.png"), bracketplot_draw,
       width = 297,height = 210,units = "mm",  dpi = 320,device = "png")


### Figure S1: Models 1 and 2 by test (and combined outcome)
myterms <- c("child_EA4","child_Cog","child_NonCog","child_WFEA",
             "mother_EA4","mother_Cog","mother_NonCog","mother_WFEA",
             "father_EA4","father_Cog","father_NonCog","father_WFEA")

names(myterms) <- c("Child's\n PGI","Child's\n PGI","Child's\n PGI","Child's\n PGI",
                    "Mother's\n PGI","Mother's\n PGI","Mother's\n PGI","Mother's\n PGI",
                    "Father's\n PGI","Father's\n PGI","Father's\n PGI","Father's\n PGI")

termsdf <- data.frame(term = myterms, label = names(myterms))

figureS1 <- results %>% 
  filter(design=="Std") %>% 
  filter(model_def %in% c("nurture_model")) %>% 
  filter(term %in% myterms) %>% 
  left_join(termsdf, "term") %>%
  mutate(outcome = factor(x = outcome, levels = c("Math5th","Math8th","Math9th",
                                                  "Reading5th","Reading8th","Reading9th",
                                                  "English5th","English8th","Combined"))) %>% 
  mutate(PGI = case_when(
    PGIdef=="Cog" ~ "Cognitive/Non-cognitive (Cog/NonCog)",
    PGIdef=="EA4" ~ "Population Educational \n attainment 4 (EA4)",
    PGIdef=="WFEA" ~"Within-family \nEducational attainment (WFEA)"
  )) %>% 
  mutate(shapefix = case_when(PGIdef=="EA4" ~ "EA4",
                              PGIdef=="WFEA" ~ "WFEA",
                              str_detect(term,"NonCog") ~ "Non-cognitive",
                              str_detect(term,"Cog") ~ "Cognitive")) %>% 
  ggplot(mapping = aes(x=label,y=estimate,ymin=conf.low,ymax=conf.high, color=shapefix, shape=shapefix)) +
  geom_pointrange(position = position_dodge2(width=0.20)) +
  geom_abline(slope = 0, intercept=0) + 
  geom_abline(slope = 0, intercept=0.08, linetype=2) + 
  ylab(theylab) + xlab(thexlab) + labs(shape = "PGI") +
  scale_color_manual(values=c("Non-cognitive"="orange","Cognitive"="blue","EA4"="black","WFEA"="black"), guide = "none") +
  facet_wrap(vars(outcome), ncol = 3) +
  theme_minimal() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(strip.text.x = element_text(size = 11)) + 
  theme(strip.text.y = element_text(size = 11)) + 
  theme(axis.text=element_text(size=9)) +
  theme(axis.title=element_text(size=11)) + 
  theme(legend.position = "bottom") + theme(panel.background = element_rect(fill="grey90", color=NA))

ggsave(filename = here("noncognitive-cognitive/final/output/figureS1-child-nurture-models.png"), figureS1,
       width = 297,height = 210,units = "mm",  dpi = 320,device = "png")
ggsave(filename = here("noncognitive-cognitive/final/output/figureS1-child-nurture-models.pdf"), figureS1,
       width = 297,height = 210,units = "mm",  dpi = 320,device = "pdf")


### Figure S2: Model 3 for each test 
makespecificplot <- function(specoutcome) {
  myterms <- c("child_EA4","child_Cog","child_NonCog","child_WFEA",
               "mother_EA4","mother_Cog","mother_NonCog","mother_WFEA",
               "father_EA4","father_Cog","father_NonCog","father_WFEA",
               "sibparent_dev_EA4","sibparent_mean_EA4", "otherparent_EA4",
               "sibparent_dev_Cog","sibparent_mean_Cog", "otherparent_Cog",
               "sibparent_dev_NonCog","sibparent_mean_NonCog", "otherparent_NonCog",
               "sibparent_dev_WFEA","sibparent_mean_WFEA", "otherparent_WFEA",
               "sibparent_EA4", "sibparent_WFEA", "sibparent_Cog", "sibparent_NonCog",
               "othersibling_EA4","othersibling_Cog","othersibling_NonCog","othersibling_WFEA"
  )
  
  childterm <- "Child's\n PGI"
  motherterm <- "Mother's\n PGI"
  fatherterm <- "Father's\n PGI"
  devterm <- "Parent's\n Deviation\n from\n Parent-Sibling\n Pair Mean PGI"
  meanterm <- "Mean of\n Parent's\n and\n Sibling's\n PGIs"
  otherterm <- "Other\n Parent's\n PGI"
  parentterm <- "Parent's\n PGI"
  uncleterm <- "Parent's\n Sibling's\n PGI"
  names(myterms) <- c(rep(childterm,4),
                      rep(motherterm,4),
                      rep(fatherterm,4),
                      c(devterm, meanterm, otherterm),
                      c(devterm, meanterm, otherterm),
                      c(devterm, meanterm, otherterm),
                      c(devterm, meanterm, otherterm),
                      rep(parentterm,4),
                      rep(uncleterm,4)
  )
  termsdf <- data.frame(term = myterms, label = names(myterms)) %>% 
    mutate(label = factor(label, levels = c(childterm, motherterm, fatherterm, meanterm, devterm, parentterm, otherterm, uncleterm)))
  
  specplot <- results %>% 
    filter(outcome==specoutcome) %>%
    filter(design=="Sib") %>% 
    filter(model_def %in% c("child_model","nurture_model","sibling_model","uncle_model")) %>% 
    mutate(model_def = case_when(model_def=="child_model" ~ "Model 1",
                                 model_def=="nurture_model" ~ "Model 2",
                                 model_def=="sibling_model" ~ "Model 3",
                                 model_def=="uncle_model" ~ "Model 4")) %>% 
    filter(term %in% myterms) %>% 
    left_join(termsdf, "term") %>%
    mutate(PGI = case_when(
      PGIdef=="Cog/NonCog" ~ "Cognitive/Non-cognitive\n (Cog/NonCog) PGI",
      PGIdef=="EA4" ~ "Population Educational\n attainment 4 (EA4) PGI",
      PGIdef=="WFEA" ~"Within-family \nEducational attainment\n (WFEA) PGI"
    )) %>% 
    mutate(PGI = factor(PGI, levels = c("Population Educational\n attainment 4 (EA4) PGI",
                                        "Cognitive/Non-cognitive\n (Cog/NonCog) PGI",
                                        "Within-family \nEducational attainment\n (WFEA) PGI"))) %>% 
    mutate(shapefix = case_when(PGIdef=="EA4" ~ "SAME",
                                PGIdef=="WFEA" ~ "SAME",
                                str_detect(term,"NonCog") ~ "Non-cognitive",
                                str_detect(term,"Cog") ~ "Cognitive")) %>% 
    ggplot(mapping = aes(x=label,y=estimate,ymin=conf.low,ymax=conf.high, color=shapefix, shape=shapefix)) +
    geom_pointrange(position = position_dodge2(width=0.20)) +
    geom_abline(slope = 0, intercept=0) + 
    geom_abline(slope = 0, intercept=0.08, linetype=2) + 
    ylab(theylab) + xlab(thexlab) + 
    scale_color_manual(values=c("Non-cognitive"="orange","Cognitive"="blue","EA4"="black","SAME"="black")) +
    facet_grid(PGI~model_def, shrink = TRUE,scales="free_x",space = "free_x", margins = ) +
    theme_minimal() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    theme(strip.text.x = element_text(size = 11)) + 
    theme(strip.text.y = element_text(size = 11)) + 
    theme(axis.text=element_text(size=9)) +
    theme(axis.title=element_text(size=11)) + 
    theme(legend.position = "none") + theme(panel.background = element_rect(fill="grey90", color=NA))
  
  thefilename <- here(str_c("noncognitive-cognitive/final/output/figureS2_",specoutcome,".png"))
  ggsave(filename = thefilename, specplot,
         width = 297,height = 210,units = "mm",  dpi = 320,device = "png")
  thefilename <- here(str_c("noncognitive-cognitive/final/output/figureS2_",specoutcome,".pdf"))
  ggsave(filename = thefilename, specplot,
         width = 297,height = 210,units = "mm",  dpi = 320,device = "pdf")
  return(0)
}
map(c("Math5th","Math8th","Math9th",
      "Reading5th","Reading8th","Reading9th",
      "English5th","English8th"), makespecificplot)



###
### Main results in tabular form
###

# for EA4
mainmods_EA4 <- c("Model 1" = all_models[["Sib:EA4:Combined:child_model"]],
                  "Model 2" = all_models[["Sib:EA4:Combined:nurture_model"]],
                  "Model 3" = all_models[["Sib:EA4:Combined:sibling_model"]],
                  "Model 4" = all_models[["Sib:EA4:Combined:uncle_model"]])
keepcoefs <- c(
  "(Intercept)" = "(Intercept)",
  "Child PGI EA4" = "child_EA4",
  "Mother PGI EA4" = "mother_EA4",
  "Father PGI EA4" = "father_EA4",
  "Parent-Sibship Deviation EA4" = "sibparent_dev_EA4",
  "Parent-Sibship Mean EA4" = "sibparent_mean_EA4",
  "Other parent PGI EA4" = "otherparent_EA4",
  "Sibling parent PGI EA4" = "sibparent_EA4",
  "Uncle/aunt PGI EA4" = "othersibling_EA4",
  "SD(Child intercepts)" = "sd__(Intercept)",
  "SD(Tests)" = "sd__Observation"
)
maintab_EA4 <-  huxreg(mainmods_EA4,
                       coefs = keepcoefs, 
                       ci_level = 0.95, error_format = "[{conf.low}, {conf.high}]", error_pos = "right")
quick_html(maintab_EA4,file = here("noncognitive-cognitive/final/output/maintab_EA4.html"), open=FALSE)

# for WFEA
mainmods_WFEA <- c("Model 1" = all_models[["Sib:WFEA:Combined:child_model"]],
                   "Model 2" = all_models[["Sib:WFEA:Combined:nurture_model"]],
                   "Model 3" = all_models[["Sib:WFEA:Combined:sibling_model"]],
                   "Model 4" = all_models[["Sib:WFEA:Combined:uncle_model"]])
keepcoefs <- c(
  "(Intercept)" = "(Intercept)",
  "Child PGI WFEA" = "child_WFEA",
  "Mother PGI WFEA" = "mother_WFEA",
  "Father PGI WFEA" = "father_WFEA",
  "Parent-Sibship Deviation WFEA" = "sibparent_dev_WFEA",
  "Parent-Sibship Mean WFEA" = "sibparent_mean_WFEA",
  "Other parent PGI WFEA" = "otherparent_WFEA",
  "Sibling parent PGI WFEA" = "sibparent_WFEA",
  "Uncle/aunt PGI WFEA" = "othersibling_WFEA",
  "SD(Child intercepts)" = "sd__(Intercept)",
  "SD(Tests)" = "sd__Observation"
)
maintab_WFEA <-  huxreg(mainmods_WFEA,
                        coefs = keepcoefs, 
                        ci_level = 0.95, error_format = "[{conf.low}, {conf.high}]", error_pos = "right")
quick_html(maintab_WFEA,file = here("noncognitive-cognitive/final/output/maintab_WFEA.html"), open=FALSE)

# for cog/noncog
mainmods_Cog <- c("Model 1" = "Sib:Cog:Combined:child_model",
                  "Model 2" = "Sib:Cog:Combined:nurture_model",
                  "Model 3" = "Sib:Cog:Combined:sibling_model",
                  "Model 4" = "Sib:Cog:Combined:uncle_model")
keepcoefs <- c(
  "(Intercept)" = "(Intercept)",
  "Child PGI Cog" = "child_Cog",
  "Child PGI NonCog" = "child_NonCog",
  "Mother PGI Cog" = "mother_Cog",
  "Mother PGI NonCog" = "mother_NonCog",
  "Father PGI Cog" = "father_Cog",
  "Father PGI NonCog" = "father_NonCog",
  "Parent-Sibship Mean Cog" = "sibparent_mean_Cog",
  "Parent-Sibship Deviation Cog" = "sibparent_dev_Cog",
  "Parent-Sibship Mean NonCog" = "sibparent_mean_NonCog",
  "Parent-Sibship Deviation NonCog" = "sibparent_dev_NonCog",
  "Sibling parent PGI Cog" = "sibparent_Cog",
  "Sibling parent PGI NonCog" = "sibparent_NonCog",
  "Other parent PGI Cog" = "otherparent_Cog",
  "Other parent PGI NonCog" = "otherparent_NonCog",
  "Uncle/aunt PGI Cog" = "othersibling_Cog",
  "Uncle/aunt PGI NonCog" = "othersibling_NonCog",
  "SD(Child intercepts)" = "sd__(Intercept)",
  "SD(Tests)" = "sd__Observation"
)
maintab_Cog <-  huxreg(all_models[names(all_models) %in% mainmods_Cog],
                       coefs = keepcoefs, 
                       ci_level = 0.95, error_format = "[{conf.low}, {conf.high}]", error_pos = "right")
quick_html(maintab_Cog,file = here("noncognitive-cognitive/final/output/maintab_Cog.html"), open=FALSE)

###
### Hypothesis tests of nurture < dynastic
###

ntest_EA4 <- multcomp::glht(model = all_models[["Sib:EA4:Combined:sibling_model"]],
                      "sibparent_mean_EA4 - sibparent_dev_EA4 <= 0")

ntest_CogNonCog <-multcomp::glht(model = all_models[["Sib:Cog:Combined:sibling_model"]],
                      "sibparent_mean_NonCog+sibparent_mean_Cog - sibparent_dev_NonCog - sibparent_dev_Cog <= 0")

ntest_WFEA <- multcomp::glht(model = all_models[["Sib:WFEA:Combined:sibling_model"]],
                      "sibparent_mean_WFEA - sibparent_dev_WFEA <= 0")

htests <- list("EA4" = ntest_EA4 ,
                "CogNonCog" = ntest_CogNonCog,
               "WFEA" = ntest_WFEA) %>% 
  map(~tidy(.x,conf.int=T)) %>% 
  bind_rows(.id="Model") %>% 
  select(-null.value,-contrast)
quick_html(htests, file = here("noncognitive-cognitive/final/output/nurture-dynastic-tests.html"),open = FALSE)

###
### Additional analysis taking relatedness into account
###
qc <- fread("/tsd/p805/data/durable/data/genetics/MoBaPsychGen_v1/MoBaPsychGen_v1-ec-eur-batch-basic-qc-rel.kin")
names(qc)

parent_ids <- complete_trios %>% select(mother_SENTRIX_ID, father_SENTRIX_ID) %>% gather(role,IID) %>% filter(IID!="") %>% 
  select(-role) %>% distinct %>% pull(IID)

# All parent relatedness measures
allpars <- filter(qc, (ID1 %in% parent_ids) & (ID2 %in% parent_ids))

cousin_parents <-allpars %>% 
  filter(PropIBD>.125) %>% 
  filter(InfType!="FS") %>% 
  select(ID1,ID2) %>% 
  gather(role, IID) %>% 
  distinct %>% 
  pull(IID)

complete_trios_unrel <- complete_trios %>% 
  filter(!(mother_SENTRIX_ID %in% cousin_parents)) %>% 
  filter(!(father_SENTRIX_ID %in% cousin_parents))

sibtrios_unrel <- filter(sibtrios, !(mother_SENTRIX_ID %in% cousin_parents) & !(father_SENTRIX_ID %in% cousin_parents))

re_controls <- c(controls,"testname")
outcomes_nonames <- outcomes %>% setNames(nm=outcomes)
stdlong_unrel <- complete_trios_unrel %>% 
  select(any_of(outcomes_nonames), child_cohort, child_female, child_w19_0634_lnr,
         child_EA4, mother_EA4, father_EA4,
         child_WFEA, mother_WFEA, father_WFEA,
         child_Cog, mother_Cog, father_Cog,
         child_NonCog, mother_NonCog, father_NonCog,
         any_of(controls)) %>% 
  gather(testname, testscore, any_of(outcomes)) %>% 
  dplyr::select(-starts_with("NP")) %>% 
  filter(!is.na(testscore)) %>%
  mutate(testscore = as.numeric(testscore)) %>%
  filter(!is.na(child_EA4)) %>%
  filter(!is.na(mother_EA4)) %>%
  filter(!is.na(father_EA4))

stdavg_unrel <- stdlong_unrel %>% 
  group_by(child_w19_0634_lnr) %>% 
  mutate(avgtest = mean(testscore)) %>% slice(1) %>% 
  ungroup

siblong_unrel <- sibtrios_unrel %>% 
  gather(testname, testscore, any_of(outcomes)) %>% 
  dplyr::select(-starts_with("NP")) %>% 
  filter(!is.na(testscore)) %>%
  mutate(testscore = as.numeric(testscore))
sibavg <- siblong_unrel %>% 
  group_by(child_w19_0634_lnr) %>% 
  mutate(avgtest = mean(testscore)) %>% slice(1) %>% 
  ungroup

# EA4 Child model
unrel_child <- lmer(data=stdlong_unrel, formula = as.formula(gen_std_child_model("EA4")))

# EA4 nurture model
unrel_nurture <- lmer(data=stdlong_unrel, formula = as.formula(gen_std_nurture_model("EA4"))) 
                  
# EA4 uncle model
unrel_uncle <- lmer(data=siblong_unrel, formula=gen_sib_uncle_model("EA4"))

# EA4 sibling model
unrel_sibling <- lmer(data=siblong_unrel, formula=gen_sib_sibling_model("EA4"))

unrel_models <- list(
  "Model 1:Full" = sib_ea4_child_model, 
  "Model 2:Full" = sib_ea4_nurture_model, 
  "Model 3:Full" = sib_ea4_sibling_model, 
  "Model 4:Full" = sib_ea4_uncle_model, 
  "Model 1:Restricted" = unrel_child, 
  "Model 2:Restricted" = unrel_nurture, 
  "Model 3:Restricted" = unrel_sibling,
  "Model 4:Restricted" = unrel_uncle
)

keepcoefs <- c(
  " Child's \n PGI" = "child_EA4",
  " Mother's\nPGI " = "mother_EA4",
  "Father's\n PGI" = "father_EA4",
  " Parent's \n Deviation\n from\n Parent-Sibling\n Pair Mean PGI" = "sibparent_dev_EA4",
  " Mean of \n Parent's\n and\n Sibling's\n PGIs" = "sibparent_mean_EA4",
  "Other\n Parent's\n PGI" = "otherparent_EA4",
  " Parent's \n PGI" = "sibparent_EA4",
  " Parent's \n Sibling's\n PGI" = "othersibling_EA4"
)
unreldef <- map(unrel_models, ~tidy(.x,conf.int=TRUE)) %>% 
  bind_rows(.id="name") %>% 
  filter(term %in% keepcoefs) %>% 
  separate(col=name, into=c("Model","Sample"), sep=":") %>% 
  left_join(data.frame(label=names(keepcoefs),term=keepcoefs)) %>% 
  group_by(Sample, Model) %>% 
  arrange(term) %>% 
  ungroup

unrelplot <- unreldef %>% 
  ggplot(mapping = aes(x=label,y=estimate,ymin=conf.low,ymax=conf.high, shape=Sample, linetype=Sample)) +
  geom_pointrange(position = position_dodge2(width=0.20)) +
  geom_abline(slope = 0, intercept=0) + 
  geom_abline(slope = 0, intercept=0.08, linetype=2, color="black") + 
  ylab(theylab) + xlab(thexlab) + 
  facet_grid(~Model, shrink = TRUE,space="free_x",scales="free_x") +
  theme_minimal() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme(strip.text.x = element_text(size = 11)) + 
  theme(strip.text.y = element_text(size = 11)) + 
  theme(axis.text=element_text(size=9)) +
  theme(axis.title=element_text(size=11)) + 
  guides(linetype = guide_legend(title = "Sample"), shape = guide_legend(title="Sample")) + 
  theme(legend.position = "right") + theme(panel.background = element_rect(fill="grey90", color=NA))
ggsave(filename = here("noncognitive-cognitive/final/output/figureS11.pdf"), unrelplot,
       width = 297,height = 210,units = "mm",  dpi = 320,device = "pdf")
ggsave(filename = here("noncognitive-cognitive/final/output/figureS11.png"), unrelplot,
       width = 297,height = 210,units = "mm",  dpi = 320,device = "png")




###
### Checks for multicollinearity in models
###

library(performance)
multicol <- map(all_models, check_collinearity) %>% 
  bind_rows(.id="model") %>% 
  filter(Term %in% keepcoefs) %>% 
  select(model, Term, VIF) %>% rename(Model = model) %>% 
  arrange()
summary(multicol$VIF)

multicol %>% filter(str_detect(Model,"Sib:EA4:Combined:")) -> multicol_out
cat(knitr::kable(multicol_out, format="html"), file = here("noncognitive-cognitive/final/output/table_multicol.html"))




