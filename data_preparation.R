library(here)
library(tidyverse)
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(broom)
library(data.table)
source(here("../my-useful-R-functions.R"))

###
### Create general trio data set
###
STEP1 = TRUE
STEP2 = TRUE
#rm(list=ls())
if (STEP1==TRUE) {
  ### Get MoBa IDs
  pregids <- fread(here("../../../data/moba/linkage/PDB2601_kobling_SSB_v12.csv"))
  parents <- fread(here("../../../data/moba/Original files/csv/PDB2601_SV_INFO_v12.csv"))
  
  ### Get basic demography data
  fasteoppl <- fread(here("../../../data/registers/original/csv/w19_0634_faste_oppl_ut.csv"))
  basic_info <- fasteoppl %>% 
    mutate(cohort = as.numeric(str_sub(foedsels_aar_mnd,1,4))) %>% 
    mutate(female = as.numeric(kjoenn)-1) %>% 
    select(w19_0634_lnr, mor_lnr, far_lnr, cohort, female, fodeland, doeds_aar_mnd)
  
  ### Get SENTRIXIDs
  child_sentrix <- fread(here("../../../data/moba/MoBaGenetics/key_CENTRIX_ID_Aug_22/PDB2601_MoBaGeneticsTot_Child_20220829.csv")) %>% 
    prefixit(.,prefix="child_",except=c("PREG_ID_2601","BARN_NR"))
  father_sentrix <- fread(here("../../../data/moba/MoBaGenetics/key_CENTRIX_ID_Aug_22/PDB2601_MoBaGeneticsTot_Father_20220829.csv"))     %>% 
    prefixit(.,prefix="father_",except=c("F_ID_2601"))
  mother_sentrix <- fread(here("../../../data/moba/MoBaGenetics/key_CENTRIX_ID_Aug_22/PDB2601_MoBaGeneticsTot_Mother_20220829.csv")) %>% 
    prefixit(.,prefix="mother_",except=c("M_ID_2601"))
  
  
  ### NAtl std. tests data in wide format
  nasjprov <- fread(here("../../../data/registers/original/csv/w19_0634_nasjonale_prover_ut.csv")) 
  nprover <-   select(nasjprov, w19_0634_lnr, PROVE, AARGANG, POENG, DELTATTSTATUS) %>% 
    filter(PROVE %in% c("NPENG08","NPREG08","NPLES08","NPENG05","NPREG09","NPLES05","NPREG05","NPLES09")) %>% 
    filter(DELTATTSTATUS=="D") %>% 
    group_by(w19_0634_lnr, PROVE) %>% 
    arrange(desc(AARGANG)) %>% 
    slice(1) %>% 
    group_by(PROVE) %>% 
    mutate(zresult = scale.default(POENG)) %>% 
    gather(info, value, -w19_0634_lnr, -PROVE) %>% 
    mutate(info = if_else(info=="zresult", "",info)) %>% 
    unite(col = varname, PROVE, info, sep="", remove=T) %>% 
    spread(varname, value)
  
  ### Link ID numbers together
  kids <- filter(pregids, rolle=="SU2PT_CHILD") %>% 
    select(-rolle) %>% rename(child_w19_0634_lnr = w19_0634_lnr)
  pars <- filter(pregids, rolle!="SU2PT_CHILD") %>% 
    select(-barn_nr) %>% 
    spread(rolle, w19_0634_lnr) %>% 
    rename(father_w19_0634_lnr = SU2PT_FATHER, mother_w19_0634_lnr = SU2PT_MOTHER)
  all_families <- kids %>%
    left_join(pars) %>% 
    left_join(parents, by=c("PREG_ID_2601")) %>% 
    rename(moba_child_cohort = FAAR) %>% 
    left_join(child_sentrix, by=c("PREG_ID_2601",barn_nr = "BARN_NR")) %>% 
    left_join(mother_sentrix, by=c("M_ID_2601")) %>% 
    left_join(father_sentrix, by=c("F_ID_2601"))
  
  # Read in scores
  ea4_scores <- fread(here("noncognitive-cognitive/revision/scores/EA4_additive_excl_23andMeauto_scores.tsv")) %>% 
    rename(EA4 = final_pred_auto) %>% select(-FID) %>% 
    mutate_at(vars(contains("EA4")), ~scale.default(.x)) %>% 
    mutate(EA4 = EA4*-1)
  names(ea4_scores)
  cog_scores <- fread(here("noncognitive-cognitive/revision/scores/Cog_GWAS_excl23andMe.txt_auto_scores.tsv")) %>% 
    rename(Cog = final_pred_auto) %>% select(-FID) %>% 
    mutate_at(vars(contains("Cog")), ~scale.default(.x)) %>% 
    mutate(Cog = Cog*-1)
  names(cog_scores)
  noncog_scores <- fread(here("noncognitive-cognitive/revision/scores/NonCog_GWAS_excl23andMe.txt_auto_scores.tsv")) %>% 
    rename(NonCog = final_pred_auto) %>% select(-FID) %>% 
    mutate_at(vars(contains("NonCog")), ~scale.default(.x)) %>% 
    mutate(NonCog = NonCog*-1)
  names(noncog_scores)
  wfea_scores <- fread(here("noncognitive-cognitive/revision/scores/WFEA.txt_auto_scores.tsv")) %>% 
    rename(WFEA = final_pred_auto) %>% select(-FID) %>% 
    mutate_at(vars(contains("WFEA")), ~scale.default(.x))
  names(wfea_scores)
  
  pcabatch <- fread(here("../../../data/genetics/MoBaPsychGen_v1/MoBaPsychGen_v1-ec-eur-batch-basic-qc-cov-noMoBaIDs.txt")) %>% 
    select(SENTRIXID, IID, FID, BATCH = genotyping_batch, num_range("PC",1:10))
  # Add scores to the mix
  all_families_with_pgis <- all_families %>% 
    left_join(prefixit(pcabatch, prefix="child_"), by=c(child_SENTRIX_ID = "child_SENTRIXID")) %>% 
    left_join(prefixit(pcabatch, prefix="mother_"), by=c(mother_SENTRIX_ID = "mother_SENTRIXID")) %>% 
    left_join(prefixit(pcabatch, prefix="father_"), by=c(father_SENTRIX_ID = "father_SENTRIXID")) %>% 
    left_join(prefixit(ea4_scores, prefix="child_"), by=c(child_SENTRIX_ID = "child_IID")) %>% 
    left_join(prefixit(ea4_scores, prefix="mother_"), by=c(mother_SENTRIX_ID = "mother_IID")) %>% 
    left_join(prefixit(ea4_scores, prefix="father_"), by=c(father_SENTRIX_ID = "father_IID")) %>% 
    left_join(prefixit(cog_scores, prefix="child_"), by=c(child_SENTRIX_ID = "child_IID")) %>% 
    left_join(prefixit(cog_scores, prefix="mother_"), by=c(mother_SENTRIX_ID = "mother_IID")) %>% 
    left_join(prefixit(cog_scores, prefix="father_"), by=c(father_SENTRIX_ID = "father_IID")) %>% 
    left_join(prefixit(noncog_scores, prefix="child_"), by=c(child_SENTRIX_ID = "child_IID")) %>% 
    left_join(prefixit(noncog_scores, prefix="mother_"), by=c(mother_SENTRIX_ID = "mother_IID")) %>% 
    left_join(prefixit(noncog_scores, prefix="father_"), by=c(father_SENTRIX_ID = "father_IID")) %>% 
    left_join(prefixit(wfea_scores, prefix="child_"), by=c(child_SENTRIX_ID = "child_IID")) %>% 
    left_join(prefixit(wfea_scores, prefix="mother_"), by=c(mother_SENTRIX_ID = "mother_IID")) %>% 
    left_join(prefixit(wfea_scores, prefix="father_"), by=c(father_SENTRIX_ID = "father_IID"))
  nrow(all_families_with_pgis)
  
  ### Link together various data components
  all_families_with_tests <- all_families_with_pgis %>% 
    left_join(nprover, by=c(child_w19_0634_lnr = "w19_0634_lnr")) %>% 
    left_join(prefixit(basic_info, prefix="child_"), by="child_w19_0634_lnr") %>% 
    left_join(prefixit(basic_info, prefix="mother_"), by="mother_w19_0634_lnr") %>% 
    left_join(prefixit(basic_info, prefix="father_"), by="father_w19_0634_lnr")
  nrow(all_families_with_tests)
  
  ### Complete sample
  complete_trios <- all_families_with_tests %>% 
    group_by(child_w19_0634_lnr) %>% 
    arrange(desc(child_EA4)) %>% 
    arrange(desc(mother_EA4)) %>% 
    arrange(desc(father_EA4)) %>% 
    slice(1) %>% 
    ungroup %>% 
    select(everything())
  nrow(complete_trios)
  
  fwrite(complete_trios, here("noncognitive-cognitive/2ndrev/data/complete-trios.csv"))
  message("Step 1. General trio data set created")
} else {
  complete_trios <- fread(here("noncognitive-cognitive/2ndrev/data/complete-trios.csv"))
  message("Step 1. Read previously generated trio data set.")
}

###
### Parent-sib data set
###
### Create sib-means and deviations in PGIs for MoBa parents

if (STEP2==TRUE) {
  outcomes <- c("NPREG05","NPREG08","NPREG09","NPLES05","NPLES08","NPLES09","NPENG05","NPENG08")
  names(outcomes) <- c("Math5th","Math8th","Math9th","Reading5th","Reading8th","Reading9th","English5th","English8th")
  
  allparents <- select(complete_trios, contains("w19_0634_lnr"), 
                       contains("mor_lnr"), 
                       contains("far_lnr"), 
                       ends_with("_EA4"), ends_with("_Cog"), 
                       ends_with("_NonCog"), ends_with("_WFEA")
                  ) %>% 
    gather(varname, value, starts_with("father_"), starts_with("mother_")) %>% 
    mutate(parent = str_sub(varname, start=1, end=6)) %>% 
    mutate(varname = str_sub(varname, start=8)) %>% 
    pivot_wider(names_from=varname, values_from=value) %>% 
    select(-starts_with("child")) %>% 
    filter(!is.na(w19_0634_lnr)) %>% 
    distinct
  
  
  sibmeans <- allparents %>% 
    arrange(mor_lnr, far_lnr) %>% 
    group_by(mor_lnr, far_lnr) %>% 
    filter(mor_lnr!="",far_lnr!="") %>% 
    mutate_at(c("EA4","Cog","NonCog","WFEA"),as.numeric) %>% 
    filter(!is.na(EA4)) %>% 
    mutate(sibgroup_n = n()) %>% 
    filter(sibgroup_n>1) %>% 
    mutate(rn = row_number()) %>% 
    mutate(othersibling_w19_0634_lnr = if_else(rn==1,lead(w19_0634_lnr),lag(w19_0634_lnr))) %>% 
    mutate(othersibling_EA4 = if_else(rn==1,lead(EA4),lag(EA4))) %>% 
    mutate(othersibling_WFEA = if_else(rn==1,lead(WFEA),lag(WFEA))) %>% 
    mutate(othersibling_Cog = if_else(rn==1,lead(Cog),lag(Cog))) %>% 
    mutate(othersibling_NonCog = if_else(rn==1,lead(NonCog),lag(NonCog))) %>% 
    select(-rn) %>% 
    mutate(sibparent_mean_EA4 = mean(EA4)) %>% 
    mutate(sibparent_mean_Cog = mean(Cog)) %>% 
    mutate(sibparent_mean_NonCog = mean(NonCog)) %>% 
    mutate(sibparent_mean_WFEA = mean(WFEA)) %>% 
    ungroup %>% 
    mutate(sibparent_dev_EA4 = EA4-sibparent_mean_EA4) %>% 
    mutate(sibparent_dev_Cog = Cog-sibparent_mean_Cog) %>% 
    mutate(sibparent_dev_NonCog = NonCog-sibparent_mean_NonCog) %>% 
    mutate(sibparent_dev_WFEA = WFEA-sibparent_mean_WFEA) %>% 
    rename(sibparent_w19_0634_lnr = w19_0634_lnr) %>% 
    rename(sibparent_sex = parent) %>% 
    rename(sibparent_EA4 = EA4) %>% 
    rename(sibparent_Cog = Cog) %>% 
    rename(sibparent_NonCog = NonCog) %>% 
    rename(sibparent_WFEA = WFEA) %>% 
    mutate(sibparent = 1) %>% 
    select(-mor_lnr, -far_lnr)
  head(sibmeans)
  nrow(sibmeans)
  
  sibtrios <- complete_trios %>% 
    select(contains("w19_0634_lnr"), one_of(outcomes), 
           ends_with("_EA4"),ends_with("_Cog"),
           ends_with("_NonCog"), ends_with("WFEA"), ends_with("SENTRIX_ID"),
           contains("BATCH"), contains("PC"), child_female, ends_with("cohort")) %>% 
    filter(!is.na(child_EA4),!is.na(father_EA4),!is.na(mother_EA4)) %>% 
    left_join(prefixit(sibmeans, prefix="father_"), by=c(father_w19_0634_lnr = "father_sibparent_w19_0634_lnr")) %>%
    left_join(prefixit(sibmeans, prefix="mother_"), by=c(mother_w19_0634_lnr = "mother_sibparent_w19_0634_lnr")) %>% 
    replace_na(replace = list(father_sibparent = 0,
                              mother_sibparent=0)) %>% 
    mutate_at(vars(contains("EA4"),contains("NonCog"),contains("Cog"), contains("WFEA")), as.numeric) %>% 
    mutate(numsibparents = as.numeric(mother_sibparent) + as.numeric(father_sibparent)) %>%
    filter(numsibparents==1) %>%
    mutate(othersibling_EA4 = if_else(mother_sibparent==1,
                                   mother_othersibling_EA4, father_othersibling_EA4)) %>%
    mutate(othersibling_WFEA = if_else(mother_sibparent==1,
                                    mother_othersibling_WFEA, father_othersibling_WFEA)) %>%
    mutate(othersibling_Cog = if_else(mother_sibparent==1,
                                   mother_othersibling_Cog, father_othersibling_Cog)) %>%
    mutate(othersibling_NonCog = if_else(mother_sibparent==1,
                                      mother_othersibling_NonCog, father_othersibling_NonCog)) %>%
    mutate(sibparent_EA4 = if_else(mother_sibparent==1,
                                   mother_EA4, father_EA4)) %>%
    mutate(sibparent_WFEA = if_else(mother_sibparent==1,
                                   mother_WFEA, father_WFEA)) %>%
    mutate(sibparent_Cog = if_else(mother_sibparent==1,
                                   mother_Cog, father_Cog)) %>%
    mutate(sibparent_NonCog = if_else(mother_sibparent==1,
                                   mother_NonCog, father_NonCog)) %>%
    mutate(sibparent_mean_EA4 = if_else(mother_sibparent==1,
                                        mother_sibparent_mean_EA4, father_sibparent_mean_EA4)) %>%
    mutate(sibparent_dev_EA4 = if_else(mother_sibparent==1,
                                       mother_sibparent_dev_EA4, father_sibparent_dev_EA4)) %>%
    mutate(otherparent_EA4 = if_else(mother_sibparent==1,
                                     father_EA4, mother_EA4)) %>% 
    mutate(sibparent_mean_Cog = if_else(mother_sibparent==1,
                                        mother_sibparent_mean_Cog, father_sibparent_mean_Cog)) %>%
    mutate(sibparent_dev_Cog = if_else(mother_sibparent==1,
                                       mother_sibparent_dev_Cog, father_sibparent_dev_Cog)) %>%
    mutate(otherparent_Cog = if_else(mother_sibparent==1,
                                     father_Cog, mother_Cog)) %>% 
    mutate(sibparent_mean_NonCog = if_else(mother_sibparent==1,
                                        mother_sibparent_mean_NonCog, father_sibparent_mean_NonCog)) %>%
    mutate(sibparent_dev_NonCog = if_else(mother_sibparent==1,
                                       mother_sibparent_dev_NonCog, father_sibparent_dev_NonCog)) %>%
    mutate(otherparent_NonCog = if_else(mother_sibparent==1,
                                     father_NonCog, mother_NonCog)) %>% 
    mutate(sibparent_mean_WFEA = if_else(mother_sibparent==1,
                                           mother_sibparent_mean_WFEA, father_sibparent_mean_WFEA)) %>%
    mutate(sibparent_dev_WFEA = if_else(mother_sibparent==1,
                                          mother_sibparent_dev_WFEA, father_sibparent_dev_WFEA)) %>%
    mutate(otherparent_WFEA = if_else(mother_sibparent==1,
                                        father_WFEA, mother_WFEA)) %>% 
    mutate(partnerrank_EA4 = percent_rank(otherparent_EA4))
  nrow(sibtrios)  
  
  fwrite(sibtrios, here("noncognitive-cognitive/2ndrev/data/sibtrios.csv"))
  message("Step 2. Parent-sib data set created")
} else {
  sibtrios <- fread(here("noncognitive-cognitive/2ndrev/data/sibtrios.csv"))
  message("Step 2. Read previously generated parent-sib dataset.")
}
