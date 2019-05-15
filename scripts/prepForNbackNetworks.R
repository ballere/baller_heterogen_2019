library(visreg)
library(mgcv)
library(tableone)
library(dplyr)
library(plm)
library(MatchIt)
library(tidyr)
library(ggplot2)
library(reshape)
library(emmeans)

##### This script grabs the subset of people who have imaging data from the original Hydra Sample  #####

#######################################################
############ READ IN, MERGE AND SUBSET DATA############
#######################################################

#read in csvs
demographics <- read.csv("/Users/eballer/BBL/from_chead/ballerDepHeterogen/data/n9498_demographics_go1_20161212.csv", header = TRUE, sep = ",") #from /data/joy/BBL/projevolumes/ballerDepHeterogen/data/n9498_demographics_go1_20161212.csv
cnb_scores <- read.csv("/Users/eballer/BBL/from_chead/ballerDepHeterogen/data/n9498_cnb_zscores_fr_20170202.csv", header = TRUE, sep = ",") #from /data/joy/BBL/projects/ballerDepHeterogen/data/n9498_cnb_zscores_fr_20170202.csv
health <- read.csv("/Users/eballer/BBL/from_chead/ballerDepHeterogen/data/n9498_health_20170405.csv", header = TRUE, sep = ",") #from /data/joy/BBL/projects/ballerDepHeterogen/data/n9498_health_20170405.csv
psych_summary <- read.csv("/Users/eballer/BBL/from_chead/ballerDepHeterogen/data/n9498_goassess_psych_summary_vars_20131014.csv", header = TRUE, sep = ",") #from /data/joy/BBL/projects/ballerDepHeterogen/data/n9498_goassess_psych_summary_vars_20131014.csv
t1_qa_scores <- read.csv("/Users/eballer/BBL/from_chead/ballerDepHeterogen/data/neuroimaging/t1struct/n1601_t1QaData_20170306.csv", header = TRUE, sep = ",") 
fc_nback_qa_scores <- read.csv("/Users/eballer/BBL/from_chead/ballerDepHeterogen/data/neuroimaging/nback/n1601_NbackQAData_20170427.csv", header = TRUE, sep = ",")

#read in Hydra cluster data -> AG (all_gender), M(males), F(females) 
hydra_AG_matched_clusters <- read.csv("/Users/eballer/BBL/from_chead/ballerDepHeterogen/hydra_output_from_cbica/hydra_matched/20181121_OLD_HYDRA_with_new_script_from_toni/CogData_age_and_sex_n1424_matched_imaging_plus_non_imaging_replication_20181121_HydraSubtypes.csv", header = TRUE, sep = ",")

#Coerce hydra values into factor format
hydra_names <- names(hydra_AG_matched_clusters[2:11])
hydra_AG_matched_clusters[hydra_names] <- lapply(hydra_AG_matched_clusters[hydra_names], factor)


############
#prepping and merging 
############

#remove people with NA for race, age, or sex.  START WITH N = 9498
demographics_noNA_race <- demographics[!is.na(demographics$race),] #everyone has a race, N = 9498
demographics_noNA_race_age <- demographics_noNA_race[!is.na(demographics_noNA_race$ageAtClinicalAssess1),] # 86 people do not have age at clinical assessment.  N = 9412
demographics_noNA_race_age_sex <- demographics_noNA_race_age[!is.na(demographics_noNA_race_age$sex),] #everyone has a sex, N = 9412
demographics_noNA_race_age_andCNBage_sex <- demographics_noNA_race_age_sex[!is.na(demographics_noNA_race_age_sex$ageAtCnb1),] #6 people do not have ageAtCnb1, N = 9406

#remove people with NA for depression or total psych score, START WITH N = 9498
psych_summary_no_NA_dep <- psych_summary[!is.na(psych_summary$smry_dep),] #take out those with NA for depression, 87 people N = 9411
psych_summary_no_NA_dep_and_smry_psych_overall <- psych_summary_no_NA_dep[!is.na(psych_summary_no_NA_dep$smry_psych_overall_rtg),] #take out those with NA for overall psych rtg, no additional people lost, N = 9411

#remove people with excluded imaging

t1_struct_include = subset.data.frame(t1_qa_scores, (t1Exclude == 0)) # lost 51, N = 1540
fc_include <- subset(fc_nback_qa_scores, nbackExclude ==0) #Lost 478, N 1123

#remove scanID from fc_include so the merge doesn't mess up
fc_include$scanid = NULL
t1_fc_include <- merge(t1_struct_include, fc_include, by = "bblid") #lost 12, new N = 1111 

#merge demographics and imaging #this is if we want to include people without full demographic data
dem_cnb <- merge(demographics_noNA_race_age_andCNBage_sex, cnb_scores, by = "bblid") #merge demographics and cnb, N = 9406
psych_health <- merge(psych_summary_no_NA_dep_and_smry_psych_overall, health, by = "bblid") #merge psych and health, N = 9411
dem_cnb_psych_health <- merge(dem_cnb, psych_health, by = "bblid") #merge all 4 csvs, lost 1 person (had demographics, but no psych ratings): N = 9405
t1_fc_dem_cnb_psych_health_merged <- merge(t1_fc_include, dem_cnb_psych_health, by = "bblid") #N = 1111 when all combined

#t1_qa_and_jlfvolume = merge(t1_struct_include, jlfvolume_summary, by = "bblid")
#dem_psych_health_qa_volume_merged <- merge (dem_cnb_psych_health, t1_qa_and_jlfvolume, by = "bblid") 
#dem_cnb_psych_health_factor_scores_merged <- merge(dem_cnb_psych_health_merged, cnb_factor_scores, by = "bblid")

#make subsets
subset_just_dep_and_no_medicalratingExclude <- subset.data.frame(t1_fc_dem_cnb_psych_health_merged, (medicalratingExclude == 0) & (smry_dep == 4), select = "bblid") 
subset_no_psych_no_medicalratingExclude <- subset.data.frame(t1_fc_dem_cnb_psych_health_merged, (medicalratingExclude == 0) & (smry_psych_overall_rtg < 4), select = "bblid") 
subset_dep_or_no_psych_and_no_medicalratingExclude <- subset.data.frame(t1_fc_dem_cnb_psych_health_merged, (medicalratingExclude == 0) & ((smry_dep == 4) | (smry_psych_overall_rtg <4))) 

#would binarize depression smry score to -1 (less than 4, not depressed) and 1 (score 4 , depressed)
dep_binarized <- ifelse(subset_dep_or_no_psych_and_no_medicalratingExclude$smry_dep == 4, 1, -1)
subset_dep_or_no_psych_and_no_medicalratingExclude_DEPBINARIZED <- cbind(subset_dep_or_no_psych_and_no_medicalratingExclude, dep_binarized) #N = 3284

#make depression and gender into factor scores
subset_dep_or_no_psych_and_no_medicalratingExclude_DEPBINARIZED$dep_binarized <- as.factor(subset_dep_or_no_psych_and_no_medicalratingExclude_DEPBINARIZED$dep_binarized)
subset_dep_or_no_psych_and_no_medicalratingExclude_DEPBINARIZED$sex <- as.factor(subset_dep_or_no_psych_and_no_medicalratingExclude_DEPBINARIZED$sex)

#divide ageAtCNB by 12 for age
subset_dep_or_no_psych_and_no_medicalratingExclude_DEPBINARIZED$age_in_years <- subset_dep_or_no_psych_and_no_medicalratingExclude_DEPBINARIZED$ageAtCnb1/12

#race binarized for plotting purposes, caucasian 1, non-caucasion 0
subset_dep_or_no_psych_and_no_medicalratingExclude_DEPBINARIZED$race_binarized <- ifelse(as.factor(subset_dep_or_no_psych_and_no_medicalratingExclude_DEPBINARIZED$race) == 1, 1, 0)

#add age squared mean-centered term for use in LM
subset_dep_or_no_psych_and_no_medicalratingExclude_DEPBINARIZED$ageSq <- as.numeric(I(scale(subset_dep_or_no_psych_and_no_medicalratingExclude_DEPBINARIZED$age_in_years, scale = FALSE, center = TRUE)^2))

#remove people with NA in their cognitive measures
subset_dep_or_no_psych_and_no_medicalratingExclude_DEPBINARIZED <- subset_dep_or_no_psych_and_no_medicalratingExclude_DEPBINARIZED[complete.cases(subset_dep_or_no_psych_and_no_medicalratingExclude_DEPBINARIZED[14:39]),]

################################################
### Merge clustering data with subsetted data ##
################################################

#AG (all gender), M (males), F (females) 
#AG Matched - n 1424, after merge- n = 325

subset_with_clusters_AG_matched <- merge(subset_dep_or_no_psych_and_no_medicalratingExclude_DEPBINARIZED, hydra_AG_matched_clusters, by = "bblid")

#get scanID lists for each cluster
scanid_TD <- data.frame(subset_with_clusters_AG_matched$scanid[which(subset_with_clusters_AG_matched$Hydra_k3 == -1)])
scanid_cluster1 <- data.frame(subset_with_clusters_AG_matched$scanid[which(subset_with_clusters_AG_matched$Hydra_k3 == 1)])
scanid_cluster2 <- data.frame(subset_with_clusters_AG_matched$scanid[which(subset_with_clusters_AG_matched$Hydra_k3 == 2)])
scanid_cluster3 <- data.frame(subset_with_clusters_AG_matched$scanid[which(subset_with_clusters_AG_matched$Hydra_k3 == 3)])

#add column indicating that all data should be included for voxelwrapper
subset_with_clusters_AG_matched$include <- c(rep(1, nrow(subset_with_clusters_AG_matched)))

#save to rds 
saveRDS(object = subset_with_clusters_AG_matched, file = "/Users/eballer/BBL/from_chead/ballerDepHeterogen/ballerDepHeterogenScripts/PrepForFunctionalNetworks/subset_with_T1_FC_and_dem_with_clusters_20181121_OLD_HYDRA_with_new_script_from_toni.rds")

