library(data.table)
library(rmeta)
annot_cell = "/n/groups/price/kushal/DiseaseNet/data/ANNOTATIONS/S2G/DeepSEA_Blood_S2G_only_Jan22"
results_cell = "/n/groups/price/kushal/DiseaseNet/data/LDSC_RESULTS/S2G/DeepSEA_Blood_S2G_only_Jan22/baselineLD_v2.1_no_Coding_no_Promoter_no_conserved"
annot_names = list.files(results_cell)
annot_idx = 1



all_traits = c('UKB_460K.body_BMIz','UKB_460K.cov_EDU_YEARS','UKB_460K.lung_FVCzSMOKE','UKB_460K.cov_SMOKING_STATUS',
               'UKB_460K.mental_NEUROTICISM','UKB_460K.blood_WHITE_COUNT','PASS_Years_of_Education2','UKB_460K.bp_SYSTOLICadjMEDz',
               'UKB_460K.body_HEIGHTz','UKB_460K.other_MORNINGPERSON','UKB_460K.body_WHRadjBMIz','UKB_460K.lung_FEV1FVCzSMOKE',
               'UKB_460K.repro_MENARCHE_AGE','UKB_460K.blood_RED_COUNT','UKB_460K.blood_PLATELET_COUNT','UKB_460K.bmd_HEEL_TSCOREz',
               'UKB_460K.blood_EOSINOPHIL_COUNT','PASS_Schizophrenia','UKB_460K.blood_RBC_DISTRIB_WIDTH','PASS_Height1','PASS_BMI1',
               'UKB_460K.disease_T2D','PASS_AgeFirstBirth','UKB_460K.disease_RESPIRATORY_ENT','UKB_460K.body_BALDING1','UKB_460K.disease_HYPOTHYROIDISM_SELF_REP',
               'UKB_460K.disease_ALLERGY_ECZEMA_DIAGNOSED','UKB_460K.disease_HI_CHOL_SELF_REP','UKB_460K.repro_MENOPAUSE_AGE','PASS_HDL','UKB_460K.pigment_SUNBURN',
               'PASS_NumberChildrenEverBorn','PASS_Anorexia','PASS_LDL','PASS_Crohns_Disease','PASS_DS','PASS_Ever_Smoked','UKB_460K.pigment_HAIR',
               'PASS_Rheumatoid_Arthritis','PASS_Type_2_Diabetes','PASS_Autism','UKB_460K.pigment_TANNING','PASS_Ulcerative_Colitis',
               'UKB_460K.disease_DERMATOLOGY','PASS_Coronary_Artery_Disease','UKB_460K.disease_AID_SURE','UKB_460K.pigment_SKIN')

blood_traits = c("UKB_460K.blood_RBC_DISTRIB_WIDTH", "UKB_460K.blood_RED_COUNT", "UKB_460K.blood_WHITE_COUNT",
                 "UKB_460K.blood_PLATELET_COUNT", "UKB_460K.blood_EOSINOPHIL_COUNT")

autoimmune_traits = c("UKB_460K.disease_AID_SURE", "PASS_Ulcerative_Colitis", "PASS_Crohns_Disease",
                      "PASS_Rheumatoid_Arthritis", "PASS_Celiac",  "PASS_Lupus", "PASS_Type_1_Diabetes",
                      "PASS_IBD", "PASS_Primary_biliary_cirrhosis")


brain_traits = c("PASS_Ever_Smoked", "UKB_460K.cov_SMOKING_STATUS", "UKB_460K.mental_NEUROTICISM", "UKB_460K.repro_MENARCHE_AGE",
                 "PASS_Years_of_Education2", "PASS_DS", "PASS_Schizophrenia", "UKB_460K.body_WHRadjBMIz",
                 "PASS_BMI1", "UKB_460K.body_BMIz")

immune_traits = c("PASS_Celiac", "PASS_Crohns_Disease", "PASS_IBD", "PASS_Lupus",
                  "PASS_Primary_biliary_cirrhosis", "PASS_Rheumatoid_Arthritis",
                  "PASS_Type_1_Diabetes", "PASS_Ulcerative_Colitis",
                  "UKB_460K.disease_ASTHMA_DIAGNOSED", "UKB_460K.disease_ALLERGY_ECZEMA_DIAGNOSED",
                  "PASS_Multiple_sclerosis", "UKB_460K.body_BMIz", "UKB_460K.disease_HYPERTENSION_DIAGNOSED",
                  "PASS_Triglycerides", "PASS_LDL", "PASS_HDL",
                  "UKB_460K.bp_DIASTOLICadjMEDz", "UKB_460K.bp_SYSTOLICadjMEDz",
                  "PASS_Alzheimer", "PASS_Anorexia", "PASS_Bipolar_Disorder",
                  "PASS_Schizophrenia", "UKB_460K.mental_NEUROTICISM", "UKB_460K.disease_AID_SURE",
                  "PASS_Type_2_Diabetes")


get_sd_annot = function(cell_path, annot_index = 1, flag=0){
  if(flag == 0){
    if(file.exists(paste0(cell_path, "/", "sd_annot_", annot_index, ".rda"))){
      sd_annot = get(load(paste0(cell_path, "/", "sd_annot_", annot_index, ".rda")))
      return(sd_annot)
    }else{
      flag = 1
    }}

  if(flag == 1){
    num = 0
    den = 0
    ll <- list.files(cell_path, pattern = ".annot.gz")
    for(m in 1:length(ll)){
      dat <- data.frame(fread(paste0("zcat ", cell_path, "/", ll[m])))
      num = num  + (nrow(dat)-1) * var(dat[,4+annot_index])
      den = den + (nrow(dat)-1)
      rm(dat)
    }
  }

  estd_sd_annot = sqrt(num/den)
  save(estd_sd_annot, file = paste0(cell_path, "/", "sd_annot_", annot_index, ".rda"))
  return(estd_sd_annot)
}


run_single_std_enrichment_analysis = function(annot_cell,
                                              results_cell,
                                              annotation,
                                              traits,
                                              index_in_results=1,
                                              flag=1){
  enrich_table = matrix(0, length(index_in_results), 3)
  cell_path = paste0(annot_cell, "/", annotation)
  sd_annot1=get_sd_annot(cell_path, annot_index=index_in_results, flag = flag)
  res = paste0(results_cell, "/", annotation, "/", traits[1], ".sumstats.results")
  tab2 = read.table(res,header=T)
  cat("Number of annotations together with baseline : ", nrow(tab2), "\n")
  annot_names = as.character(tab2$Category[index_in_results])
  Mref = 5961159
  for(id in 1:length(index_in_results)){
    meta_enr        = NULL;
    meta_enrstat    = NULL;
    for(trait_id in 1:length(traits)){
      result.file=paste0(results_cell, "/", annotation, "/", traits[trait_id], ".sumstats.results")
      res=read.table(result.file,header=T)
      logfile = paste(results_cell, "/", annotation, "/", traits[trait_id],".sumstats.log", sep="")
      log = read.table(logfile,h=F,fill=T)
      h2g = as.numeric(as.character(log[which(log$V4=="h2:"),5]))
      myenrstat    = (sd_annot1*h2g/Mref)*((res[index_in_results[id],3]/res[index_in_results[id],2])-(1-res[index_in_results[id],3])/(1-res[index_in_results[id],2])) #step1
      myenrstat_z  = qnorm(res[index_in_results[id],7]/2) #step2
      myenrstat_sd = myenrstat/myenrstat_z #step3
      meta_enrstat = rbind(meta_enrstat   , c(myenrstat, myenrstat_sd));
      meta_enr     = rbind(meta_enr, c(sd_annot1*res[index_in_results[id],5],
                                       sd_annot1*res[index_in_results[id],6] ));
    }
    test_eni1=meta.summaries(meta_enr[,1], meta_enr[,2],method="random")
    test_eni2=meta.summaries(meta_enrstat[,1], meta_enrstat[,2],method="random")
    cat("Printing enrichment results for annotation:", annot_names[id], "\n")
    cat(test_eni1$summary, " ",test_eni1$se.summary, " ", 2*pnorm(-abs(test_eni2$summary/test_eni2$se.summary)), "\n");
    enrich_table[id, ] = c(test_eni1$summary, test_eni1$se.summary, 2*pnorm(-abs(test_eni2$summary/test_eni2$se.summary)))
  }
  rownames(enrich_table) = annot_names
  return(enrich_table)
}

out3 = c()
for(m in 1:length(annot_names)){
  out3 = rbind(out3, run_single_std_enrichment_analysis(annot_cell, results_cell, annotation = annot_names[m],
                                                        traits = c(blood_traits, autoimmune_traits),
                                                        index_in_results = 1, flag=1))
}

rownames(out3) = annot_names

