traits = as.character(read.delim("/n/groups/price/kushal/singlecellLDSC/data/traits.txt", header=F)[,1])
traits2 = as.character(sapply(traits, function(x) return(strsplit(x, ".sumstats")[[1]][1])))
annot_cell = "/n/groups/price/kushal/singlecellLDSC/data/ANNOTATIONS/HVI_gene_score_Feb1_plus_Healthy_plus_All"
results_cell = "/n/groups/price/kushal/singlecellLDSC/data/LDSC_RESULTS/HVI_gene_score_Feb1_plus_Healthy_plus_All/baselineLD_MAF"
annot_modules = list.files(results_cell)

tt= array(0, c(length(annot_modules), length(traits2), 3, 4))
for(aa in 1:length(annot_modules)){
  for(kk in 1:length(traits2)){
    temp = read.delim(paste0(results_cell, "/", annot_modules[aa], "/",
                             traits2[kk], "_ldsc_postprocess.txt"), header=T)
    outt = as.matrix(temp[match(c("100kb", "5kb", "Roadmap_Enhancer"), rownames(temp)), c(1, 3, 4, 6)])
    tt[aa, kk, , ] = outt
  }
}
dimnames(tt)[[1]] = annot_modules
dimnames(tt)[[2]] = traits2
dimnames(tt)[[3]] = c("100kb", "5kb", "Roadmap_Enhancer")
dimnames(tt)[[4]] = c("taustar", "p-taustar", "E", "p.E")

save(tt, file = "/n/groups/price/kushal/singlecellLDSC/output/HVI_gene_score_Feb1_plus_Healthy_plus_All.rda")
