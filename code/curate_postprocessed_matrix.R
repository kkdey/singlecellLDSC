traits = as.character(read.delim("/n/groups/price/kushal/singlecellLDSC/data/traits.txt", header=F)[,1])
traits2 = unique(as.character(sapply(traits, function(x) return(strsplit(x, ".sumstats")[[1]][1]))))
annot_cell = "/n/groups/price/kushal/singlecellLDSC/data/ANNOTATIONS/Modules/healthy/celltype_enriched/blood"
results_cell = "/n/groups/price/kushal/singlecellLDSC/data/LDSC_RESULTS/Modules/healthy/celltype_enriched/blood/baselineLD_v2.1"
annot_modules = list.files(results_cell)

tt= array(0, c(length(annot_modules), length(traits2), 4, 6))
for(aa in 1:length(annot_modules)){
  for(kk in 1:length(traits2)){
    temp = read.delim(paste0(results_cell, "/", annot_modules[aa], "/",
                             traits2[kk], "_ldsc_postprocess.txt"), header=T)
    outt = as.matrix(temp[match(c("100kb", "5kb", "ABCT", "RoadmapT"), rownames(temp)), 1:6])
    tt[aa, kk, , ] = outt
  }
}
dimnames(tt)[[1]] = annot_modules
dimnames(tt)[[2]] = traits2
dimnames(tt)[[3]] = c("100kb", "5kb", "ABC", "Roadmap")
dimnames(tt)[[4]] = c("taustar", "se-taustar", "p-taustar", "E", "se.E", "p.E")

save(tt, file = "/n/groups/price/kushal/singlecellLDSC/output/healthy_celltype_enriched_blood_output_Mar4_2020_v2.1.rda")
