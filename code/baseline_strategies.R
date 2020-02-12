

#######################  Baseline strategies   #######################################

for(numchr in 1:22){
  base1 = data.frame(fread(paste0("zcat /n/groups/price/kushal/singlecellLDSC/data/ANNOTATIONS/Baselines/",
                                  "baselineLD_v2.1", "/",
                                  "baselineLD", ".", numchr, ".annot.gz")))
  base2 = base1[,c(1:4, 5, 60:71, 72)]
  if(!dir.exists(paste0("/n/groups/price/kushal/singlecellLDSC/data/ANNOTATIONS/Baselines/", "baselineLD_MAF"))){
    dir.create(paste0("/n/groups/price/kushal/singlecellLDSC/data/ANNOTATIONS/Baselines/", "baselineLD_MAF"))
  }
  write.table(base2,
              file = gzfile(paste0("/n/groups/price/kushal/singlecellLDSC/data/ANNOTATIONS/Baselines/",
                                   "baselineLD_MAF", "/",
                                   "baselineLD", ".",
                                   numchr, ".annot.gz")),
              quote=FALSE, row.names=FALSE)
  cat("We are at chrom:", numchr, "\n")
}

