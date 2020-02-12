

########################  Log fold increase   #####################################


library(data.table)

df2 = data.frame(fread("/Users/kushaldey/Documents/singlecellLDSC/data/healthyvinflamed_score.csv"))
pval_logfold2 = apply(df2[,-1], 2, function(x) return(pnorm(x, 0, 1, lower.tail = F)))

#logfold[order(pval_logfold[,1], decreasing = F)[1:10], 1]

qq_logfold2 = apply(pval_logfold2, 2, function(y){
  z = -2*log(y+1e-08)
  PR = (z - min(z))/(max(z) - min(z))
  return(PR)
})


for(mm in 1:ncol(qq_logfold)){
  df = cbind.data.frame(df2[,1], qq_logfold2[, mm])
  write.table(df, file = paste0("/Users/kushaldey/Documents/singlecellLDSC/data/HealthyVInflamed_gene_score_Feb1/",
              colnames(qq_logfold)[mm], ".txt"), row.names = F, col.names = F, sep = "\t", quote=F)

  # genes = df[,1][order(qq_logfold[,mm], decreasing = T)[1:2000]]
  # df2 = cbind.data.frame(genes, 1)
  # write.table(df2, file = paste0("/Users/kushaldey/Documents/singlecellLDSC/data/Healthy_gene_score_Top10_Feb1/",
  #                               colnames(qq_logfold)[mm], ".txt"),
  #             row.names = F, col.names = F, sep = "\t", quote=F)

}

genes = df2[,1]
tab = cbind.data.frame(genes, 1)
write.table(tab, file = paste0("/Users/kushaldey/Documents/singlecellLDSC/data/All_IBD_genes.txt"),
            row.names = F, col.names = F, sep = "\t", quote=F)


