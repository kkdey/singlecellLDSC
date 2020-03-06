
############################  Gene sets for pathway analysis  ######################################################

library(data.table)
input_cell = "/Users/kushaldey/Documents/singlecellLDSC/data/Modules/healthy/celltypeenriched"
output_cell = "/Users/kushaldey/Documents/singlecellLDSC/data/pathdb_brain_healthy2"
tissuename = "brain"
filename=paste0(input_cell, "/", tissuename, "_score.csv")
pfilename = paste0(input_cell, "/", tissuename, "_pval.csv")
pvaldf = data.frame(fread(pfilename))

library(data.table)

df2 = data.frame(fread(filename))
pdf2 = data.frame(fread(pfilename))
temp = df2[,-1]
temp[temp < 0] = 0

pval_logfold2 = apply(temp, 2, function(x) return(2*pnorm(abs(x), 0, 1, lower.tail = F)))


for(mm in 2:ncol(pvaldf)){
  idx = which(pval_logfold2[,(mm-1)] < 1e-08)
  genes = pvaldf[idx, 1]
  name = colnames(pvaldf)[mm]
  write.table(genes, file = paste0(output_cell, "/", name, ".txt"), col.names = F, row.names = F, sep = "\t", quote=F)
}

mm=1
tissuename = "brain"
output_cell="/n/groups/price/kushal/singlecellLDSC/data/Gene_Scores/Modules/healthy/celltype_enriched"
temp = "/n/groups/price/kushal/singlecellLDSC/data/pathdb_brain_healthy3"
files = list.files(paste0(output_cell, "/", tissuename), full.names = T)
files2 = list.files(paste0(output_cell, "/", tissuename))


for(mm in 1:length(files)){
  tt = read.delim(files[mm], header=F)
  idx = which(tt[,2] > 0.8)
  genes = tt[idx, 1]
  write.table(genes, paste0(temp, "/", files2[mm]), col.names = F, row.names = F, sep = "\t", quote=F)
}

temp1 = read.delim("/n/groups/price/kushal/singlecellLDSC/data/Gene_Scores/Modules/healthy/celltype_enriched/brain/LAMP5.txt",
                   header=F)
temp2 = read.delim("/n/groups/price/kushal/singlecellLDSC/data/pathdb_brain_healthy3/LAMP5.txt", header=F)

temp1[match(temp2[,1], temp1[,1]), ]

length(intersect(temp1[order(temp1[,2], decreasing = T)[1:2000], 1], temp2[,1]))



