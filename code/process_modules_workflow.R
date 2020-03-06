tissuename = "brain"
input_cell="/n/groups/price/kushal/singlecellLDSC/Modules/healthy/celltypeenriched"
output_cell="/n/groups/price/kushal/singlecellLDSC/data/Gene_Scores/Modules/healthy/celltype_enriched"
filename=paste0(input_cell, "/", tissuename, "_score.csv")
pfilename = paste0(input_cell, "/", tissuename, "_pval.csv")

library(data.table)

df2 = data.frame(fread(filename))
pdf2 = data.frame(fread(pfilename))
temp = df2[,-1]
temp[temp < 0] = 0

pval_logfold2 = apply(temp, 2, function(x) return(2*pnorm(abs(x), 0, 1, lower.tail = F)))

qq_logfold2 = apply(pval_logfold2, 2, function(y){
  z = -2*log(y+1e-08)
  PR = (z - min(z))/(max(z) - min(z))
  return(PR)
})


if(!dir.exists(paste0(output_cell, "/", tissuename))){
  dir.create(paste0(output_cell, "/", tissuename))
}

for(mm in 1:ncol(qq_logfold2)){
  df = cbind.data.frame(df2[,1], qq_logfold2[, mm])
  write.table(df, file = paste0(output_cell, "/", tissuename,  "/",
                                colnames(qq_logfold2)[mm], ".txt"), row.names = F, col.names = F, sep = "\t", quote=F)

}

mm=1
tissuename="brain"
tt = read.delim(paste0(output_cell, "/", tissuename,  "/",
                         colnames(qq_logfold2)[mm], ".txt"), header=F)
