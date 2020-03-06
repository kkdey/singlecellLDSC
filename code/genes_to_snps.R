
library(data.table)

options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
print(args)
numchr = as.numeric(toString(args[1]))


#########################################   Genes to SNPs   ###################################################################

if(!dir.exists("/n/groups/price/kushal/Gene2SNP")){
  dir.create("/n/groups/price/kushal/Gene2SNP")
}

df = data.frame(fread(paste0("/n/groups/price/kushal/Mouse_Humans/data/ABC9/",
                             "AllPredictions.AvgHiC.ABC0.015.minus150.withcolnames.ForABCPaper.txt.gz")))

df2 = data.frame(fread(paste0("/n/groups/price/kushal/Mouse_Humans/data/ABC9/",
                              "AllPredictions.AvgHiC.ABC0.015.minus150.Blood.ForABCPaper.txt.gz")))

base = data.frame(fread(paste0("/n/groups/price/kushal/Mouse_Humans/data/ANNOTATIONS/Baselines/",
                                  "baselineLD_v2.1", "/", "baselineLD", ".", numchr, ".annot.gz")))


df_chr = df[which(df$chr == paste0("chr", numchr)), ]
gene_names = unique(df_chr$TargetGene)

for(gg in gene_names){
  temp = df_chr[which(df_chr$TargetGene == gg), ]
  starts = temp$start
  ends = temp$end
  idx2 = c()
  for(cc in 1:length(starts)){
    idx2 = c(idx2, which(base$BP > starts[cc] & base$BP < ends[cc]))
  }
  rsids = base$SNP[idx2]
  if(!dir.exists(paste0("/n/groups/price/kushal/Gene2SNP", "/", gg))){
    dir.create(paste0("/n/groups/price/kushal/Gene2SNP", "/", gg))
  }
  write.table(rsids, file = paste0("/n/groups/price/kushal/Gene2SNP", "/", gg, "/", "ABC.txt"),
              row.names = F, col.names = F, sep = "\t", quote = F)

  cat("Finished processing gene:", gg, "\n")
}

