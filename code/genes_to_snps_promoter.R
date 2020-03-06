
library(data.table)

options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
print(args)
numchr = as.numeric(toString(args[1]))


df = data.frame(fread(paste0("/n/groups/price/kushal/Mouse_Humans/data/ABC9/",
                             "AllPredictions.AvgHiC.ABC0.015.minus150.withcolnames.ForABCPaper.txt.gz")))

cc = paste0(df$TargetGene, "::", df$TargetGeneTSS, "::", df$chr)
cc2 = unique(cc)
gene_names = as.character(sapply(cc2, function(x) return(strsplit(x, "::")[[1]][1])))
tss = as.numeric(sapply(cc2, function(x) return(strsplit(x, "::")[[1]][2])))
chr = as.character(sapply(cc2, function(x) return(strsplit(x, "::")[[1]][3])))


df2 = cbind.data.frame("Gene" = gene_names,
                       "chr" = chr,
                       "start" = tss - 5000,
                       "end" = tss + 5000)

base = data.frame(fread(paste0("/n/groups/price/kushal/Mouse_Humans/data/ANNOTATIONS/Baselines/",
                               "baselineLD_v2.1", "/", "baselineLD", ".", numchr, ".annot.gz")))


df_chr = df2[which(df2$chr == paste0("chr", numchr)), ]
gene_names = as.character(unique(df_chr$Gene))

for(gg in gene_names){
  temp = df_chr[which(df_chr$Gene == gg), ]
  start = temp$start
  end = temp$end
  idx =  which(base$BP > start & base$BP < end)
  rsids = base$SNP[idx]
  if(!dir.exists(paste0("/n/groups/price/kushal/Gene2SNP", "/", gg))){
    dir.create(paste0("/n/groups/price/kushal/Gene2SNP", "/", gg))
  }
  write.table(rsids, file = paste0("/n/groups/price/kushal/Gene2SNP", "/", gg, "/", "Promoter.txt"),
              row.names = F, col.names = F, sep = "\t", quote = F)

  cat("Finished processing gene:", gg, "\n")
}





