library(data.table)
library(R.utils)
library(xgboost)

options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
print(args)
annot_dir <- toString(args[1])
genescore_dir <- toString(args[2])
bed_dir <- toString(args[3])
annot_name <- toString(args[4])

score_file = paste0(genescore_dir, "/", annot_name, ".txt")
gene_scores = read.delim(score_file, header=F)

if(!dir.exists(paste0(annot_dir, "/", annot_name))){
  dir.create(paste0(annot_dir, "/", annot_name))
}

if(!dir.exists(paste0(bed_dir, "/", annot_name))){
  dir.create(paste0(bed_dir, "/", annot_name))
}

source("/n/groups/price/kushal/singlecellLDSC/code/ALL_bedgraph_calc.R")

gene_scores = read.delim(score_file, header=F)

scores = gene_scores[,2]
names(scores) = gene_scores[,1]

out4 = Roadmap_Enhancer_bedgraph_calc(scores,
                                      output_cell = paste0(bed_dir, "/", annot_name),
                                      output_bed = "Roadmap_Enhancer.bed")

df = read.table("/n/groups/price/kushal/Mouse_Humans/data/Gene_5kb.txt", header=T)
df[which(df[,2] < 0), 2] = 0
score = gene_scores[match(df$gene, gene_scores[,1]), 2]
score[is.na(score)] = 0

temp = cbind.data.frame(df[,1:3], score)
temp2 = temp[which(temp$score != 0),]
write.table(temp2, file = paste0(bed_dir, "/", annot_name, "/", "5kb.bed"),
            quote=F, sep = "\t", row.names=F, col.names=F)


df = read.table("/n/groups/price/kushal/Mouse_Humans/data/Gene_100kb.txt", header=T)
df[which(df[,2] < 0), 2] = 0
score = gene_scores[match(df$gene, gene_scores[,1]), 2]
score[is.na(score)] = 0

temp = cbind.data.frame(df[,1:3], score)
temp2 = temp[which(temp$score != 0),]
write.table(temp2, file = paste0(bed_dir, "/", annot_name, "/", "100kb.bed"),
            quote=F, sep = "\t", row.names=F, col.names=F)


