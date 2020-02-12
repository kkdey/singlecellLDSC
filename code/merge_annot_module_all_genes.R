

library(data.table)
library(R.utils)
library(xgboost)

options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
print(args)
annot_dir <- toString(args[1])
annot_name <- toString(args[4])

####################  Merge S2G of modules with all genes   ######################################

ll = list.files("/n/groups/price/kushal/singlecellLDSC/data/ANNOTATIONS/All_genes/All_IBD_genes/100kb")
annot_dir="/n/groups/price/kushal/singlecellLDSC/data/ANNOTATIONS/Healthy_gene_score_Feb1"
annot_name="Best4..Enterocytes"

for(numchr in 1:22){
  for(mm in 1:length(ll)){
    df =
  }
}
