library(data.table)
library(R.utils)
library(xgboost)

options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
print(args)
annot_dir <- toString(args[1])
annot_name <- toString(args[2])
output_dir <- toString(args[3])

####################  Merge S2G of modules with all genes   ######################################

#annot_dir="/n/groups/price/kushal/singlecellLDSC/data/ANNOTATIONS/Healthy_gene_score_Feb1"
#annot_name="Best4..Enterocytes"
#output_dir = "/n/groups/price/kushal/singlecellLDSC/data/ANNOTATIONS/Healthy_gene_score_Feb1_plus_All_IBD"
ll = list.files(paste0(annot_dir, "/", annot_name))


for(numchr in 1:22){
  for(mm in 1:length(ll)){
    df = data.frame(fread(paste0(annot_dir, "/",
                                 annot_name, "/",
                                 ll[mm], "/",
                                 ll[mm], ".", numchr, ".annot.gz")))
    temp = data.frame(fread(paste0("/n/groups/price/kushal/singlecellLDSC/data/ANNOTATIONS/All_genes/All_IBD_genes/",
                                   ll[mm], "/",
                                   ll[mm], ".", numchr, ".annot.gz")))
    newdf = cbind.data.frame(df[,1:4], temp[,5], df[,5])
    colnames(newdf) = c(colnames(df)[1:4], paste0("All_", ll[mm]), paste0("Healthy_", ll[mm]))

    if(!dir.exists(paste0(output_dir, "/", annot_name))){
      dir.create(paste0(output_dir, "/", annot_name))
    }
    if(!dir.exists(paste0(output_dir, "/", annot_name, "/", ll[mm]))){
      dir.create(paste0(output_dir, "/", annot_name, "/", ll[mm]))
    }
    write.table(newdf,
                file = gzfile(paste0(output_dir, "/",
                                     annot_name, "/",
                                     ll[mm], "/",
                                     ll[mm], ".",
                                     numchr, ".annot.gz")),
                quote=FALSE, row.names=FALSE)
  }
  cat("We are at chr:", numchr, "\n")
}


