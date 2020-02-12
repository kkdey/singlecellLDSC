ABC_bedgraph_calc <- function(scores,
                              output_cell = "/n/groups/price/kushal/Mouse_Humans/data/BEDFILES/ABC_Bedgraphs2_MAX",
                              output_bed = "temp.bed"){

  filenames = paste0(as.character(read.table("/n/groups/price/steven/h2gene/DATA/ABC/ABC.listQC.txt")[,1]),
                     ".OriginalPeaks2.bed")

  bed1 = c()

  for(mm in 1:length(filenames)){
    if(file.info(paste0("/n/groups/price/steven/h2gene/DATA/ABC/abc/", filenames[mm]))$size > 0){
      df = data.frame(fread(paste0("/n/groups/price/steven/h2gene/DATA/ABC/abc/", filenames[mm])))
      df2 = df[which(df[,6] == "distal"),]
      common_genes = intersect(names(scores), df2[,4])
      idx = which(!is.na(match(df2[,4], common_genes)))
      df3 = df2[idx,]
      grades = scores[match(df3$V4, names(scores))]
      temp = cbind(df3, grades)
      bed1 = rbind(bed1, temp)
    }
    cat("Finished for file number :", mm, "\n")
  }

  bed2 = bed1[,c(1:3, 7)]
  write.table(bed2, paste0(output_cell, "/", output_bed),
              sep = "\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
}



Roadmap_Enhancer_bedgraph_calc = function(scores,
                                          output_cell = "/n/groups/price/kushal/Mouse_Humans/data/BEDFILES/ROAD_Enhancer_Bedgraphs",
                                          output_bed = "temp.bed"){
  options(scipen = 999)

#  ll = list.files("/n/groups/price/kushal/DiseaseNet/data/RoadmapLinks", pattern = "_7_2.5.txt")
#  Enhancer = c()
#  for(mm in 1:length(ll)){
#    temp =  read.table(paste0("/n/groups/price/kushal/DiseaseNet/data/RoadmapLinks/",
#                       ll[mm]), header=F)
#    Enhancer = rbind(Enhancer, temp[,1:4])
#    cat("We processed file:", mm, "\n")
#  }
#  gene_df = read.delim("/n/groups/price/kushal/DiseaseNet/data/ensembl_to_gene_names.txt", header = F)
#  ff = gene_df[match(Enhancer[,4], gene_df[,2]), 1]
#  tmp = cbind.data.frame(Enhancer, ff, 1)
#  dff = tmp[,c(1:3, 5, 6)]
#  write.table(dff, file = "/n/groups/price/kushal/Mouse_Humans/data/R2G/Roadmap_Enhancers_All.txt",
#              col.names = F, row.names = F, quote=F, sep = "\t")


  Enhancer = data.frame(fread("/n/groups/price/kushal/Mouse_Humans/data/R2G/Roadmap_Enhancers_All.txt",
                        header=F))
  matched_ids = match(Enhancer[,4], names(scores))
  temp = as.numeric(scores)[matched_ids]
  temp[is.na(temp)] = 0

  df3 = cbind(Enhancer[,1:3], temp)

  write.table(df3, paste0(output_cell, "/", output_bed),
              sep = "\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
}
