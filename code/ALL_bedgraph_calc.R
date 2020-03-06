ABC_full_bedgraph_calc <- function(scores,
                              output_cell,
                              output_bed = "temp.bed"){

  df_pre = data.frame(fread(paste0("/n/groups/price/kushal/Mouse_Humans/data/ABC9/",
                                   "AllPredictions.AvgHiC.ABC0.015.minus150.withcolnames.ForABCPaper.txt.gz")))
  df = df_pre
  df2 = cbind.data.frame(df$chr, df$start, df$end, df$TargetGene)
  colnames(df2) = c("chr", "start", "end", "TargetGene")
  matched_ids = match(df2$TargetGene, names(scores))
  temp = as.numeric(scores)[matched_ids]
  temp[is.na(temp)] = 0
  final_bed = cbind(df2[,c(1:3)], temp)
  write.table(final_bed, paste0(output_cell, "/", output_bed),
              sep = "\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
}



Roadmap_full_bedgraph_calc = function(scores,
                                          output_cell,
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
#  geneanno = read.csv("/n/groups/price/kushal/Enhancer_MasterReg/data/geneanno_plus.txt", sep = "\t")
# # gene_df = read.delim("/n/groups/price/kushal/DiseaseNet/data/ensembl_to_gene_names.txt", header = F)
#  ff = geneanno$symbol[match(Enhancer[,4], geneanno$id)]
#  tmp = cbind.data.frame(Enhancer, ff, 1)
#  dff = tmp[,c(1:3, 5, 6)]
#  write.table(dff, file = "/n/groups/price/kushal/Mouse_Humans/data/R2G/Roadmap_Enhancers_All_Intergenic.txt",
#              col.names = F, row.names = F, quote=F, sep = "\t")


  Enhancer1 = data.frame(fread("/n/groups/price/kushal/Enhancer_MasterReg/data/R2G//Roadmap_Enhancers_All_Genic.txt",
                        header=F))
  Enhancer2 = data.frame(fread("/n/groups/price/kushal/Enhancer_MasterReg/data/R2G//Roadmap_Enhancers_All_Intergenic.txt",
                               header=F))
  Enhancer = rbind.data.frame(Enhancer1, Enhancer2)
  matched_ids = match(Enhancer[,4], names(scores))
  temp = as.numeric(scores)[matched_ids]
  temp[is.na(temp)] = 0

  df3 = cbind(Enhancer[,1:3], temp)

  write.table(df3, paste0(output_cell, "/", output_bed),
              sep = "\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
}



ABC_tissue_bedgraph_calc <- function(scores,
                              output_cell,
                              tissuename = "BRN",
                              output_bed = "temp.bed"){

  df_pre = data.frame(fread(paste0("/n/groups/price/kushal/Mouse_Humans/data/ABC9/",
                                   "AllPredictions.AvgHiC.ABC0.015.minus150.withcolnames.ForABCPaper.txt.gz")))
  df_pre = df_pre[which(df_pre$class == "intergenic" | df_pre$clas == "genic"), ]
  if(tissuename == "BRN"){
    tissuenames2 = c("neur", "Neur", "astro", "spinal", "Brain", "brain")
  }
  if(tissuename == "PANC"){
    tissuenames2 = c("pancrea")
  }
  if(tissuename == "BLD"){
    tissuenames2 = as.character(read.table("/n/groups/price/steven/h2gene/DATA/ABC8/ABC.listbloodQC.txt", header=F)[,1])
  }
  if(tissuename == "LNG"){
    tissuenames2 = c("lung")
  }

  tissue_ids = as.numeric(unlist(sapply(tissuenames2, function(x) return(grep(x, df_pre$CellType)))))
  df = df_pre[tissue_ids, ]
  df2 = cbind.data.frame(df$chr, df$start, df$end, df$TargetGene)
  colnames(df2) = c("chr", "start", "end", "TargetGene")
  matched_ids = match(df2$TargetGene, names(scores))
  temp = as.numeric(scores)[matched_ids]
  temp[is.na(temp)] = 0
  final_bed = cbind(df2[,c(1:3)], temp)
  write.table(final_bed, paste0(output_cell, "/", output_bed),
              sep = "\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
}



Roadmap_tissue_bedgraph_calc = function(scores,
                                      output_cell,
                                      tissuename = "BRN",
                                      output_bed = "temp.bed"){
  options(scipen = 999)
  # tissuename = "PANC"
  # tissuename = "BLD"
  # tissuename = "LNG"
  roadmap_meta = read.delim("/n/groups/price/kushal/DiseaseNet/data/Roadmap_map_EID_names.txt",
                            header=F)

  Road_ids =  unique(as.character(roadmap_meta[unlist(sapply(tissuename,
                                                             function(x) return(grep(x, roadmap_meta[,2])))), 1]))

  print(unique(as.character(roadmap_meta[unlist(sapply(tissuename,
                                                       function(x) return(grep(x, roadmap_meta[,2])))), 3])))

  Enhancer = c()
  for(ee in Road_ids){
    temp = read.table(paste0("/n/groups/price/kushal/DiseaseNet/data/RoadmapLinks/",
                            "links_", ee, "_7_2.5.txt"))
    temp2 = read.table(paste0("/n/groups/price/kushal/DiseaseNet/data/RoadmapLinks/",
                              "links_", ee, "_6_2.5.txt"))
    Enhancer = rbind(Enhancer, temp[,1:4], temp2[,1:4])
    cat("We processed file:", ee, "\n")
  }

  geneanno = read.csv("/n/groups/price/kushal/Enhancer_MasterReg/data/geneanno_plus.txt", sep = "\t")
  ff = geneanno$symbol[match(Enhancer[,4], geneanno$id)]
  tmp = cbind.data.frame(Enhancer, ff, 1)
  dff = tmp[,c(1:3, 5, 6)]

  matched_ids = match(dff[,4], names(scores))
  temp = as.numeric(scores)[matched_ids]
  temp[is.na(temp)] = 0

  outdf = cbind(dff[,1:3], temp)

  write.table(outdf, paste0(output_cell, "/", output_bed),
              sep = "\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

}
