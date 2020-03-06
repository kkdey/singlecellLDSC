

folder = "/n/groups/price/kushal/singlecellLDSC/data/Gene_Scores/Modules/healthy/celltype_enriched/colon"
ll = list.files(folder)
geneanno = read.csv("/n/groups/price/kushal/Enhancer_MasterReg/data/geneanno_plus.txt", sep = "\t")

size_vec = c()
for(mm in 1:length(ll)){
  temp = data.frame(fread(paste0(folder, "/", ll[mm])))
  matched_ids = match(temp[,1], geneanno$symbol)
  temp2 = temp[which(!is.na(matched_ids)), ]
  size_vec = c(size_vec, mean(temp2[,2]))
}
names(size_vec) = ll


