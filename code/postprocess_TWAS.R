
##############################  Postprocess  TWAS   ###################################################

ll = list.files('/n/groups/price/kushal/singlecellLDSC/data/TWAS-Hub/HSQ')

list_genes = vector(mode = 'list', length(ll))
for(mm in 1:length(ll)){
  temp = data.frame(data.table::fread(paste0('/n/groups/price/kushal/singlecellLDSC/data/TWAS-Hub/HSQ/', ll[mm])))
  temp_mat = temp[,-1]
  rownames(temp_mat) = temp[,1]
  list_genes[[mm]] = temp[,1]
  cat("we are at tissue:", mm, "\n")
}

all_genes = Reduce(union, list_genes)

final_df = matrix(0, length(all_genes), 74*length(ll))

for(mm in 1:length(ll)){
  temp = data.frame(data.table::fread(paste0('/n/groups/price/kushal/singlecellLDSC/data/TWAS-Hub/HSQ/', ll[mm])))
  temp_mat = temp[,-1]
  rownames(temp_mat) = temp[,1]
  final_df[match(rownames(temp_mat), all_genes), ((mm-1)*74+1):(mm*74)] = as.matrix(temp_mat)
  cat("we are at tissue:", mm, "\n")
}
rownames(final_df) = all_genes
colnames(final_df) = apply(expand.grid(colnames(temp_mat), ll), 1, paste, collapse=".data.")
fwrite(data.frame(final_df), file = '/n/groups/price/kushal/singlecellLDSC/data/TWAS-Hub/Genes_by_X_TWAS_HSQ_Celltypes.txt',
       quote=F, sep = "\t", row.names=T, col.names=T)


final_df = matrix(0, length(all_genes), length(ll))
for(mm in 1:length(ll)){
  temp = data.frame(data.table::fread(paste0('/n/groups/price/kushal/singlecellLDSC/data/TWAS-Hub/HSQ/', ll[mm])))
  vec = apply(as.matrix(temp[,-1]), 1, max)
  final_df[match(temp[,1], all_genes), mm] = vec
  cat("we are at tissue:", mm, "\n")
}
rownames(final_df) = all_genes
colnames(final_df) = ll
fwrite(data.frame(final_df), file = '/n/groups/price/kushal/singlecellLDSC/data/TWAS-Hub/Genes_by_X_TWAS_HSQ.txt',
       quote=F, sep = "\t", row.names=T, col.names=T)





cc = data.frame(data.table::fread('/n/groups/price/kushal/singlecellLDSC/data/TWAS-Hub/Genes_by_X_TWAS_HSQ_Celltypes.txt'))


apply(expand.grid(colnames(temp_mat), ll), 1, paste, collapse=".data.")


