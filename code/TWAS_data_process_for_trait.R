options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
print(args)
cellname <- toString(args[1])
traitname <- toString(args[2])

################################  TWAS analysis per trait   #####################################################


require(data.table)
geneanno = data.table::fread("/n/groups/price/kushal/ExPecto/resources/geneanno.csv")
gene_names = geneanno$symbol


dat2 = data.frame(fread(paste0(cellname, "/", traitname)))
dat2$TWAS.Z.ABS = abs(dat2$TWAS.Z)
tab = xtabs(TWAS.Z.ABS~ID+PANEL,aggregate(TWAS.Z.ABS~ID+PANEL,dat2,max))
tab_frame = as.data.frame.matrix(tab)
write.table(tab_frame, file = gzfile(paste0(cellname, "/Zmats/", traitname, ".txt.gz")),
            sep = "\t", quote=F)

tab2 = xtabs(HSQ~ID+PANEL,aggregate(HSQ~ID+PANEL,dat2,max))
tab2_frame = as.data.frame.matrix(tab2)
write.table(tab2_frame, file = gzfile(paste0(cellname, "/HSQ/", traitname, ".txt.gz")),
            sep = "\t", quote=F)


#dat2 = data.frame(fread(paste0("/n/groups/price/kushal/singlecellLDSC/data/TWAS-Hub/CD_deLange2017.dat"))
