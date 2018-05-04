################################################################################################################################################
### the shapefit tool is released under the GNU General Public License (see http://www.gnu.org/licenses/gpl.html).
################################################################################################################################################
### whole trajectory aligned on first frame, ligand and protein coordinates extracted seperately
### for better comparison of one ligand in two different protein, alignment of both first frames and the according trajectories is suggested
### complementary POVME maps created for protein (cutoff 0.59 suggested) and "normal" POVME maps for the ligand (two different cutoffs e.g. 0.59 and 1.59) seperately with the modified POVME script (grid point distance of 0.5 suggested; with option: always same number of points per frame file; if other than the atoms parametrized in the POVME file are present in the protein or the ligand, those need to be parametrized there before map generation)
### for later pairing with other frame parameters, make sure ordering of the input data according to the MD timeline is maintained
#################################################################################################################################################
### inputs:
### 1.) extracted ligand coordinates (after alignment) over whole trajectory in pdb format
### 2.) starting point map used for POVME in pdb format
### 3.) complementary POVME map of the protein over the whole trajectory in pdb format
### 4.) two POVME ligand maps with different cutoffs over whole trajectory in pdb format
#################################################################################################################################################
### outputs:
###	1.) fixed width file of a protein-ligand overlap ratio (shapefit) over the whole trajectory -> "overlapRatio.dat"
###	2.) fixed width file of per (heavy) atom shapefit -> "overlapRatiobyatom.dat" (atom order as specified in "colnamesOverlapRatiobyatom.dat")
#################################################################################################################################################

workingDirectory <- "/some/working/directory"
setwd(workingDirectory)

###specification of input files

ligandPDB <- "pdbfile_of_ligandconfs/over_MD/without_headline.pdb" 
startingPoints <- "file_of_starting_point_map/as_created_by_POVME.pdb"
povmeOutputDir <- "directory/of/povme/output_files/"
filePatternProtein <- "pattern/according/to/povmeoutput"     #something like "trjproteinpart.+inLigOrProt.+\\.pdb" to specifiy the POVME output files of the points inside the protein 
filePatternLiginner <-  #something like "trjligandpart[0-9][0-9]inner059frame\\_[0-9]+\\.pdb" to specifiy the POVME output files of the points surrounding the ligand with distance to heavy atoms of at least 0.59 Angstrom
filePatternLigouter <-    #something like "trjligandpart[0-9][0-9]innerframe\\_[0-9]+\\.pdb" to specifiy the POVME output files of the points surrounding the ligand with distance to heavy atoms of at least 1.59 Angstrom

###load required packages 

library(doSNOW)
library(R.utils)
library(gdata)
library(RANN)


###load the ordering of atoms for the ligand 

ligatomfilelength <- countLines(ligandPDB)
ligatomlength <- (grep("END", readLines(ligandPDB))[[1]]-1)
numberofframes <- (ligatomfilelength)/(ligatomlength+1)
ligatomlist <- list()
for (i in 1:numberofframes){ligatomlist[[i]] <- read.fwf(ligandPDB, widths = c(-8,3,5,-14,8,8,8,-23,1), skip=((i-1)*(ligatomlength+1)), nrow = ligatomlength)}
saveRDS(ligatomlist,"ligatomlist.rds")

###load starting point map

filelength <- countLines(startingPoints)
pointMaplist <- read.fwf(startingPoints, widths = c(-30,8,8,8,-24), nrow = filelength)
pointMapmatrix <- do.call(cbind, pointMaplist)
saveRDS(pointMapmatrix,"pointMapmatrix.rds")

###load point maps of the protein in the binding area

filelistProt <- list.files(path = povmeOutputDir, pattern = filePatternProtein)
filepathlistProt <- paste(povmeOutputDir, filelistProt, sep="")

datalistProt <- list()
sizeProt <- list()
for (i in 1:length(filepathlistProt)){sizeProt[[i]] <- countLines(filepathlistProt[[i]])}

cluster <- makeCluster(3)
registerDoSNOW(cluster)

datalistProt <- foreach(x=1:length(filepathlistProt), .errorhandling="remove") %dopar% read.fwf(filepathlistProt[[x]], skip=1, widths=c(-30,8,8,8,-24), nrow=sizeProt[[x]]-2)

stopCluster(cluster)

saveRDS(datalistProt, "datalistProt.rds")

datalistnot0Prot <- list()
for (i in 1:length(datalistProt)) {datalistnot0Prot[[i]] <- datalistProt[[i]][apply(datalistProt[[i]], MARGIN=1,function(x) !all(x==0)), ]}
saveRDS(datalistnot0Prot,"datalistnot0Prot.rds")


###load ligand surrounding point maps with at least 1.59 Angstrom distance to the ligand

filelistInvlig <- list.files(path = povmeOutputDir, pattern = filePatternLigouter) 
filepathlistInvlig <- paste(povmeOutputDir, filelistInvlig, sep="")

datalistInvlig <- list()
sizeInvlig <- list()
for (i in 1:length(filepathlistInvlig)){sizeInvlig[[i]] <- countLines(filepathlistInvlig[[i]])}

cluster <- makeCluster(3)
registerDoSNOW(cluster)

datalistInvlig <- foreach(x=1:length(filepathlistInvlig), .errorhandling="remove") %dopar% read.fwf(filepathlistInvlig[[x]], skip=2, widths=c(-30,8,8,8,-24), nrow=sizeInvlig[[x]]-3)

stopCluster(cluster)

saveRDS(datalistInvlig, "datalistInvlig.rds")

datalistnot0Invlig <- list()
 for (i in 1:length(datalistInvlig)) {datalistnot0Invlig[[i]] <- datalistInvlig[[i]][apply(datalistInvlig[[i]], MARGIN=1,function(x) !all(x==0)), ]}
saveRDS(datalistnot0Invlig,"datalistnot0Invlig.rds")


###load ligand surrounding point maps with at least 0.59 Angstrom distance to the ligand

filelistInvliginner <- list.files(path = povmeOutputDir, pattern = filePatternLiginner) 
filepathlistInvliginner <- paste(povmeOutputDir, filelistInvliginner, sep="")

datalistInvliginner <- list()
sizeInvliginner <- list()
for (i in 1:length(filepathlistInvliginner)){sizeInvliginner[[i]] <- countLines(filepathlistInvliginner[[i]])}

cluster <- makeCluster(3)
registerDoSNOW(cluster)

datalistInvliginner <- foreach(x=1:length(filepathlistInvliginner), .errorhandling="remove") %dopar% read.fwf(filepathlistInvliginner[[x]], skip=2, widths=c(-30,8,8,8,-24), nrow=sizeInvliginner[[x]]-3)

stopCluster(cluster)

saveRDS(datalistInvliginner, "datalistInvliginner.rds")

datalistnot0Invliginner <- list()
 for (i in 1:length(datalistInvliginner)) {datalistnot0Invliginner[[i]] <- datalistInvliginner[[i]][apply(datalistInvliginner[[i]], MARGIN=1,function(x) !all(x==0)), ]}
saveRDS(datalistnot0Invliginner,"datalistnot0Invliginner.rds")

###from both ligand surrounding maps calculate the ligand surrounding surface map points (overlap)

notduplicatedPoints <- list()
datalistLigIDs <- list()
datalistLig <- list()
for (t in 1:length(datalistnot0Invliginner)){notduplicatedPoints[[t]] <- !duplicated(rbind (datalistnot0Invlig[[t]], datalistnot0Invliginner[[t]]))
		     mapStart <- dim(datalistnot0Invlig[[t]])[1]+1
                     mapEnd <- dim(rbind(datalistnot0Invlig[[t]], datalistnot0Invliginner[[t]]))[1]
                     datalistLigIDs[[t]] <- notduplicatedPoints[[t]][mapStart:mapEnd]
		     datalistLig[[t]] <- datalistnot0Invliginner[[t]][datalistLigIDs[[t]],]}


saveRDS(datalistLig, "datalistLig.rds")

###assign ligand atoms to near parts of the ligand surrounding point map

ligatomlistnoHs <- list()

for (l in 1:length(ligatomlist)){ligatomlistnoHs[[l]] <- ligatomlist[[l]][!(ligatomlist[[l]]["V6"]=="H"),]}


nnresult<-list()

for (l in 1:length(ligatomlist)){
	nnresult[[l]] <- nn2(ligatomlistnoHs[[l]][,3:5], datalistLig[[l]][,1:3], k=1, treetype="kd")
}

resultlist <- list()

for (l in 1:length(ligatomlist)){
results <- list()
	for (n in 1:length(nnresult[[l]][[1]])){
		results[[n]] <- ligatomlistnoHs[[1]][nnresult[[l]][[1]][n],2]
	}
resultlist[[l]] <- unlist(results)
}

for (l in 1:length(ligatomlist)){
datalistLig[[l]]["atom"] <- resultlist[[l]]
}

saveRDS(datalistLig, "datalistLigAtom.rds")

###find overlap of ligand surrounding and protein surrounding maps

duplicatedPoints <- list()
overlapMapIDs <- list()
overlapMap <- list()
for (t in 1:length(datalistnot0Prot)){duplicatedPoints[[t]] <- duplicated(rbind (datalistnot0Prot[[t]], datalistLig[[t]][1:3]))
                     mapStart <- dim(datalistnot0Prot[[t]])[[1]]+1
                     mapEnd <- length(duplicatedPoints[[t]])
                     overlapMapIDs[[t]] <- duplicatedPoints[[t]][mapStart:mapEnd]
		     overlapMap[[t]] <- datalistLig[[t]][overlapMapIDs[[t]],]}


saveRDS(overlapMap, "overlapMap.rds")

overlapLength <- list()
for (i in 1:length(overlapMap)){overlapLength[[i]] <- (dim(overlapMap[[i]])[[1]])}
saveRDS(overlapLength, "overlapLength.rds")

write.fwf((as.data.frame(unlist(overlapLength))), file="overlapLength.dat",rownames=FALSE, colnames=FALSE, quote=FALSE, justify="right", width=10)

overlapRatio <- list()
for (i in 1:length(overlapLength)){overlapRatio[[i]] <- (overlapLength[[i]]/(dim(datalistLig[[i]])[1]))}
saveRDS(overlapRatio, "overlapRatio.rds")

write.fwf((as.data.frame(unlist(overlapRatio))), file="overlapRatio.dat",rownames=FALSE, colnames=FALSE, quote=FALSE, justify="right", width=12)

ratiobyatom <- list()
for(l in 1:length(ligatomlist)){
ratiobyatom[[l]] <- table(overlapMap[[l]]$atom)/table(datalistLig[[l]]$atom)
}
saveRDS(ratiobyatom, "overlapRatiobyatom.rds")


df <- lapply(ratiobyatom,as.data.frame)
tdf <- lapply(df,t)
tdf2 <- list()
for (l in 1:length(tdf)){
tdf2[[l]] <- as.numeric(tdf[[l]][2,])}
exp <- do.call(rbind,tdf2)
saveRDS(exp, "plotdata.rds")
names <- unlist(df[[1]][1])

write.fwf(exp, file="overlapRatiobyatom.dat",rownames=FALSE, colnames=FALSE, quote=FALSE, justify="right", width=12)
write.fwf(t(df[[1]][1]), file="colnamesOverlapRatiobyatom.dat",rownames=FALSE, colnames=FALSE, quote=FALSE, justify="right", width=12)

