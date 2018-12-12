#########################################
#
# NDA17 Release 1.1 General GFA Script
#
#########################################

library(GFA)

# support function for data chunks from:
# https://stackoverflow.com/questions/5577221/how-can-i-load-an-object-into-a-variable-name-that-i-specify-from-an-r-data-file

loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

# Directories with the data set:
# mydir <- paste0("/home/librad.laureateinstitute.org/mpaulus/ABCD_data/NDA17_1.1/")
mydir <- paste0("/Users/mpaulus/Dropbox (Personal)/Private/RDataAnalysis/ABCD_Data/NDA1.1/")
# mydir <- paste0("/Users/mpaulus/Dropbox (Personal)/Private/RDataAnalysis/ABCD_Data/NegReinforcement/")

# File components:
myfile <- c("nda17_1.1_")
GFAtext <- c("preGFA_")
dateext <- c("_11.30.2018")

# Steps:
# Read in Data Chunks

# Free Surfer Data
chunkfile <- paste0(mydir,myfile,"FS_cort_subcort",".RDATA")
MPPfs <- loadRData(chunkfile)
fsnames <- names(MPPfs)
# Cortical Thickness
thickvars <- fsnames[grep("thick_cort",fsnames)]
thickvars <- thickvars[-grep("mean",thickvars)]
# Sulcal Depth
sulcvars <- fsnames[grep("sulc_cort",fsnames)]
sulcvars <- sulcvars[-grep("mean",sulcvars)]
# Gray Matter Volume
volvars <- fsnames[grep("vol_cort",fsnames)]
volvars <- volvars[-grep("total",volvars)]
# Subcortical Volumes
subcortvars <- fsnames[grep("vol_subcort",fsnames)]
subcortvars <- subcortvars[-grep("volume",subcortvars)]
subcortvars <- subcortvars[-grep("wholebrain",subcortvars)]
subcortvars <- subcortvars[-grep("ventricles",subcortvars)]
subcortvars <- subcortvars[-grep("_csf",subcortvars)]

# MID variables
chunkfile <- paste0(mydir,myfile,"mid",".RDATA")
MPPmid <- loadRData(chunkfile)
fsnames <- names(MPPmid)
# reduces to only the beta values
fsnames <- fsnames[grep("beta",fsnames)]
# reduces to only anticipation
fsnames <- fsnames[grep("antic",fsnames)]
# reduces to use neutral as comparion
fsnames <- fsnames[grep("neutral",fsnames)]
# matches for both loss and gain
toMatch <- c("antic.loss","antic.reward")
fsnames <- fsnames[grepl(paste(toMatch, collapse="|"), fsnames)]
# uses the average of run1 and run2
fsnames <- fsnames[grep("_all_",fsnames)]
# divides into loss and reward variables
midrewardvars <- fsnames[grep(".reward.",fsnames)]
midlossvars <- fsnames[grep(".loss.",fsnames)]

# names(frame) <- sub(".*\\.", "", names(frame))

# Symptom Data - Parent
chunkfile <- paste0(mydir,myfile,"CBCL",".RDATA")
MPPcbcl <- loadRData(chunkfile)
fsnames <- names(MPPcbcl)
cbclvars <- fsnames[grep("cbcl_scr",fsnames)]


# Symptom Data - Youth
chunkfile <- paste0(mydir,myfile,"bis",".RDATA")
MPPbis <- loadRData(chunkfile)
fsnames <- names(MPPbis)
bisvars <- fsnames[grep("bis_ss",fsnames)]

chunkfile <- paste0(mydir,myfile,"upps",".RDATA")
MPPupps <- loadRData(chunkfile)
fsnames <- names(MPPupps)
uppsvars <- fsnames[grep("upps_ss",fsnames)]

# Cognition Data
chunkfile <- paste0(mydir,myfile,"cog",".RDATA")
MPPcog <- loadRData(chunkfile)
fsnames <- names(MPPcog)
cogvars <- fsnames[-grep("src_",fsnames)]


# Labels for different variable sets:
cbclabels <- c("Anxious/Depressed","Withdrawn","Somatic Sx","Social Problems","Thought Problems","Attention Problems","Rule Breaking","Aggressive Behavior","Internalizing","Externalizing","Total Problems")
bislabels <- c("BIS Total","BAS Reward Reactivity","BAS Drive","BAS Fun Seeking")
uppslabels <- c("Negative Urgency","Lack of Planning","Sensation Seeking","Positive Urgency","Lack of Perseverance")
coglabels <- c("Picture Vocabulary","Flanker Tes","List Sorting","Card Sorting","Pattern Comparison","Picture Sequence","Oral Reading Recog","Fluid Composite","Crystallized Composite","Cognition Total","RAVLT Short Delay","RAVLT Long Delay","WISC-V Matrix Reasoning")

# Get the labels from the desikan atlas csv file

# Desikan labels:
myall <- paste0(mydir,"ABCD_FS_labels",".csv")
desikan <- read.csv(myall,header=TRUE)

# Reorder variables to fit location and side and label:
thickvars_r <-thickvars[c(desikan$Order_orig)]
thicklabels <- paste(desikan[,2],desikan[,4],"TH",sep=" ")
sulcvars_r <- sulcvars[c(desikan$Order_orig)]
sulclabels <- paste(desikan[,2],desikan[,4],"SD",sep=" ")
volvars_r <- volvars[c(desikan$Order_orig)]
vollabels <- paste(desikan[,2],desikan[,4],"GV",sep=" ")

# Combine Data File
MPPall <- merge(MPPfs,MPPcog,by="src_subject_id",all=FALSE)
MPPall <- merge(MPPall,MPPcbcl,by="src_subject_id",all=FALSE)
MPPall <- merge(MPPall,MPPbis,by="src_subject_id",all=FALSE)
MPPall <- merge(MPPall,MPPupps,by="src_subject_id",all=FALSE)
MPPall <- merge(MPPall,MPPmid,by="src_subject_id",all=FALSE)

# Select Variables for GFA
fspass_sid <- MPPall[which(MPPall$fsqc_qc=='pass'),"src_subject_id"]

# Select the IDs for which we have good data and keep the data with subject
# that pass QC:
MPPmidall <- subset(MPPall,(MPPall$fsqc_qc=="pass") &
         (MPPall$tfmri_mid_beh_perform.flag==1) & 
         (MPPall$tfmri_mid_all_beta_dof>200) &
         (MPPall$tfmri_mid_all_sem_dof>200))

# Prep Data for GFA
indepvars <- c(thickvars,sulcvars,volvars,subcortvars,cogvars,cbclvars,bisvars,uppsvars)

#Original Variables
THICK <- as.matrix(MPPall[which(MPPall$fsqc_qc=='pass'),c(thickvars_r)])
SULC <- as.matrix(MPPall[which(MPPall$fsqc_qc=='pass'),c(sulcvars_r)])
VOL <- as.matrix(MPPall[which(MPPall$fsqc_qc=='pass'),c(volvars_r)])
SUBCORT <- as.matrix(MPPall[which(MPPall$fsqc_qc=='pass'),c(subcortvars)])

# Subtract the means of each subject:
rowmeans <- rowMeans(THICK)
dTHICK <- sweep(THICK,1,rowmeans)
rowmeans <- rowMeans(SULC)
dSULC <- sweep(SULC,1,rowmeans)
rowmeans <- rowMeans(VOL)  
dVOL <- sweep(VOL,1,rowmeans)
rowmeans <- rowMeans(SUBCORT)  
dSUBCORT <- sweep(SUBCORT,1,rowmeans)

# For the general GFA with sMRI
# COG <- as.matrix(MPPall[which(MPPall$fsqc_qc=='pass'),c(cogvars)])
# CBL <- as.matrix(MPPall[which(MPPall$fsqc_qc=='pass'),c(cbclvars)])
# BIS <- as.matrix(MPPall[which(MPPall$fsqc_qc=='pass'),c(bisvars)])
# UPPS <- as.matrix(MPPall[which(MPPall$fsqc_qc=='pass'),c(uppsvars)])

# MID variables:
MIDR <- as.matrix(MPPmidall[,c(midrewardvars)])
MIDL <- as.matrix(MPPmidall[,c(midlossvars)])

COG <- as.matrix(MPPmidall[,c(cogvars)])
CBL <- as.matrix(MPPmidall[,c(cbclvars)])
BIS <- as.matrix(MPPmidall[,c(bisvars)])
UPPS <- as.matrix(MPPmidall[,c(uppsvars)])

corM <- cor(MIDR, method='spearman', use='pair')
corrplot::corrplot(corM, title = "MID reward", type = "upper",method='ellipse',
                   addgrid.col = NA,tl.cex = .3,tl.col = "black",mar=c(0,0,1,0))

corM <- cor(MIDL, method='spearman', use='pair')
corrplot::corrplot(corM, title = "MID loss", type = "upper",method='ellipse',
                   addgrid.col = NA,tl.cex = .3,tl.col = "black",mar=c(0,0,1,0))

tmp <- cbind(MIDR,MIDL)
corM <- cor(tmp, method='spearman', use='pair')
corrplot::corrplot(corM, title = "MID reward/loss", method='ellipse',addgrid.col = NA, tl.cex = .3,tl.col = "black",mar=c(0,0,1,0))

# correlation plot of all variables (for testing purposes only)
#tmp <- cbind(THICK,SULC,VOL,SUBCORT,COG,CBL,BIS,UPPS)
#tmp <- cbind(dTHICK,dSULC,dVOL,dSUBCORT,COG,CBL,BIS,UPPS)
tmp <- cbind(MIDR,MIDL,COG,CBL,BIS,UPPS)

corM <- cor(tmp, method='spearman', use='pair')
corrplot::corrplot(corM, title = "Combined Variables", method='ellipse',addgrid.col = NA, tl.cex = .3,tl.col = "black",mar=c(0,0,1,0))

# form a list of variables from the data set
# MY <- list(THICK,SULC,VOL,SUBCORT,COG,CBL,BIS,UPPS)
# MY <- list(dTHICK,dSULC,dVOL,dSUBCORT,COG,CBL,BIS,UPPS)
MY <- list(MIDR,MIDL,COG,CBL,BIS,UPPS)

# Normalize all features
mynorm <- normalizeData(MY, type="scaleFeatures")

# set up the GFA defaults

# Get the default options
opts <- getDefaultOpts()
opts$vrbose <- 0

# number of data for posterior vector:
opts$iter.saved = 100

startK = length(c(indepvars))

# Run GFA
# GFAtext <- c("THICK_SULC_VOL_SUBCORT_COG_CBL_BIS_UPPS_")
# GFAtext <- c("dTHICK_dSULC_dVOL_dSUBCORT_COG_CBL_BIS_UPPS_")
GFAtext <- c("MIDR_MIDL_COG_CBL_BIS_UPPS_")
dateext <- c("_12.06.2018")
totalruns <- 1

set.seed(123);
res <- list()
for(i in 1:totalruns){
  print(i)
  myall <- paste0(mydir,myfile,GFAtext,i,dateext,".RData",sep="")
  print(myall)
  
  res[[i]] <- gfa(mynorm$train, K=startK, opts=opts)
  
  # Save as an interim variable
  myres <- res[[i]]
  
  # Save the GFA results
  # Write the result to a file
  myall <- paste0(mydir,myfile,GFAtext,i,dateext,".RData",sep="")
  save(myres, file=myall)
}

# Save robust components
rob <- robustComponents(res)
myall <- paste0(mydir,myfile,GFAtext,"Robust",dateext,".RData",sep="")
save(rob, file=myall)
