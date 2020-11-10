###################################################################################
################# Mass Univariate Analysis for UKBB data ##########################
###################################################################################
###################################################################################
rm(list = ls(all = TRUE))  
#start with an empty environment

setwd("~/cardiac/Experiments_of_Maria/HCM_expressivity/statistical_analysis/")
library(data.table)
library(multtest)
library(mutools3D)
library(RcppArmadillo) # required in mulitply.cpp
library(RcppEigen) # required in mulitply.cpp
library(Rcpp) # required in permFL_fast
library(float) # required in permFL_fast

func_dir<-"~/cardiac/Experiments_of_Maria/HCM_expressivity/statistical_analysis/"
#folder directory of the code. Will be used for the functions as: source(paste(dir,"functions/xxx", sep="/")) 
#in the permFL_fast function.

sourceCpp(paste(func_dir,"functions/multiply.cpp",sep=""))
#include a cpp functions for faster multiplication of matrices

# source(paste(func_dir,"functions/murq.R",sep=""))
#include the mur function with QR decomposition (slightly faster!), HC4m = FALSE
#or include the mur function from mutools3D package

#include the functions for TFCE from mutools3D package

source(paste(func_dir,"functions/permFL_fast.R",sep=""))
#include the permFL_fast functions with murq and HC4m = FALSE

#CLINICAL DATA MATRIX
#NCOL = N COVARIATES UNDER STUDY
#NROW = N PATIENTS

func.dir<-"~/cardiac/DL_segmentation/3D_HCM_2020_Experiments/Mass_univariate_analysis/"

inputClinical <- readRDS(paste(func.dir,"data/WTedLVclinicaldata.rds", sep=""))
#inputClinical<-as.data.frame(read.csv(paste(func.dir,"data/Phenotypes_13k_imp.csv", sep="")))
summary_table<-as.data.frame(read.table(paste(func.dir,"data/summarytable_meta_2_UKB_HCM.txt", sep=""), header=T))
# HCM summary table
head(summary_table)


pos1<-match(inputClinical$eid_18545,summary_table$eid_18545)
pos1<-na.omit(pos1)
summary_tablenew<-summary_table[pos1,]

pos<-match(summary_tablenew$eid_18545,inputClinical$eid_18545)
pos<-na.omit(pos)
inputClinical<-inputClinical[pos,c(1:6)]

nr1<-which(inputClinical$Race==-1) 
inputClinical<-inputClinical[-nr1,] # remove -1 (Do not know) from Race
nr3<-which(inputClinical$Race==-3) 
inputClinical<-inputClinical[-nr3,] # remove -3 (Prefer not to say) from Race

summary_tablenew<-summary_tablenew[-nr1,]
summary_tablenew<-summary_tablenew[-nr3,]

inputClinical$Age <- scale(inputClinical$Age)
inputClinical$BSA <- scale(inputClinical$BSA)
inputClinical$SBP <- scale(inputClinical$SBP)
#inputClinical$DBP <- scale(inputClinical$DBP)

head(inputClinical)
X <- data.matrix(inputClinical[,-c(1)])
head(X)
X <- cbind(1, X)
#X<-X[,c(1:6)]
dim(X)
head(X)

## HCM steps
setwd("~/cardiac/DL_segmentation/3D_HCM_2020_Experiments/Mass_univariate_analysis/output/HCM")
l<-c("step1", "step2a","step2a_only","step2b","step2b_only","step3a","step3a_only","step3b","step3b_only","step3ab","step4","step4_13", "step3ab_min_step2b")
# create folders to save the output for each step. No need to do this at each step! 
for (i in 1:length(l)){
  dir.create(l[i])}

#step 0 vs 1
{
  p0<-which(summary_tablenew$HCMstep0==1)
  p1<-which(summary_tablenew$HCMstep1==1)
  
  pcomb<-c(p0,p1)
  pcomb<-unique(pcomb)
  summary_new1<-summary_tablenew[pcomb,]
  ps<-match(summary_new1$eid_18545,inputClinical$eid_18545)
  
  Xn<-as.data.frame(cbind(X[ps,], summary_new1$HCMstep1))
  Xn<-as.matrix(Xn)
  colnames(Xn)[8]<-c("HCM")
  length(which(Xn[,8]==1))
  setwd(paste(func.dir,"output/HCM/step1", sep=""))
}
#step 0 vs 2a
{
  p0<-which(summary_tablenew$HCMstep0==1)
  p2a<-which(summary_tablenew$HCMstep2A==1)
  p2b<-which(summary_tablenew$HCMstep2B==1)
    p2a <-p2a[!(p2a%in% p2b)]
  
  pcomb<-c(p0,p2a)
  summary_new2a<-summary_tablenew[pcomb,]
  ps<-match(summary_new2a$eid_18545,inputClinical$eid_18545)
  s<-which(is.na(ps))
  ps<-na.omit(ps)
  
  Xn<-as.data.frame(cbind(X[ps,], summary_new2a$HCMstep2A))
  Xn<-as.matrix(Xn)
  colnames(Xn)[8]<-c("HCM")
  length(which(Xn[,8]==1))
  setwd(paste(func.dir,"output/HCM/step2a", sep=""))
}
#step 0 vs 2b
{
  p0<-which(summary_tablenew$HCMstep0==1)
  p2b<-which(summary_tablenew$HCMstep2B==1)
  p3a<-which(summary_tablenew$HCMstep3A==1)
  p3b<-which(summary_tablenew$HCMstep3B==1)
  p3ab<-c(p3a,p3b)
  p2b <-p2b[!(p2b%in% p3ab)]
  
  pcomb<-c(p0,p2b)
  summary_new2b<-summary_tablenew[pcomb,]
  ps<-match(summary_new2b$eid_18545,inputClinical$eid_18545)
  s<-which(is.na(ps))
  ps<-na.omit(ps)
  
  Xn<-as.data.frame(cbind(X[ps,], summary_new2b$HCMstep2B))
  Xn<-as.matrix(Xn)
  colnames(Xn)[7]<-c("HCM")
  length(which(Xn[,7]==1))
  setwd(paste(func.dir,"output/HCM/step2b_only", sep=""))
}
#step 0 vs 3a
{
  p0<-which(summary_tablenew$HCMstep0==1)
  p3a<-which(summary_tablenew$HCMstep3A==1)
  p3b<-which(summary_tablenew$HCMstep3B==1)
  p3a <-p3a[!(p3a%in% p3b)]
  
  pcomb<-c(p0,p3a)
  summary_new3a<-summary_tablenew[pcomb,]
  ps<-match(summary_new3a$eid_18545,inputClinical$eid_18545)
  s<-which(is.na(ps))
  ps<-na.omit(ps)
  
  Xn<-as.data.frame(cbind(X[ps,], summary_new3a$HCMstep3A))
  Xn<-as.matrix(Xn)
  colnames(Xn)[8]<-c("HCM")
  length(which(Xn[,8]==1))
  setwd(paste(func.dir,"output/HCM/step3a_only", sep=""))
}
#step 0 vs 3b
{
  p0<-which(summary_tablenew$HCMstep0==1 )
  p3b<-which(summary_tablenew$HCMstep3B==1)
  p3a<-which(summary_tablenew$HCMstep3A==1)
  p3b <-p3b[!(p3b%in% p3a)]
  
  
  pcomb<-c(p0,p3b)
  summary_new3b<-summary_tablenew[pcomb,]
  ps<-match(summary_new3b$eid_18545,inputClinical$eid_18545)
  s<-which(is.na(ps))
  ps<-na.omit(ps)
  
  Xn<-as.data.frame(cbind(X[ps,], summary_new3b$HCMstep3B))
  Xn<-as.matrix(Xn)
  colnames(Xn)[8]<-c("HCM")
  length(which(Xn[,8]==1))
  setwd(paste(func.dir,"output/HCM/step3b_only", sep=""))
}
#step 0 vs 4
{
  p0<-which(summary_tablenew$HCMstep0==1 )
  p4<-which(summary_tablenew$HCMstep4==1)
  hasHCM<-which(summary_tablenew$hasHCM_k == "TRUE")
  p4 <-p4[!(p4%in% hasHCM)]
  
  
  pcomb<-c(p0,p4)
  summary_new4<-summary_tablenew[pcomb,]
  ps<-match(summary_new4$eid_18545,inputClinical$eid_18545)
  s<-which(is.na(ps))
  ps<-na.omit(ps)
  
  Xn<-as.data.frame(cbind(X[ps,], summary_new4$HCMstep4))
  Xn<-as.matrix(Xn)
  length(which(Xn[,7]==1))
  colnames(Xn)[7]<-c("HCM")
  setwd(paste(func.dir,"output/HCM/step4_13", sep=""))
}
#step 0 vs 3a and 3b
{
  p0<-which(summary_tablenew$HCMstep0==1)
  p3a<-which(summary_tablenew$HCMstep3A==1)
  p3b<-which(summary_tablenew$HCMstep3B==1)
  
  p3ab<-c(p3a,p3b)
  p3ab<-unique(p3ab)
  p3ab<-p3ab[order(p3ab)]

  pcomb<-c(p0,p3ab)
  pcomb<-unique(pcomb)
  
  merge_col3ab<-as.integer(ifelse(summary_tablenew$HCMstep3A==1|summary_tablenew$HCMstep3B==1,1,0))
  summary_tablenew$HCMstep3A_3B<-merge_col3ab
  summary_new3ab<-summary_tablenew[pcomb,]
  length(which(summary_new3ab$HCMstep3A_3B==1))
  ps<-match(summary_new3ab$eid_18545,inputClinical$eid_18545)
  s<-which(is.na(ps))
  ps<-na.omit(ps)
  
  Xn<-as.data.frame(cbind(X[ps,], summary_new3ab$HCMstep3A_3B))
  Xn<-as.matrix(Xn)
  length(which(Xn[,8]==1))
  colnames(Xn)[8]<-c("HCM")
  setwd(paste(func.dir,"output/HCM/step3ab", sep=""))
}

##step 3A+B minus 2B
{
  p0<-which(summary_tablenew$HCMstep0==1)
  p2b<-which(summary_tablenew$HCMstep2B==1)
  p3a<-which(summary_tablenew$HCMstep3A==1)
  p3b<-which(summary_tablenew$HCMstep3B==1)
  
  p3ab<-c(p3a,p3b)
  p3ab<-unique(p3ab)
  p3ab<-p3ab[order(p3ab)]
  
  length(p2b)
  p2b <-p2b[!(p2b%in% p3ab)]
  length(p2b)
  
  pcomb<-c(p0,p2b)
  summary_new2b<-summary_tablenew[pcomb,]
  ps<-match(summary_new2b$eid_18545,inputClinical$eid_18545)
  s<-which(is.na(ps))
  ps<-na.omit(ps)
  
  Xn<-as.data.frame(cbind(X[ps,], summary_new2b$HCMstep2B))
  Xn<-as.matrix(Xn)
  colnames(Xn)[8]<-c("HCM")
  length(which(Xn[,8]==1))
  setwd(paste(func.dir,"output/HCM/step3ab_min_step2b", sep=""))
}


# These lines are for producing the meshCoordinates.txt file for this dataset
Xcoord <- readRDS(paste(func.dir,"data/WTedXcoordinate.rds", sep=""))
Ycoord <- readRDS(paste(func.dir,"data/WTedYcoordinate.rds", sep=""))
Zcoord <- readRDS(paste(func.dir,"data/WTedZcoordinate.rds", sep=""))

{
  Xcoordd<-Xcoord[pos,]
  Xcoordd<-Xcoordd[-nr1,]
  Xcoordd<-Xcoordd[-nr3,]

  Xcoordd<-Xcoordd[ps,]
  dim(Xcoordd)

  Ycoordd<-Ycoord[pos,]
  Ycoordd<-Ycoordd[-nr1,]
  Ycoordd<-Ycoordd[-nr3,]

  Ycoordd<-Ycoordd[ps,]
  dim(Ycoordd)

  Zcoordd<-Zcoord[pos,]
  Zcoordd<-Zcoordd[-nr1,]
  Zcoordd<-Zcoordd[-nr3,]

  Zcoordd<-Zcoordd[ps,]
  dim(Zcoordd)
}

#IMAGING DATA MATRIX
#NCOL = N POINTS ON THE ATLAS
#NROW = N PATIENTS
Yw <- readRDS(paste(func.dir,"data/WTedLV.rds", sep=""))
#DATA PRE-PROCESSING

dim(Yw)
Ywd<-Yw[pos,]
Ywd<-Ywd[-nr1,] # remove -1 (Do not know) from Race
Ywd<-Ywd[-nr3,] # remove -3 (Prefer not to say) from Race

#Ywd<-Yw
Ywd<-Ywd[ps,]
Ywd <- scale(Ywd)
dim(Ywd)

#IMAGING DATA MATRIX
#NCOL = N POINTS ON THE ATLAS
#NROW = N PATIENTS
Ys <- readRDS(paste(func.dir,"data/S2SedLV.rds", sep=""))
#DATA PRE-PROCESSING

Ysd<-Ys[pos,]
Ysd<-Ysd[-nr1,] # remove -1 (Do not know) from Race
Ysd<-Ysd[-nr3,] # remove -3 (Prefer not to say) from Race

Ysd<-Ysd[ps,]
Ysd <- scale(Ysd)
dim(Ysd)

#NUMBER OF CORES TO USE
# All core detected here
nofCores = 48

# 1000 permutations
nPermutations=1000

whichEE <- 3
#1 endo, 2 epi, 3 full shape
endoEpi <- read.table(paste(func.dir,"data/endo_epi.txt", sep=""))
vert2print <- list(which(endoEpi[,4]==0),which(endoEpi[,4]==1),1:length(endoEpi[,4]))
Ywd<-Ywd[,vert2print[[whichEE]]]
Ysd<-Ysd[,vert2print[[whichEE]]]

##READ THE NNLIST ASSOCIATED TO EACH VERTEX
NNmatrix <- readRDS(paste(func.dir,"data/redNNmatrix.rds", sep=""))
##READ AREAS ASSOCIATED TO EACH VERTEX
A <- readRDS(paste(func.dir,"data/LVarea.rds", sep=""))

## Produce a new mean of meshCoordinates

# x_Coordinates <- Xcoordd[,vert2print[[whichEE]]]
# y_Coordinates <- Ycoordd[,vert2print[[whichEE]]]
# z_Coordinates <- Zcoordd[,vert2print[[whichEE]]]
# 
# X_coordinates <- colMeans(x_Coordinates)
# Y_coordinates <- colMeans(y_Coordinates)
# Z_coordinates <- colMeans(z_Coordinates)
# mesh_Coordinates <- data.frame(X_coordinates, Y_coordinates, Z_coordinates)
# colnames(mesh_Coordinates) <- c("x", "y", "z")

## or use the meshCoordinates from the mean of 13k I have created
mesh_Coordinates <- read.table(paste(func.dir,"data/meshCoordinates_13k.txt", sep=""), quote="\"", comment.char="", header=T)
colnames(mesh_Coordinates) <- c("x", "y", "z")


extractNames<-colnames(Xn)
length(extractNames)
extractNames
extract=7
{
  resultw <- murq(Xn,Ywd,extract)
  results <- mur(Xn,Ysd,extract)
  #extract only the pheno betas
  
  #MULTIPLE TESTING CORRECTION
  correctedw <- mt.rawp2adjp(resultw[,3], proc=c("BH"), na.rm = FALSE)
  pvalueADJ5tsbh <- array(dim = length(resultw[,3]))
  BHpvaluesw <- correctedw$adjp[order(correctedw$index),][,2]
  correcteds <- mt.rawp2adjp(results[,3], proc=c("BH"), na.rm = FALSE)
  pvalueADJ5tsbh <- array(dim = length(results[,3]))
  BHpvaluess <- correcteds$adjp[order(correcteds$index),][,2]
  
  meshCoordinates <- cbind(mesh_Coordinates,99999)
  #PRINT OUTPUT
  meshCoordinates[vert2print[[whichEE]],4] <- resultw[,1]
  write.table(meshCoordinates, paste(extractNames[iSNP],"_beta_WT.txt",sep=""), col.names = FALSE, row.names = FALSE)
  meshCoordinates[vert2print[[whichEE]],4] <- resultw[,3]
  write.table(meshCoordinates, paste(extractNames[iSNP],"_pvalues_WT.txt",sep=""), col.names = FALSE, row.names = FALSE)
  meshCoordinates[vert2print[[whichEE]],4] <- BHpvaluesw
  write.table(meshCoordinates, paste(extractNames[iSNP],"_BHpvalues_WT.txt",sep=""), col.names = FALSE, row.names = FALSE)
  
  meshCoordinates[vert2print[[whichEE]],4] <- results[,1]
  write.table(meshCoordinates, paste(extractNames[iSNP],"_beta_S2S.txt",sep=""), col.names = FALSE, row.names = FALSE)
  meshCoordinates[vert2print[[whichEE]],4] <- results[,3]
  write.table(meshCoordinates, paste(extractNames[iSNP],"_pvalues_S2S.txt",sep=""), col.names = FALSE, row.names = FALSE)
  meshCoordinates[vert2print[[whichEE]],4] <- BHpvaluess
  write.table(meshCoordinates, paste(extractNames[iSNP],"_BHpvalues_S2S.txt",sep=""), col.names = FALSE, row.names = FALSE)
  length(which(BHpvaluesw<0.05))
  length(which(BHpvaluess<0.05))
  
  {
    start.time <- Sys.time()
    signifw<-permFL_fast(Xn, Ywd, extract, A, NNmatrix, nPermutations, E=0.5, H=2, func_dir)
    signifs<-permFL_fast(Xn, Ysd, extract, A, NNmatrix, nPermutations, E=0.5, H=2, func_dir)
    end.time <- Sys.time()
    time.taken<-end.time - start.time
    time.taken
  }
  
  signw<-signifw$pvalues # get p-values
  signs<-signifs$pvalues # get p-values
  
  pfdr5TSBHw <- mt.rawp2adjp(signw[,1], proc=c("BH"), na.rm = FALSE)
  pvalueADJ5tsbh <- array(dim = length(signw[,1]))
  BHpvaluesTFCEw <- pfdr5TSBHw$adjp[order(pfdr5TSBHw$index),][,2]
  pfdr5TSBHs <- mt.rawp2adjp(signs[,1], proc=c("BH"), na.rm = FALSE)
  pvalueADJ5tsbh <- array(dim = length(signs[,1]))
  BHpvaluesTFCEs <- pfdr5TSBHs$adjp[order(pfdr5TSBHs$index),][,2]
  # 
  # #PRINT OUTPUT
  meshCoordinates[vert2print[[whichEE]],4] <- signw[,1]
  write.table(meshCoordinates, paste(extractNames[iSNP],"_pvaluesTFCE_WT.txt",sep=""), col.names = FALSE, row.names = FALSE)
  meshCoordinates[vert2print[[whichEE]],4] <- BHpvaluesTFCEw
  write.table(meshCoordinates, paste(extractNames[iSNP],"_BHpvaluesTFCE_WT.txt",sep=""), col.names = FALSE, row.names = FALSE)
  meshCoordinates[vert2print[[whichEE]],4] <- signs[,1]
  write.table(meshCoordinates, paste(extractNames[iSNP],"_pvaluesTFCE_S2S.txt",sep=""), col.names = FALSE, row.names = FALSE)
  meshCoordinates[vert2print[[whichEE]],4] <- BHpvaluesTFCEs
  write.table(meshCoordinates, paste(extractNames[iSNP],"_BHpvaluesTFCE_S2S.txt",sep=""), col.names = FALSE, row.names = FALSE)
  
  length(which(BHpvaluesTFCEw<0.05))
  length(which(BHpvaluesTFCEs<0.05))
}

# END
