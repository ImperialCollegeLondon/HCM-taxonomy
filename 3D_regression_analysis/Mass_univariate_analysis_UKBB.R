###################################################################################
################# Mass Univariate Analysis for UKBB data ##########################
###################################################################################
###################################################################################
rm(list = ls(all = TRUE))  
#start with an empty environment

install.packages("data.table")
install.packages("BiocManager")
BiocManager::install("multtest")
install.packages("devtools")
install_github("UK-Digital-Heart-Project/mutools3D", build_vignettes = TRUE)
install.packages("Rcpp")
install.packages("RcppArmadillo")
install.packages("RcppEigen")
install.packages("float")
install.packages("doParallel")
install.packages("foreach")
install.packages("igraph")
install.packages("plyr")
# install all packages required

setwd("~/cardiac/Experiments_of_Maria/HCM_expressivity/3D_regression_analysis/")
library(data.table)
library(multtest)
library(mutools3D)
library(Rcpp) # required in permFL_fast
library(RcppArmadillo) # required in mulitply.cpp
library(RcppEigen) # required in mulitply.cpp
library(float) # required in permFL_fast

func_dir<-"~/cardiac/Experiments_of_Maria/HCM_expressivity/3D_regression_analysis/"
#folder directory of the code. Will be used for the functions as: source(paste(dir,"functions/xxx", sep="/")) 
#in the permFL_fast function.

sourceCpp(paste(func_dir,"functions/multiply.cpp",sep=""))
#include a cpp functions for faster multiplication of matrices

# source(paste(func_dir,"functions/murq.R",sep=""))
# include the mur function with QR decomposition (slightly faster!), HC4m = FALSE
#or include the mur function from mutools3D package

#include the functions for TFCE from mutools3D package

source(paste(func_dir,"functions/permFL_fast.R",sep=""))
#include the permFL_fast functions with murq and HC4m = FALSE

#CLINICAL DATA MATRIX
#NCOL = N COVARIATES UNDER STUDY
#NROW = N PATIENTS

inputClinical <- readRDS(paste(func_dir,"data/X.rds",sep=""))

# HCM summary table
summary_table<-as.data.frame(read.table(paste(func_dir,"data/summarytable_meta_UKB_HCM.txt", sep=""), header=T))
head(summary_table)

inputClinical$Age <- scale(inputClinical$Age)
inputClinical$BSA <- scale(inputClinical$BSA)
inputClinical$SBP <- scale(inputClinical$SBP)
#inputClinical$DBP <- scale(inputClinical$DBP)

X <- data.matrix(inputClinical[,-c(1)])
X <- cbind(1, X)
X<-X[,c(1:6)]
dim(X)
head(X)

## HCM steps
setwd("output")
l<-c("Sarc-ve", "Sarc","Sarc+ve")
# create folders to save the output for each step. No need to do this at each step! 
for (i in 1:length(l)){
  dir.create(l[i])}

#step 0 vs Sarc-ve
{
  p0<-which(summary_table$HCMstep0==1)
  p1<-which(summary_table$HCMstep1==1)
  
  pcomb<-c(p0,p1)
  pcomb<-unique(pcomb)
  summary_new1<-summary_table[pcomb,]
  ps<-match(summary_new1$eid,inputClinical$eid)
  
  Xn<-as.data.frame(cbind(X[ps,], summary_new1$HCMstep1))
  Xn<-as.matrix(Xn)
  colnames(Xn)[7]<-c("HCM")
  length(which(Xn[,7]==1))
  setwd("output/step1")
}
#step 0 vs Sarc
{
  p0<-which(summary_table$HCMstep0==1)
  p2b<-which(summary_table$HCMstep2B==1)
  p3a<-which(summary_table$HCMstep3A==1)
  p3b<-which(summary_table$HCMstep3B==1)
  p3ab<-c(p3a,p3b)
  p2b <-p2b[!(p2b%in% p3ab)]
  
  pcomb<-c(p0,p2b)
  summary_new2b<-summary_table[pcomb,]
  ps<-match(summary_new2b$eid,inputClinical$eid)
  ps<-na.omit(ps)
  
  Xn<-as.data.frame(cbind(X[ps,], summary_new2b$HCMstep2B))
  Xn<-as.matrix(Xn)
  colnames(Xn)[7]<-c("HCM")
  length(which(Xn[,7]==1))
  setwd(paste(func.dir,"output/HCM/step2b_only", sep=""))
}
#step 0 vs Sarc+ve
{
  p0<-which(summary_table$HCMstep0==1)
  p3a<-which(summary_table$HCMstep3A==1)
  p3b<-which(summary_table$HCMstep3B==1)
  
  p3ab<-c(p3a,p3b)
  p3ab<-unique(p3ab)
  p3ab<-p3ab[order(p3ab)]

  pcomb<-c(p0,p3ab)
  pcomb<-unique(pcomb)
  
  merge_col3ab<-as.integer(ifelse(summary_table$HCMstep3A==1|summary_table$HCMstep3B==1,1,0))
  summary_table$HCMstep3A_3B<-merge_col3ab
  summary_new3ab<-summary_table[pcomb,]
  length(which(summary_new3ab$HCMstep3A_3B==1))
  ps<-match(summary_new3ab$eid,inputClinical$eid)
  ps<-na.omit(ps)
  
  Xn<-as.data.frame(cbind(X[ps,], summary_new3ab$HCMstep3A_3B))
  Xn<-as.matrix(Xn)
  length(which(Xn[,7]==1))
  colnames(Xn)[7]<-c("HCM")
  setwd("output/step3ab")
}

#IMAGING DATA MATRIX
#NCOL = N POINTS ON THE ATLAS
#NROW = N PATIENTS
Y <- readRDS(paste(func_dir,"data/Y.rds", sep=""))
#DATA PRE-PROCESSING
Y<-Y[ps,] # make Y the same length as X and in the same order.
dim(Y)

#NUMBER OF CORES TO USE
# All core detected here
# nofCores = 48

# 1000 permutations
nPermutations=1000

whichEE <- 3
#1 endo, 2 epi, 3 full shape
endoEpi <- read.table(paste(func_dir,"data/endo_epi.txt", sep=""))
vert2print <- list(which(endoEpi[,4]==0),which(endoEpi[,4]==1),1:length(endoEpi[,4]))
Y<-Y[,vert2print[[whichEE]]]

##READ THE NNLIST ASSOCIATED TO EACH VERTEX
NNmatrix <- readRDS(paste(func_dir,"data/redNNmatrix.rds", sep=""))
##READ AREAS ASSOCIATED TO EACH VERTEX
A <- readRDS(paste(func_dir,"data/LVarea.rds", sep=""))

# use the meshCoordinates from the mean of 13k I have created
mesh_Coordinates <- read.table(paste(func.dir,"data/meshCoordinates.txt", sep=""), quote="\"", comment.char="", header=T)
colnames(mesh_Coordinates) <- c("x", "y", "z")

extractNames<-colnames(Xn)
length(extractNames)
extractNames
extract=7

# Run the 3D Mass univariate regression analysis
{
  #extract only the pheno betas
  result <- mur(Xn,Y, extract)
  
  #MULTIPLE TESTING CORRECTION
  corrected <- mt.rawp2adjp(result[,3], proc=c("BH"), na.rm = FALSE)
  pvalueADJ5tsbh <- array(dim = length(result[,3]))
  BHpvalues <- corrected$adjp[order(corrected$index),][,2]

  meshCoordinates <- cbind(mesh_Coordinates,99999)
  #PRINT OUTPUT
  meshCoordinates[vert2print[[whichEE]],4] <- result[,1]
  write.table(meshCoordinates, paste(extractNames[extract],"_beta_WT.txt",sep=""), col.names = FALSE, row.names = FALSE)
  meshCoordinates[vert2print[[whichEE]],4] <- result[,3]
  write.table(meshCoordinates, paste(extractNames[extract],"_pvalues_WT.txt",sep=""), col.names = FALSE, row.names = FALSE)
  meshCoordinates[vert2print[[whichEE]],4] <- BHpvalues
  write.table(meshCoordinates, paste(extractNames[extract],"_BHpvalues_WT.txt",sep=""), col.names = FALSE, row.names = FALSE)
  
  signif<-permFL_fast(Xn, Y, extract, A, NNmatrix, nPermutations, E=0.5, H=2, func_dir)
   
  sign<-signif$pvalues # get p-values
  
  pfdr5TSBH <- mt.rawp2adjp(sign[,1], proc=c("BH"), na.rm = FALSE)
  pvalueADJ5tsbh <- array(dim = length(sign[,1]))
  BHpvaluesTFCE <- pfdr5TSBH$adjp[order(pfdr5TSBH$index),][,2]
  # 
  # #PRINT OUTPUT
  meshCoordinates[vert2print[[whichEE]],4] <- sign[,1]
  write.table(meshCoordinates, paste(extractNames[extract],"_pvaluesTFCE_WT.txt",sep=""), col.names = FALSE, row.names = FALSE)
  meshCoordinates[vert2print[[whichEE]],4] <- BHpvaluesTFCE
  write.table(meshCoordinates, paste(extractNames[extract],"_BHpvaluesTFCE_WT.txt",sep=""), col.names = FALSE, row.names = FALSE)
}

# END
