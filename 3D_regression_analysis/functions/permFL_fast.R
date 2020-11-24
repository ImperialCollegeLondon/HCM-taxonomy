permFL_fast<-function(X, Y, extract, A, NNmatrix, nPermutations, E=0.5, H=2, func_dir){
  library(doParallel)
  library(parallel)
  library(foreach)
  library(mutools3D)
  registerDoParallel(detectCores())
  #or detect the number of cores available and use them
  library(igraph)
  require(multtest)
  library(Rcpp)
  library(plyr)
  library(float)
  
  sourceCpp(paste(func_dir,"functions/multiply.cpp",sep=""))
  # include a cpp functions for faster multiplication of matrices!
  
  # source(paste(func_dir,"functions/murq.R",sep=""))
  #include the mur function with QR decomposition (slightly faster!), HC4m = FALSE
  #source(paste(func_dir,"functions/mur.R",sep=""))
  #or include the mur function, HC4m = FALSE
  
  #source(paste(func_dir,"functions/TFCE.R",sep=""))
  #include the functions for TFCE
  
  set.seed(1234)
  #set seed for reproducibility
  
  Z <- X[,-extract]
  #compute Z (nuisance matrix)
  
  Zpinv<-Z %*%solve(t(Z) %*% Z) %*% t(Z)
  Ypr<-eigenMapMatMult(Zpinv,Y) # faster multiplication!
  Yp<-Y-Ypr
  
  Rz<-fl(Yp) # float for precision matrix
  
  #parallelization for iF=1:10 in a loop
  nPerm<-nPermutations/10
  permresP<-vector(mode = "list", length=nPerm) # a list of nPerm length
  for (iR in 1:nPerm){
    resP <- foreach(iF=1:10, .combine=rbind)%dopar%{
      Yper<- Rz@Data[sample(1:nrow(Rz@Data)),]
      #Y permuted for the Freedman and Lane procedure
      
      resMUR <- mur(X,Yper,extract)
      rm(Yper)
      computed <- matrix(0, ncol=ncol(Y), length(extract))
      
      computed <- TFCE(h=round(resMUR[,2],2), A=A, NNmatrix=NNmatrix, E=E, H=H)
      #compute TFCE
      return(computed)
    }
    permresP[[iR]]<-resP
    
  }
  
  resP<-ldply(permresP,rbind)
  #for each element of a list, apply function then combine results into the resP data frame
  closeAllConnections()
  
  significance <- matrix(0, ncol=length(extract)*2, nrow=ncol(Y))
  #TFCE derived p-values 
 
  results <- matrix(0, ncol=3*length(extract), nrow=ncol(Y))
  results <- mur(X,Y,extract)
  tfceScores <- list()
  #compute the residual matrix of Z
  iEx<-1
  tfceScores[[iEx]] <- TFCE(results[,2+(iEx-1)*3], A=A, NNmatrix=NNmatrix, E=E, H=H)
  #list of TFCE scores to analyse
  TFCEmatrix <- resP[seq(1,nrow(resP), by=length(extract)),]
  
  
  minimum = sort(apply(TFCEmatrix,1,min))
  if (length(which(minimum<0)>0)) { 
    thrMin = minimum[ceiling(0.05*nrow(TFCEmatrix))]
  } else {
    thrMin = 0
  } 
  
  maximum = sort(apply(TFCEmatrix,1,max))
  if (length(which(maximum>0)>0)) {
    thrMax = maximum[floor(0.95*nrow(TFCEmatrix))]
  }else{
    thrMax = 0
  } 
  
  for(a in 1:ncol(Y)){
    if(tfceScores[[iEx]][a]>=0){
      significance[a,1+(iEx-1)*2]  <- length(which(TFCEmatrix[,a]> tfceScores[[iEx]][a]))/nPermutations
      if(tfceScores[[iEx]][a] > thrMax) significance[a,2+(iEx-1)*2] = 1
    }
    
    if(tfceScores[[iEx]][a]<0){
      significance[a,1+(iEx-1)*2]  <- length(which(TFCEmatrix[,a]< tfceScores[[iEx]][a]))/nPermutations
      if(tfceScores[[iEx]][a] < thrMin) significance[a,2+(iEx-1)*2] = 1
    }
    
    if(significance[a,1+(iEx-1)*2]==0) significance[a,1+(iEx-1)*2] <- 1/nPermutations #minimum pvalue achievable.
  }
  
  
  
  length(which(significance[,1]<0.05))
  TFCEresults = list("pvalues" = significance, "TFCEmatrix" = TFCEmatrix, "tfceScores" = tfceScores)
  
  #else TFCEresults = significance
  
  rm(significance)
  rm(TFCEmatrix)
  rm(tfceScores)
  
  return(TFCEresults)    
}
