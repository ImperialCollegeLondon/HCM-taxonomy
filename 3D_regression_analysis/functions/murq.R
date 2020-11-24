#' Mass Univariate Regression.
#'
#' Fit a linear regression model at each vertex of a 3D atlas. Inputs of the function are a NxC matrix X modeling an effect under study (N = number of subjects, C = number of variables + intercept term),
#' a NxV imaging matrix Y containing the values of a 3D phenotype at each atlas vertex (V = number of vertices in the 3D mesh), and an array (extract) containg the positions of
#' variables in X of which extract the informations of interest. The output is a Vx(3xlength(extract)) matrix containing the regression coefficient, its related t-statistic and the p-value at each vertex of the computational model for each variable specifiec in extract.
#' @param X is the design matrix. Number of rows = number of subjects in the study, number of columns = number of vertices in the atlas. Numerical varable must be normalized to 0-mean and unit-standard deviation. Categorical variables must be coded using dummy coding. The first column should contain the intercept (all 1s).
#' @param Y is the imaging matrix. Number of rows = N. Number of columns = V.
#' @param extract is an array expressing which covariates in X you want to extract.
#' @keywords mur regression
#' @export
#' @examples
#' extract <- c(1,3) #extract the first and third covariate.
#' result <- murq(X, Y, extract)
#' betas <- result[,1]
#' tstatistics <- result[,2]
#' pvalues <- result[,3]

murq <- function(X, Y, extract){

  nPoints <- ncol(Y)
  #number of points in the mesh/vertexes under study

  tMUR <- matrix(0, ncol=3*length(extract), nrow=nPoints)
  #beta and p obtained at each point

  dofC <- nrow(X) - ncol(X)
  
  xTxInv <- solve(t(X) %*% X)
  #(tX X)^-1
  
  #QR decomposition of X matrix into an orthogonal matrix and a triangular matrix
  QR <- qr(X)
  for(y in 1:nPoints){
    #do the regression for all the vertexes of the atlas
    #the for cycle is less time consuming than a foreach in this cas
    beta <- solve.qr(QR, Y[,y])
    #regression coefficients
    resid <-  Y[,y] - (X %*% beta)
    #residuals

    for(iEx in 1:length(extract)){
      se <- (t(resid) %*% resid)/dofC
      se <-  sqrt(as.numeric(se) * xTxInv[extract,extract])

      t <- beta[extract[iEx]]/se

      tMUR[y,1+(iEx-1)*3] <- beta[extract[iEx]] #beta
      tMUR[y,2+(iEx-1)*3] <- t #t
      tMUR[y,3+(iEx-1)*3] <- 2 * pt(-abs(t), df=nrow(Y)-ncol(X)) # p-value
    }

  }

  return(tMUR)

}
