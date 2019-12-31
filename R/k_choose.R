#' @title  Selecting the optimal ridge regression parameters by GCV
#' @description  Selecting the optimal ridge regression parameters by GCV
#' @param X the independent variable
#' @param Y the response variable
#' @param kmax the maximum value of the ridge regression parameter
#' @param N the number of cross validation
#' @return the optimal ridge regression parameter
#' @examples
#' \dontrun{
#' library(MASS)
#' data <- read.table("E:/Rcpp/data.txt", header = TRUE)
#' x <- as.matrix(data[,2:21])
#' y <- data[,1]
#' k_best1 <- k_choose1(x, y, kmax = 2, N = 2000)
#' print(k_best1)
#'  #compare with the optimal parameter gotten by existing ridge package 
#' select(lm.ridge(y~x, data = data, lambda = seq(0,2,0.001)))}
#' @export
k_choose1 <- function(X, Y, kmax, N){
  K <- seq(0, kmax, kmax/N)
  GCV <- numeric(0)
  for (k in K) {
    # get the reponse variable with parameter k
    Y_k <-  X %*% solve(t(X) %*% X + k * diag(ncol(X))) %*% t(X) %*% Y
    # compute the value of generalized cross validation with parameter k
    GCV_k <- sum((Y - Y_k)^2) / (nrow(X) - sum(matrix((svd(X)$d)^2/((svd(X)$d)^2 
              + rep(k, rep(length(svd(X)$d), 1))), length(svd(X)$d))))^2
    GCV <- c(GCV, GCV_k)
  }
  # return the optimal parameter k which makes the GCV minimum
  return (K[which.min(GCV)])
}

#' @title  Selecting the optimal ridge regression parameters by minimuming GIC
#' @description Getting the optimal ridge regression parameters by minimuming GIC
#' @param X the independent variable
#' @param Y the response variable
#' @return the optimal ridge regression parameter
#' @examples
#' \dontrun{
#' data <- read.table("E:/Rcpp/data.txt", header = TRUE)
#' x <- as.matrix(data[,2:21])
#' y <- data[,1]
#' k_best2 <- k_choose2(x, y)
#' @export
k_choose2 <- function(X, Y){ 
  beta <- solve(t(X) %*% X) %*% t(X) %*% Y
  ng <- length(Y)
  numzero <- length(which(beta == 0))
  lenbeta <- length(beta)
  temp <- 0
  for(i in 1:ng)
  {
    temp <- temp + 1/ng * (Y[i] - sum(X[i,] * beta))^2   
  }
  GIC <- log(temp) + (lenbeta - numzero - 1) * log(log(ng))/ng * log(max(lenbeta, ng))
  return(GIC)
}