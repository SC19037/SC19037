#' @title  Verify the Law of Large Numbers by figure
#' @description  Verify the law of large Numbers by figure
#' @param r the distribution
#' @param mean the mean of this distribution
#' @param n the number of samples
#' @return the figure
#' @examples
#' \dontrun{
#' n <- 1000
#' LLN(rbinom(n, 1, 0.3), 0.3, n)
#' LLN(runif(n,0,1),0.5,n)
#' LLN(rpois(n, 1), 1, n)
#' LLN(rexp(n, 1), 1, n)
#' }
#' @export
LLN <- function(r, mean, n){
  y <- numeric(n)
  for (i in 1:n){
    y[i] <- mean(sample(r,i,replace=TRUE))
  }
  lln_data <- data.frame(1:n, y)
  plot(lln_data, xlab = "Sample Size", ylab = "Average Value", col=" light green")
  abline(h = mean, col = "red")
}


#' @title  Verify the Central Limit Theorem by figure
#' @description  Verify the Central Limit Theorem by figure
#' @param r the distribution
#' @param parameter the parameter of this distribution
#' @param mean the mean of this distribution
#' @param var the variance of this distribution
#' @param n the number of samples
#' @param N the repetition times
#' @return the figure
#' @examples
#' \dontrun{
#' n <- c(1,10,30,100)
#' par(mfrow=c(4,4))
#' CLT(rbinom, parameter = c(1, 0.3), 0.3, 0.21, n)
#' CLT(rpois, parameter = 1, 1, 1, n)
#' CLT(runif, parameter = c(0,1), 0.5, 1/12, n)
#' CLT(rexp, parameter = 1, 1, 1, n)
#' }
#' @export 
CLT <- function (r, parameter, mean, var, n, N=1000){
  
  if (length(parameter) == 1){
    for (i in n) {
      x <- matrix(r(i*N, parameter), nrow = N, ncol = i)
      x <- (apply(x, 1, sum) - i * mean )/(sqrt(i) * sqrt(var))
      hist(x, col = 'light yellow', probability = T, main = paste("n = ",i),
           ylim = c(0, max(0.5, density(x)$y)))
      lines(density(x), col='red', lwd = 3)
      curve(dnorm(x), col='blue', lwd = 3, add = T)
    }
  }
  
  if (length(parameter) == 2){
    p1 <- parameter[1]
    p2 <- parameter[2]
    for (i in n) {
      x <- matrix(r(i*N, p1, p2 ), nrow = N, ncol = i) 
      x <- (apply(x, 1, sum) - i * mean )/(sqrt(i) * sqrt(var))
      hist(x, col = 'light yellow', probability = T, main = paste("n = ",i),
           ylim = c(0, max(0.5, density(x)$y)))
      lines(density(x), col='red', lwd = 3)
      curve(dnorm(x), col='blue', lwd = 3, add = T)
    }
  }
}

