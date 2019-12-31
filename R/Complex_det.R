#' @title The determinant of a complex matrix
#' @description The determinant of a complex matrix
#' @param x the complex matrix
#' @return the determinant of x\code{n}
#' @examples
#' \dontrun{
#' x <- matrix(c(1+1i,2,3+4i,2i,3+7i,1-3i,2+4i,-1i,9),3,3)
#' Complex_det(x)
#' }
#' @export
Complex_det <- function(x){
  if (length(x) == 1){
    return (x)
  }
  c <- numeric(nrow(x))
  det <- 0
  n <- nrow(x)
  for (i in 1:n){
    y <- x[-i,-1]
    c[i] <- (-1)^(i+1) * (x[i,1] %*% Complex_det(y))
    det <- det + c[i]
  }
  return(det)
}