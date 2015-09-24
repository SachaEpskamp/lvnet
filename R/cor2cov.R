# H. Seltman, May 2014

# Goal: convert a correlation matrix and variance vector
#       into the corresponding covariance matrix
#
# Input: 
#   'corMat' is a square matrix with 1's on the diagonal
#      and valid correlations on the off-diagonal
#   'varVec' is a valid variance vector, with length
#      matching the dimension of 'covMat'.  A single
#      row or single column matrix is also allowed.
# Output:
#   the covariance matrix
# 
# A warning is given if the covariance matrix is not
#   positive definite.
#
cor2cov = function(corMat, varVec) {
  # test the input
  if (!is.matrix(corMat)) stop("'corMat must be a matrix")
  n = nrow(corMat)
  if (ncol(corMat) != n) stop("'corMat' must be square")
  if (mode(corMat) != "numeric") stop("'corMat must be numeric")
  if (mode(varVec) != "numeric") stop("'varVec must be numeric")
  if (!is.null(dim(varVec))) {
    if (length(dim(varVec)) != 2) stop("'varVec' should be a vector")
    if (any(dim(varVec)==1)) stop("'varVec' cannot be a matrix")
    varVec = as.numeric(varVec) # convert row or col matrix to a vector
  }
  if (!all(diag(corMat) == 1)) stop("correlation matrices have 1 on the diagonal")
  if (any(corMat < -1 | corMat > +1)) 
    stop("correlations must be between -1 and 1")
  if (any(varVec <= 0)) stop("variances must be non-negative")
  if (length(varVec) != n) stop("length of 'varMat' does not match 'corMat' size")
  
  # Compute the covariance
  sdMat = diag(sqrt(varVec))
  rtn = sdMat %*% corMat %*% t(sdMat)
  if (det(rtn)<=0) warning("covariance matrix is not positive definite")
  return(rtn)
}
