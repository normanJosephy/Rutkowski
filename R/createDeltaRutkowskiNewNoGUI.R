#' Computes small delta and portfolio G,H given a (u,d) pair
#'
#' Calls collectGHRutkowski() and rutkowskiDelta() and returns a list
#' consisting of delta, portfolio matrices GMatrix, HMatrix and
#' matrix pathMatrix, all for one (u,d) pair.
#'
#' @param u u parameter
#' @param d d parameter
#' @export createDeltaRutkowskiNew
#' @return outputList A list consisting of delta, GMatrix, HMatrix, pathMatrix
createDeltaRutkowskiNew = function(u,d) {
#  unpackList(myEnv)
#  unpackList(computedEnv)
  pathMatrix = computedEnv$paths
  lambda = myEnv$lambda
  mu     = myEnv$mu
  r      = myEnv$r
  K      = myEnv$K
  #
  nPaths = ncol(pathMatrix)
  nTimes = nrow(pathMatrix)
  rNames = paste('time-',0:(nTimes-1),sep='')
  cNames = paste('path-',1:nPaths,sep='')
  HMatrix = matrix(NA,nrow=nTimes,ncol=nPaths,dimnames=list(rNames,cNames))
  GMatrix = matrix(NA,nrow=nTimes,ncol=nPaths,dimnames=list(rNames,cNames))
  colnames(pathMatrix) = paste('path-',1:nPaths,sep='')
  rownames(pathMatrix) = paste('time-',0:(nTimes-1),sep='')
  delta = matrix(NA,nrow=nTimes-1,ncol=nPaths)
  rownames(delta) = paste('time-',1:(nTimes-1),sep='')
  colnames(delta) = cNames
  for (i in 1:nPaths) {
    path        = pathMatrix[,i]
    oList       = collectGHRutkowski(path,lambda,mu,r,u,d,K)
    G           = oList$hgInitial[,'G']
    H           = oList$hgInitial[,'H']
    GMatrix[,i] = G
    HMatrix[,i] = H
    delta[,i]   = rutkowskiDelta(path,H,G,lambda,mu,r)
  }
  # Prepare delta
  del = rbind(rep(NA,ncol(delta)),delta)
  rownames(del) = NULL # paste('path-',1:nPaths,sep='')
  colnames(del) = rep('deltaRut',ncol(delta)) # c('time-1','time-2')
  outputList = list(delta=delta,GMatrix=GMatrix,HMatrix=HMatrix,pathMatrix=pathMatrix)
  invisible(outputList)
}

