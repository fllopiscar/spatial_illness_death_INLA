'inla.rgeneric.multiLeroux3d.model' <- function(
  cmd = c("graph", "Q", "mu", "initial", "log.norm.const",
          "log.prior", "quit"),
  theta = NULL) {
  
  #Internal function
  interpret.theta <- function() {
    
    lambda = 1 / (1 + exp(-theta[1L]))
    prec.var = exp(theta[2L:4L]) 
    rho=2*exp(theta[5L:7L])/(1+exp(theta[5L:7L]))-1
    
    prec.matrix.var        <- diag(1/prec.var)
    prec.matrix.var[2,1] <- rho[1]/sqrt(prec.var[1]*prec.var[2])
    prec.matrix.var[3,1] <- rho[2]/sqrt(prec.var[1]*prec.var[3])
    prec.matrix.var[1,2:3] <- prec.matrix.var[2:3,1]
    prec.matrix.var[2,3]   <- rho[3]/sqrt(prec.var[2]*prec.var[3])
    prec.matrix.var[3,2]   <- prec.matrix.var[2,3]
    prec.matrix.var <- solve(prec.matrix.var)
    
    return(list(lambda=lambda, prec.var=prec.var,
                rho=rho,prec.matrix.var=prec.matrix.var))
  }
  
  graph <- function(){
    require(Matrix)
    
    prec.matrix.var <-matrix(nrow=3,ncol=3,1)
    prec.matrix.sp <- Diagonal(nrow(W), x = 1) + W
    
    return(kronecker(prec.matrix.var,prec.matrix.sp))
  }
  
  Q <- function() {
    require(Matrix)
    
    param <- interpret.theta()
    D <- Diagonal(nrow(W), apply(W,1,sum))
    
    prec.matrix.sp <- (1-param$lambda)*Diagonal(nrow(W), x = 1) + param$lambda*(D-W)
    
    return( kronecker(param$prec.matrix.var, prec.matrix.sp) )
  }
  
  mu <- function()
  {
    return(numeric(0))
  }
  
  log.norm.const <- function() {
    return(numeric(0))
    
  }
  
  log.prior <- function() {
    param = interpret.theta()
    
    
    res <-log(1) +log(param$lambda) + log(1 - param$lambda)+
          log(MCMCpack::dwish(W=param$prec.matrix.var,v,S=diag(rep(d.elem,3))))+
          sum(log(param$prec.var))+
          sum(log(2) + theta[5:7]-2*log(1 + exp(theta[5:7])))
      
      
    return(res)
  }
  
  initial <- function() {
    return(c(0,0,0,0,0,0,0))
  }
  
  quit <- function() {
    return(invisible())
  }
  
  res <- do.call(match.arg(cmd), args = list())
  return(res)
}
