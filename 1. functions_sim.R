# Codes were taken from https://github.com/klgrantham/compare-SC-SW

# Generate basic staircase design matrix
SCdesmat <- function(S, K=1, pre=1, post=1) {
  # Inputs:
  #  S - number of unique treatment sequences
  #  K - number of times each sequence is repeated/number of clusters per sequence
  #  pre - number of pre-switch measurement periods;default set to 1
  #  post - number of post-switch measurement periods;default set to 1
  # Output:
  #  Design matrix
  Xsc <- matrix(data=NA, nrow=S, ncol=(S+pre+post-1))
  for(i in 1:S) {
    Xsc[i,i:(i+pre-1)] <- 0
    Xsc[i,(i+pre):(i+pre+post-1)] <- 1
  }
  XscK <- Xsc[sort(rep(1:S, K)), ]
  return(XscK)
}

# Calculate staircase treatment effect variance

VarSClin <- function(S,K,m, rho0, r){
  a <- (1 + (m-1)*rho0)/m
  b <- r*rho0
  
  vartheta <- (2*((S^2+2)*a - (S^2-4)*b))/(K*S*(S^2-1))
  return(vartheta)
}

VarSCcat <- function(S,K,m, rho0, r){
  a <- (1 + (m-1)*rho0)/m
  b <- r*rho0
  
  fracterm <- ((a + sqrt(a^2-b^2))^S - b^S)/((a + sqrt(a^2-b^2))^S + b^S)
  vartheta <- (2*(a-b)^2)/(K*(S*(a-b)-sqrt(a^2-b^2)*fracterm))
  return(vartheta)
}

pow <- function(vars, effsize, siglevel=0.05){
  z <- qnorm(siglevel/2)
  pow <- pnorm(z + sqrt(1/vars)*effsize)
  return(pow)
}

