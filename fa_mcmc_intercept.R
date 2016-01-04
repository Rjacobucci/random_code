library(lavaan)
library(truncnorm)
hs = data.frame(scale(HolzingerSwineford1939[7:12]))
p = ncol(hs)
N = nrow(hs)
SampCov = cov(hs) + colMeans(hs) %*% t(colMeans(hs))


start <- c(rep(0.5,6),rep(0.5,6),rep(0,6))
pars = start


likelihood <- function(pars){
  lambda = pars[1:6]
  psi <- pars[7:12]
  int <- pars[13:18] #pars[13:18]
  
  c <- (2 * pi)^-(p/2) * (det(diag(psi)))^-0.5
  

    p.lik = exp(-0.5*
                     as.numeric(colMeans(hs) - int -  lambda) %*% 
                     solve(diag(psi)) %*% 
                     (colMeans(hs) - int -  lambda)
    )

  (c *  p.lik)
}



prior <- function(pars){
  lam.prior1 = pars[1:6]
  psi.prior1 = pars[7:12]
  int.prior1 = pars[13:18]
  lam.prior2 = rep(NA,6)
  psi.prior2 = rep(NA,6)
  int.prior2 = rep(NA,6)
  
  for(i in 1:length(lam.prior1)){
    lam.prior2[i] = dnorm(lam.prior1[i],mean=0.5, sd=1, log = T)
  }
  
  for(i in 1:length(psi.prior1)){
    psi.prior2[i] = dnorm(psi.prior1[i],mean=0.5, sd=1, log = T)
  }
  
  for(i in 1:length(int.prior1)){
    int.prior2[i] = dnorm(int.prior1[i],mean=0, sd=1, log = T)
  }

  return(sum(lam.prior2)+sum(psi.prior2) + sum(int.prior2))
}


posterior <- function(pars){
  return (likelihood(pars) + prior(pars))
}




proposalfunction <- function(pars){
  return(rtruncnorm(18,a=0,mean = pars, sd= c(rep(0.1,6),rep(0.1,6),rep(0.1,6))))
}


run_metropolis_MCMC <- function(start, iterations){
  chain = array(dim = c(iterations+1,18))
  chain[1,] = start
  for (i in 1:iterations){
    proposal = proposalfunction(chain[i,])
    
  
    probab = exp(posterior(proposal) - posterior(chain[i,]))
    
    out = try(runif(1) < probab,silent=T)
    if(inherits(out, "try-error") | is.na(out)==TRUE) { 
        chain[i+1,] = chain[i,]
    }else{
      if (runif(1) < probab){
        chain[i+1,] = proposal
      }else{
        chain[i+1,] = chain[i,]
      }
    }

  }
  return(chain)
}


chain = run_metropolis_MCMC(start, 10000)

burnIn = 5000
acceptance = 1-mean(duplicated(chain[-(1:burnIn),]))


colMeans(chain[-c(1:5000),1:18])


mod <- '
f =~ NA*x1 + x2 + x3 + x4 + x5 + x6
f~~1*f
'
out = cfa(mod,hs,meanstructure = T)
summary(out)
