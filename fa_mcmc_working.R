library(lavaan)
library(truncnorm)
hs = HolzingerSwineford1939[7:12]
p = ncol(hs)
N = nrow(hs)
SampCov = cov(hs)


start <- c(rep(0.5,6),rep(0.5,6))
#pars = start


likelihood <- function(pars){
  lambda = pars[1:6]
  psi <- pars[7:12]
  #int <- pars[13:18] #pars[13:18]
  
  ImpCov <- lambda %*% t(lambda) + diag(psi)
  
  fit_ret <- function(ImpCov,SampCov,p){
    fit=  0.5*(log(det(ImpCov)) + sum(diag(SampCov %*% solve(ImpCov))) - log(det(SampCov))  - p)
    fit
  }
  
  fit = fit_ret(ImpCov,SampCov,p)
  
  log.lik <- function(N,p,SampCov,fit){
    c <- N*p/2 * log(2 * pi)
    logl_sat = -c -(N/2) * log(det(SampCov)) - (N/2)*p
    logl = -N * (fit- logl_sat/N)
    logl
  }
  
  log.lik(N,p,SampCov,fit)
}



prior <- function(pars){
  lam.prior1 = pars[1:6]
  psi.prior1 = pars[7:12]
  lam.prior2 = rep(NA,6)
  psi.prior2 = rep(NA,6)
  
  for(i in 1:length(lam.prior1)){
    lam.prior2[i] = dnorm(lam.prior1[i],mean=0.8, sd=1, log = T)
  }
  
  for(i in 1:length(psi.prior1)){
    psi.prior2[i] = dnorm(psi.prior1[i],mean=0.1, sd=1, log = T)
  }
  

  return(sum(lam.prior2)+sum(psi.prior2))
}


posterior <- function(pars){
  return (likelihood(pars) + prior(pars))
}




proposalfunction <- function(pars){
  return(rtruncnorm(12,a=0,mean = pars, sd= c(rep(0.1,6),rep(0.1,6))))
}


run_metropolis_MCMC <- function(start, iterations){
  chain = array(dim = c(iterations+1,12))
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


colMeans(chain[-c(1:5000),])


mod <- '
f =~ NA*x1 + x2 + x3 + x4 + x5 + x6
f~~1*f
'
out = cfa(mod,hs)
summary(out)
