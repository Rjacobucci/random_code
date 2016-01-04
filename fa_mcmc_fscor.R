library(lavaan)
hs = HolzingerSwineford1939[7:12]
p = ncol(hs)
N = nrow(hs)
SampCov = cov(hs)


start <- c(rep(0.5,6),rep(0.5,6),rep(0,301))
pars = start


likelihood <- function(pars){
  lambda = pars[1:6]
  psi <- pars[7:12]
  int <- colMeans(hs) #pars[13:18]
  fscor = pars[13:313]
  
  c <- (2 * pi)^-(p/2) * (det(diag(psi)))^-0.5
  
  p.lik = rep(NA,N)
  for(i in 1:N){
    p.lik[i] = exp(-0.5*
          as.numeric(hs[i,] - int -  lambda * fscor[i]) %*% 
          solve(diag(psi)) %*% 
          t(hs[i,] - int -  lambda * fscor[i])
        )
  }
 (c *  sum(p.lik))
  
}




prior <- function(pars){
  lam.prior1 = pars[1:6]
  psi.prior1 = pars[7:12]
  fscor.prior1 = pars[13:313]
  
  lam.prior2 = rep(NA,6)
  psi.prior2 = rep(NA,6)
  fscor.prior2 = rep(NA,N)
  
  for(i in 1:length(lam.prior1)){
    lam.prior2[i] = dnorm(lam.prior1[i],mean=0.8, sd=1, log = T)
  }
  
  for(i in 1:length(psi.prior1)){
    psi.prior2[i] = dnorm(psi.prior1[i],mean=0.1, sd=1, log = T)
  }
  
  for(i in 1:length(fscor.prior1)){
    fscor.prior2[i] = dnorm(fscor.prior1[i],mean=0, sd=1, log = T)
  }

  return(sum(lam.prior2)+sum(psi.prior2)+sum(fscor.prior2))
}


posterior <- function(pars){
  return (likelihood(pars) + prior(pars))
}




proposalfunction <- function(pars){
  return(rnorm(length(start),mean = pars, sd= c(rep(0.1,6),rep(0.1,6),rep(0.1,N))))
}


run_metropolis_MCMC <- function(start, iterations){
  chain = array(dim = c(iterations+1,length(start)))
  chain[1,] = start
  for (i in 1:iterations){
    proposal = proposalfunction(chain[i,])
    
    probab = exp(posterior(proposal) - posterior(chain[i,]))
    if (runif(1) < probab){
      chain[i+1,] = proposal
    }else{
      chain[i+1,] = chain[i,]
    }
  }
  return(chain)
}


chain = run_metropolis_MCMC(start, 100)

burnIn = 5000
acceptance = 1-mean(duplicated(chain[-(1:burnIn),]))


colMeans(chain[-c(1:5000),1:6])


mod <- '
f =~ NA*x1 + x2 + x3 + x4 + x5 + x6
f~~1*f
'
out = cfa(mod,hs)
summary(out)
