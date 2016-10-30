library(lavaan)
library(semTools)
library(GPArotation)
library(OpenMx)
library(ggplot2)
d = function(chisq,df,N) (chisq -df)/(N-1)
rmsea = function(ncp,df) sqrt(ncp/df)


hs <- HolzingerSwineford1939[,7:15]
hs.scale <- data.frame(scale(hs))

dat.train <- hs.scale[1:150,]
cov.train <- cov(dat.train)
dat.test <- hs.scale[151:301,]
cov.test <- cov(dat.test)

#vals = as.numeric(coef(efa.lav))

# max number of factors per # of items
# 6 factors for 9 items


ret <- list()
ret2 <- list()
floads.list <- list()
influ <- list()
conv <- list()
num.facs <- list()
SHRINK = 0
count = 0
counts=50
#res2 <- data.frame(matrix(NA,counts,3))
#coefs = rep(1,14)

nfac = 1
efa.lav = efaUnrotate(dat.train,nfac)
coefs = as.numeric(coef(efa.lav))


#own = matrix(NA,counts+1,(nfac*9 + 9))
starts = as.numeric(coefs[-c((9*nfac + 1):(nfac*9 + 9))])



while(count < counts){

  count = count + 1
  SHRINK <- 0.0007*(count-1) # 0.01 works well
  num.facs[[count]] <- nfac


  fitReg <- regsem(efa.lav,lambda=SHRINK,type="ridge",gradFun="none",
                   optMethod="optimx",hessFun="none",parallel="no",
                   Start=starts,subOpt="L-BFGS-B",fac.type="efa")#own[count,]

  conv[[count]] = fitReg$out$convcode

  floads <- matrix(coef(fitReg$out)[c(1:(nfac*9))],9,nfac)
  floads.list[[count]] <- floads
  hess <- fitReg$hess

  #zeros <- sum(floads < .0001)
  #floads2 <- GPForth(floads)$loadings
  resids <- matrix(coef(fitReg$out)[c((nfac*9 +1):(nfac*9 + 9))],9,1)

  ImpCov <- floads%*%t(floads) + vec2diag(resids)

  fit = log(det(ImpCov)) + trace(cov.train %*% ginv(ImpCov)) - log(det(cov.train))  - 9
  chisq = 150*fit
  zeros <- sum(floads < .0001)
  #df=efa.lav@Fit@test[[1]]$df# + zeros;
  df=trace(hess)
  N=150

  ncp = d(chisq,df,N)
  RMSEA = rmsea(ncp,df)

  # from Kenny website
  k = as.numeric(fitMeasures(efa.lav)["npar"])
  k2 = k - zeros
  BIC = chisq + log(150)*((k2*(k2+1)/2) - df)

  ret[[count]] = BIC
  ret2[[count]] = RMSEA


  # par influence
  optPars <- as.vector(coef(fitReg$out))
  #round(optPars <- as.vector(coef(fitReg$out)),3)
  # Lee & Wang (1996) Vector method
  U <- diag(sqrt(abs(optPars))/2)
  pca <- prcomp(U %*% hess %*% U,center=F,scale=F)
  #round(pca$sdev,3)
  omega1 <- pca$rotation[,1]
  omega18 <- pca$rotation[,18]
  par.up <- optPars + U %*% omega18
  par.low <- optPars + U %*% omega1

  influ[[count]] <- (par.up[1:9] - par.low[1:9])/optPars[1:9]



  fload.sum <- apply(abs(floads),2,sum)

  if(any(fload.sum < 1.2)){ # 1.2 works well
    nfac = nfac -1
    efa.lav = efaUnrotate(dat.train,nfac)
    coefs = as.numeric(coef(efa.lav))
    starts = as.numeric(coefs[-c((9*nfac + 1):(nfac*9 + 9))])
  }else {
    nfac = nfac
    coeff = as.numeric(coef(fitReg$out))
    starts = coeff
  }

}
#qplot(seq(1,150),unlist(ret),geom="dotplot")
#dotchart(unlist(ret),color=unlist(num.facs))
#plot(unlist(ret),color=unlist(num.facs))

par(mar = c(5, 4, 4, 4) + 0.3)
plot(unlist(ret),col=unlist(num.facs),ylab="BIC",xlab="Iteration")
par(new=T)
plot(unlist(ret2),col=unlist(num.facs),axes = F, xlab = "", ylab = "",pch=20)
axis(side=4)
mtext("RMSEA", side=4, line=3)




fload.mat <- matrix(unlist(floads.list),9,50,byrow=F)
fload.diff <- apply(fload.mat,1,sd)


fload.diffs <- matrix(0, 9,51)
for(i in 2:50){
  fload.diffs[,i] <- fload.mat[,i] - fload.mat[,i-1]
}


fload.rank = rank(apply(fload.mat,1,sd))

influ.mat <- matrix(unlist(influ),9,50,byrow=F)

influ.row <- rowMeans(influ.mat)
cor(fload.rank,influ.row)

# p.66 ESL
eigens <- eigen(cov.train)$vectors[,1]

cors <- rep(NA,50)
cors2 <- rep(NA,50)
cors3 <- rep(NA,50)
for(i in 1:50){
  cors[[i]] <- cor(fload.mat[,i],influ.mat[,i])
  cors2[[i]] <- cor(fload.diffs[,i+1],influ.mat[,i])
  cors3[[i]] <- cor(fload.diff,influ.mat[,i])
}
cors
#cors2
cors3




##### old ######
#plot(unlist(ret))
#par(new=TRUE)
#lines(unlist(ret2))
#lines(unlist(ret2))
#plot(unlist(num.facs),unlist(ret))

#palette()
# http://vis.supstat.com/2013/04/plotting-symbols-and-color-palettes/
