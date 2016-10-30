library(glmnet);library(covTest);library(MASS)
library(dplyr);library(reshape2);library(data.table)

Ns <- c(40,60,100,500)
cond <- c(1,2,3)
count <- 0
niter=100
type1.res <- matrix(NA,niter*length(Ns)*length(cond),10)
type2.res <- matrix(NA,niter*length(Ns)*length(cond),10)






for(i in 1:niter){
  for(j in 1:length(Ns)){
    for(k in 1:length(cond)){
  
      count <- count + 1

M <- 31
N <- Ns[j]
dat <- data.frame(matrix(rnorm(N*M,mean=0,sd=1), N, M))

colnames(dat) <- c("y","x1","x2","x3","x4","x5",
                   "x6","x7","x8","x9","x10",
                   "x11","x12","x13","x14","x15",
                   "x16","x17","x18","x19","x20",
                   "x21","x22","x23","x24","x25",
                   "x26","x27","x28","x29","x30")

u <- rnorm(N)
r1 <- r2 <- r3 <- r4 <- r5 <- r6 <- 0

if(cond[k] == 1){
  r1 <- 0.3
  vars <- 1:5
}else if(cond[k] == 2){
  r1 <- 1
  vars <- 1:5
}else if(cond[k] == 3){
  r1 <- 1
  r2 <- 0.3
  vars <- c(1:10)
}

comb <- c(rep(r1,5),rep(r2,5),rep(r3,5),rep(r4,5),rep(r5,5),rep(r6,5))

dat$y <- r1*dat$x1 + r1*dat$x2 + r1*dat$x3 + r1*dat$x4 + r1*dat$x5 + 
     r2*dat$x6 + r2*dat$x7 + r2*dat$x8 + r2*dat$x9 + r2*dat$x10 + 
     r3*dat$x11 + r3*dat$x12 + r3*dat$x13 + r3*dat$x14 + r3*dat$x15 + 
     r4*dat$x16 + r4*dat$x17 + r4*dat$x18 + r4*dat$x19 + r4*dat$x20 + 
     r5*dat$x21 + r5*dat$x22 + r5*dat$x23 + r5*dat$x24 + r5*dat$x25 + 
     r6*dat$x26 + r6*dat$x27 + r6*dat$x28 + r6*dat$x29 + r6*dat$x30 +
     u



out <- lm(y ~ ., dat)
p_s <- summary(out)$coefficients[,"Pr(>|t|)"] < 0.05
p_s <- p_s[-1]

type1.1 <- sum((p_s == TRUE) & (comb == 0))
type2.1 <- sum((p_s == FALSE) & (comb != 0))

# bonferroni correction

p_s2 <- summary(out)$coefficients[,"Pr(>|t|)"] < 0.0017083
p_s2 <- p_s2[-1]

type1.2 <- sum((p_s2 == TRUE) & (comb == 0))
type2.2 <- sum((p_s2 == FALSE) & (comb != 0))

# lasso

x <- data.matrix(dat[,2:31]); y <- data.matrix(dat[,1])
lasso <- cv.glmnet(x,y,intercept=FALSE)

lambda1 <- lasso$lambda.min
lambda2 <- lasso$lambda.1se

lam.coef1 <- coef(lasso,lambda1) == 0
lam.coef2 <- coef(lasso,lambda2) == 0

lam.coef1 <- lam.coef1@x[2:31]
lam.coef2 <- lam.coef2@x[2:31]

#nwrong2 <- 30 - sum(lam.coef1 == c(comb == 0))
#nwrong3 <- 30 - sum(lam.coef2 == c(comb == 0))

type1.3 <- sum((lam.coef1 == FALSE) & (comb == 0))
type2.3 <- sum((lam.coef1 == TRUE) & (comb != 0))

type1.4 <- sum((lam.coef2 == FALSE) & (comb == 0))
type2.4 <- sum((lam.coef2 == TRUE) & (comb != 0))


# covTest


a=lars(x,y,intercept=FALSE)
cov.out <- covTest(a,x,y)
res1 <- cov.out$results

val <- 1:30 %in% res1[,1]

if(any(val == FALSE)){
  type1.5 <- type1.6 <- type2.5 <- type2.6 <- NA
}else{
  p_s3 <- cov.out$results[,"P-value"] < 0.05
  
  
  type1.5 <- sum((p_s3 == TRUE) & (comb[res1[,1]] == 0))
  type2.5 <- sum((p_s3 == FALSE) & (comb[res1[,1]] != 0))
  
  # bonferroni correction
  
  p_s4 <- cov.out$results[,"P-value"] < 0.0017083
  
  
  type1.6 <- sum((p_s4 == TRUE) & (comb[res1[,1]] == 0))
  type2.6 <- sum((p_s4 == FALSE) & (comb[res1[,1]] != 0))
}


# stepwise selection

step.out <- stepAIC(out,direction="both",trace=0)

p_s5 <- summary(step.out)$coefficients[,"Pr(>|t|)"] < 0.05
p_s5 <- p_s5[-1]

variables <- rownames(summary(step.out)$coefficients)[-1]

inds <- as.numeric(gsub("x","",variables))

type1.7 <- sum((p_s5 == TRUE) & (comb[inds] == 0))
type2.7 <- sum((p_s5 == FALSE) & (comb[inds] != 0))

# bonferroni correction

p_s6 <- summary(step.out)$coefficients[,"Pr(>|t|)"] < 0.0017083
p_s6 <- p_s6[-1]

type1.8 <- sum((p_s6 == TRUE) & (comb[inds] == 0))
type2.8 <- sum((p_s6 == FALSE) & (comb[inds] != 0))


type1.res[count,] <- c(Ns[j],k,type1.1,type1.2,type1.3,type1.4,type1.5,type1.6,type1.7,type1.8)
type2.res[count,] <- c(Ns[j],k,type2.1,type2.2,type2.3,type2.4,type2.5,type2.6,type2.7,type2.8)



rm("x");rm("y");rm("dat")
print((count/c(niter*length(Ns)*length(cond)))*100)
    }
  }
}

type1.res2 <- na.omit(type1.res)
type2.res2 <- na.omit(type2.res)
type1.res2 <- data.frame(type1.res2)
type2.res2 <- data.frame(type2.res2)
names(type1.res2) <- c("N","mod","lmP","lmPbonf","lassoMin","lasso1SE","lassoP","lassoPbonf","stepP","stepPbonf")
names(type2.res2) <- c("N","mod","lmP","lmPbonf","lassoMin","lasso1SE","lassoP","lassoPbonf","stepP","stepPbonf")


library(data.table)
DT <- data.table(type1.res2)

(DT[, lapply(.SD, mean), by = N])
(DT[, lapply(.SD, mean), by = mod])

DT2 <- data.table(type2.res2)
(DT2[, lapply(.SD, mean), by = N])
(DT2[, lapply(.SD, mean), by = mod])

library(ggplot2)
qplot(x=lassoP,data=DT)

ggplot(type1.res, aes(N, fill=workshop ) ) +
  geom_bar(position="dodge")


