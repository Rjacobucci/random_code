library(lavaan)
library(caret)

samp = 500
iters=5000
count=0
val = rep(NA,iters)
type1.1 = rep(NA,iters)
type1.2 = rep(NA,iters)

for(i in 1:iters){
count=count+1
mod <- "
y ~ 0*x
y~~0.5*z
#x ~~0.5*z
"
dat <- simulateData(mod,sample.nobs=samp)

ids = sample(1:nrow(dat),nrow(dat)*0.5,prob=pnorm(dat$z))

dat$group = rep(0,nrow(dat))
dat[ids,"group"] = 1

#glm.out = glm(group ~ z,dat,family="binomial")
#wgts = 1/predict(glm.out,type="response")

#glm.out = glmnet(as.matrix(dat[,c("x","z")]),dat$group,family="binomial")
glm.out2 = cv.glmnet(as.matrix(dat[,c("x","z")]),dat$group,family="binomial",alpha=0)
wgts=predict(glm.out2,as.matrix(dat[,c("x","z")]),glm.out2$lambda.min,type="response")


out = lm(y ~ x + z,dat[ids,])
#summary(out)

coef1 = summary(out)$coefficients[2,"Pr(>|t|)"]

out2 = lm(y ~ x,dat[ids,],weights=wgts[ids])
#summary(out2)


coef2 = summary(out2)$coefficients[2,"Pr(>|t|)"]

val[count] = as.numeric(coef1 < coef2)
type1.1[count] =  as.numeric(coef1 < 0.05)
type1.2[count] =  as.numeric(coef2 < 0.05)
}
summary(val);summary(type1.1);summary(type1.2)
