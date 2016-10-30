
# read in class probabilities
# data_cprob.dat

dat <- read.table("/Users/RJacobucci/Documents/Github/semtree_mixture_compare/ecls-k_analyses/data_cprob.dat")
head(dat)

probs1 <- dat$V16
probs2 <- dat$V17
class1 <- dat$v18
class2 <- dat$v19


library(nnet)

wgt = probs1*probs2



out <- multinom(V19 ~ V8 + V9 + V10 + V11 + V12 + V13 + V14 + V15,dat)
summary(out)
38
-0.185/-0.192
0.094/0.098



library(mlogit)
