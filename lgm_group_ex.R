library(lavaan)


model <- ' i =~ 1*t1 + 1*t2 + 1*t3 + 1*t4
           s =~ 0*t1 + 1*t2 + 2*t3 + 3*t4 '

dat1 <- simulateData(model,model.type="growth")
dat1$group <- 1



model2 <- ' i =~ 1*t1 + 1*t2 + 1*t3 + 1*t4
           s =~ 3*t1 + 2*t2 + 1*t3 + 0*t4 '

dat2 <- simulateData(model2,model.type="growth")
dat2$group <- 2

comb = rbind(dat1,dat2)

model9 <- ' i =~ 1*t1 + 1*t2 + 1*t3 + 1*t4
           s =~ 0*t1 + 1*t2 + 2*t3 + 3*t4 '

fit <- growth(model9, data=dat1)
summary(fit)
fitmeasures(fit)
predict(fit)
