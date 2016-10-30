library(lavaan)

HS.model <- ' visual =~ x1 + x2 + x3
textual =~ x4 + x5 + x6
speed =~ x7 + x8 + x9 '
fit <- cfa(HS.model, data=HolzingerSwineford1939, se="standard",information="observed")
str(fit)
vcov = round(fit@vcov$vcov,4)
summary(fit)

lavInspect(fit,"se")

library(Matrix)

hess = lavInspect(fit,"hessian")
rankMatrix(hess) # equals number of parameters



info = lavInspect(fit,"information")
info.o = lavInspect(fit,"information.observed")

grad = lavInspect(fit,"gradient")


lavInspect(fit,"information.first.order")

lavInspect(fit,"inverted.information")


lavInspect(fit,"gradient")

round((1/301)* solve(hess),4)


1/3.59


(1/301) * info.o


sqrt(diag((2/301) * hess))


grad2 = grad %*% t(grad)

four = solve(hess) %*% grad2 %*% solve(hess)



