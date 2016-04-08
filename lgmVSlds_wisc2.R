# comparison of latent growth and latent difference score models

wisc <- read.table("/Users/RJacobucci/Documents/Github/random_code/wisc4vpe.dat")
names(wisc)<- c("V1","V2","V4","V6","P1","P2","P4", "P6", "Moeducat")

library(lavaan)
library(semPlot)

model <- ' i =~ 1*V1 + 1*V2 + 1*V4 + 1*V6
           s =~ 1*V1 + 2*V2 + 4*V4 + 6*V6 
#residuals equal
V1 ~~ resid*V1; V2 ~~ resid*V2; V4 ~~ resid*V4; V6 ~~ resid*V6;'
fit <- growth(model, data=wisc)
summary(fit)
fitMeasures(fit)

model2 <- ' i =~ 1*V1 + 1*V2 + 1*V4 + 1*V6
           s =~ 0*V1 + l1*V2 + l2*V4 + 5*V6 
#residuals equal
V1 ~~ resid*V1; V2 ~~ resid*V2; V4 ~~ resid*V4; V6 ~~ resid*V6;'
fit2 <- growth(model2, data=wisc)
summary(fit2)


lds_1 <- "

#latent variables 
lV1 =~ 1*V1
lV2 =~ 1*V2
lV4 =~ 1*V4
lV6 =~ 1*V6



#autoregressions
lV2 ~ 1*lV1; lV4 ~ 1*lV2; lV6 ~ 1*lV4

#change - delta; d
dV1 =~ 1*lV2; dV2 =~ 1*lV4; dV3 =~ 1*lV6

#intercept and slope
inV =~ 1*lV1;
#slope =~ alpha*dV1 + alpha*dV2 + alpha*dV3
# match lgm
slope =~ 1*dV1 + 2*dV2 + 2*dV3
#slope =~ 1.17*dV1 + alpha1*dV2 + alpha2*dV3

	#manifest means @0
	V1 ~ 0*1; V2 ~0*1; V4 ~ 0*1; V6 ~ 0*1

	#slope and intercept means
slope ~ 1; 
inV ~ 1;

#Latent variances and covariance
slope ~~ slope;
inV ~~ inV;
slope ~~ inV;

#means and vars @0
lV1 ~ 0*1; lV2 ~0*1; lV4 ~ 0*1; lV6 ~ 0*1
dV1 ~ 0*1; dV2 ~0*1; dV3 ~ 0*1

lV1 ~~ 0*lV1; lV2 ~~ 0*lV2; lV4 ~~ 0*lV4; lV6 ~~ 0*lV6
dV1 ~~ 0*dV1; dV2 ~~ 0*dV2; dV3 ~~ 0*dV3

#auVo-proportions
dV1 ~ beta*lV1; dV2 ~ beta*lV2; dV3 ~ beta*lV4;

#residuals equal
V1 ~~ resid*V1; V2 ~~ resid*V2; V4 ~~ resid*V4; V6 ~~ resid*V6;
"

fit.lds <- lavaan(lds_1, data=wisc)
summary(fit.lds)



#### plotting the lds #####
wisc.verb <- wisc[,c(1:4,9)]
ntot <- nrow(wisc.verb)    # total number of observations
wisc.verb.sel <- wisc.verb[sample(ntot, 50), ]



wisc.long <- reshape(wisc.verb, varying = c("V1", "V2", "V4", "V6"), v.names = "verbal",
                     times = c(1, 2, 4, 6), direction = "long")

wisc.long.sel <- reshape(wisc.verb.sel, varying = c("V1", "V2", "V4", "V6"),
                         v.names = "verbal", times = c(1, 2, 4, 6), direction = "long")
head(wisc.long)
names(wisc.long)[2] <- "grade"
names(wisc.long.sel)[2] <- "grade"


library(ggplot2)
(p = qplot(grade,verbal,group=id,data=wisc.long.sel,alpha=I(1/2),
      geom = c("line","point"),xlab = "Grade of Testing", ylab = "Verbal[t]")) 


# additive change (linear model, alpha=1*year)
alpha1 = 1; alpha2=2; alpha3 = 2
y1.1 = 19.824
y2.1 = y1.1 + alpha1*4.67
y3.1 = y2.1 + alpha2*4.67
y4.1 = y3.1 + alpha3*4.67
#lines(c(1,2,4,6),c(y1.1,y2.1,y3.1,y4.1),lwd=3,col="red")


df1 <- data.frame(x=c(1,2,4,6), y = c(y1.1,y2.1,y3.1,y4.1),id=c(999,999,999,999))
(p2 = p + geom_line(data=df1,aes(x=x,y=y,id=id),colour="red",size=1))




alpha1 = 1; alpha2=2; alpha3 = 2
beta= .432

y1.2 = 19.17
y2.2 = y1.2 + alpha1*-1.84 + beta*y1.2
y3.2 = y2.2 + alpha2*-1.84 + beta*y2.2
y4.2 = y3.2 + alpha3*-1.84 + beta*y3.2
#lines(c(1,2,4,6),c(y1.2,y2.2,y3.2,y4.2),lwd=3,col="red")

df2 <- data.frame(x=c(1,2,4,6), y = c(y1.2,y2.2,y3.2,y4.2),id=c(9999,9999,9999,9999))
(p3 = p2 + geom_line(data=df2,aes(x=x,y=y,id=id),colour="green",size=1))




######### try to do latent difference score model with nlme ########

library(nlme)
head(wisc.long)

# first, match growth model
nlme.out1 = lme(verbal ~ grade,random= ~ grade | id,wisc.long,method="ML")
summary(nlme.out1)
coef1 = coef(nlme.out1)



nlme.outDif = lme(verbal ~ grade* ,random= ~ grade | id,wisc.long,method="ML")
summary(nlme.outDif)



if(grade == 1){
 verbal ~ 1
}else if(grade > 1){
 verbal ~ grade
}
h_n1 = g_n1;
traject = h_n1;

DO t = 2 to 7;
h_n[t] = h_n[t-1] + (g_n2 + beta*h_n[t-1]);
IF grade = t+1 THEN traject = h_n[t];
END;







model2 <- ' i =~ 1*V1 + 1*V2 + 1*V4 + 1*V6
           s =~ 1*V1 + 2*V2 + 4*V4 + 6*V6 
;'
fit2 <- growth(model2, data=wisc)
summary(fit2)
fitMeasures(fit2)

nlme.out2 = lme(verbal ~ grade,random=list(id=pdSymm(~grade)),weights=varIdent(form=~1|grade),
                wisc.long,method="ML")
summary(nlme.out2)



model3 <- ' i =~ 1*V1 + 1*V2 + 1*V4 + 1*V6
           s =~ 1*V1 + s1*V2 + s3*V4 + 6*V6 
'
fit3 <- growth(model3, data=wisc)
summary(fit3)


nlme.out3 = lme(verbal ~ grade,
                random= ~ 1 | id,
                weights=varIdent(form=~1|grade),
                data=wisc.long,method="ML")
summary(nlme.out3)



### try time varying loadings
head(wisc.long)


wisc.long2 = cbind(wisc.long, age=round(wisc.long$grade + rnorm(nrow(wisc.long),5),0))


nlme.out02 = lme(verbal ~ 1 + age,random=list(verbal~age|id,~1|id),
                 wisc.long2,method="ML")
summary(nlme.out02)
coef1 = coef(nlme.out01)



# create comparison dataset
wisc44 = reshape(wisc.long2,
                #varying = c("verbal","age"),
                timevar = "grade",
                direction = "wide")

#wisc44[,"id",c("verbal.1","verbal.2","verbal.4","verbal.6"),c("age.1","age.2","age.3","age.4"))]

age = wisc44[,c("age.1","age.2","age.4","age.6")]
verbal = wisc44[,c("verbal.1","verbal.2","verbal.4","verbal.6")]
attach(age)
ageDev1 = age.2 - age.1
ageDev2 = age.4 - age.2
ageDev3 = age.6 - age.4

ageDev = cbind(ageDev1,ageDev2,ageDev3)
write.table(cbind(age,verbal,ageDev),"/Users/RJacobucci/Desktop/wisc44.dat",
            col.names=FALSE,row.names=FALSE,na=".")






wisc44 = data.frame(age,verbal)
colnames(wisc44) = c("age1","age2","age4","age6","V1","V2","V4","V6")
attach(wisc44)
wisc44$mult1 = age1*V1
wisc44$mult2 = age2*V2
wisc44$mult4 = age4*V4
wisc44$mult6 = age6*V6

lds_1 <- "

#latent variables 
lV1 =~ 1*mult1
lV2 =~ 1*mult2
lV4 =~ 1*mult4
lV6 =~ 1*mult6



#autoregressions
lV2 ~ 1*lV1; lV4 ~ 1*lV2; lV6 ~ 1*lV4

#change - delta; d
dV1 =~ 1*lV2; dV2 =~ 1*lV4; dV3 =~ 1*lV6

#intercept and slope
inV =~ 1*lV1;
#slope =~ alpha*dV1 + alpha*dV2 + alpha*dV3
# match lgm
slope =~ 1*dV1 + 2*dV2 + 2*dV3
#slope =~ 1.17*dV1 + alpha1*dV2 + alpha2*dV3

#manifest means @0
mult1 ~ 0*1; mult2 ~0*1; mult4 ~ 0*1; mult6 ~ 0*1

#slope and intercept means
slope ~ start(-.8)*1; 
inV ~ start(4)*1;

#Latent variances and covariance
slope ~~ start(0.6)*slope;
inV ~~ start(2.2)*inV;
slope ~~ start(-0.3)*inV;

#means and vars @0
lV1 ~ 0*1; lV2 ~0*1; lV4 ~ 0*1; lV6 ~ 0*1
dV1 ~ 0*1; dV2 ~0*1; dV3 ~ 0*1

lV1 ~~ 0*lV1; lV2 ~~ 0*lV2; lV4 ~~ 0*lV4; lV6 ~~ 0*lV6
dV1 ~~ 0*dV1; dV2 ~~ 0*dV2; dV3 ~~ 0*dV3

#auVo-proportions
#dV1 ~ beta*lV1; dV2 ~ beta*lV2; dV3 ~ beta*lV4;

#residuals equal
mult1 ~~ start(35)*resid*mult1; mult2 ~~ resid*mult2; 
mult4 ~~ resid*mult4; mult6 ~~ resid*mult6;
"

fit.lds <- lavaan(lds_1, data=wisc44,missing="fiml")
summary(fit.lds)
