# comparison of latent growth and latent difference score models
library(blavaan)


dat <- read.table("/Users/RJacobucci/Documents/Github/semtree_mixture_compare/ecls-k_analyses/ecls_mix_semtree_20160210.dat",
                  na.strings=".")

colnames(dat) <- c("id", "fine", "gross", "learn", "control", "interp", "ext",
                   "int", "math1", "gk1", "read1","read2","read3","read4",
                   "read5","read6","read7")
dat2 = dat[,-c(9)]

library(lavaan)
library(semPlot)


model <- ' i =~ 1*read1 + 1*read2 + 1*read3 + 1*read4 + 1*read5 + 1*read6 + 1*read7
           s =~ 0*read1 + s1*read2 + s2*read3 + s3*read4 + s4*read5 + s5*read6 + 1*read7
#residuals equal
#read1 ~~ resid*read1; read2 ~~ resid*read2; read3~~resid*read3;
#read4 ~~ resid*read4; read5~~resid*read5; read6 ~~ resid*read6;
#read7~~resid*read7
'
fit <- growth(model, data=dat2)
summary(fit)


model2 <- ' i =~ 1*V1 + 1*V2 + 1*V4 + 1*V6
           s =~ 0*V1 + l1*V2 + l2*V4 + 5*V6 
#residuals equal
V1 ~~ resid*V1; V2 ~~ resid*V2; V4 ~~ resid*V4; V6 ~~ resid*V6;'
fit2 <- growth(model2, data=wisc)
summary(fit2)


lds_1 <- "

#latent variables 
lx1 =~ 1*read1
lx2 =~ 1*read2
lx3 =~ 1*read3
lx4 =~ 1*read4
lx5 =~ 1*read5
lx6 =~ 1*read6
lx7 =~ 1*read7


#autoregressions
lx2 ~ 1*lx1; lx3 ~ 1*lx2; lx4 ~ 1*lx3; lx5 ~ 1*lx4;
lx6 ~ 1*lx5; lx7 ~ 1*lx6

#change - delta; d
dx1 =~ 1*lx2; dx2 =~ 1*lx3; dx3 =~ 1*lx4
dx4 =~ 1*lx5; dx5 =~ 1*lx6; dx6 =~ 1*lx7

#intercept and slope
int =~ 1*lx1;

slope =~ 1*dx1 + 1*dx2 + 1*dx3 + 1*dx4 + 1*dx5 + 1*dx6


#manifest means @0
read1 ~ 0*1; read2 ~0*1; read3 ~ 0*1; read5 ~ 0*1
read6 ~ 0*1; read7 ~ 0*1

#slope and intercept means
slope ~ start(3)*1; 
int ~ start(33)*1;

#Latent variances and covariance
slope ~~ start(12)*slope;
int ~~ start(122)*int;
slope ~~ start(-33)*int;

#means and vars @0
lx1 ~ 0*1; lx2 ~0*1;lx3 ~ 0*1; lx4 ~ 0*1; 
lx5 ~0*1; lx6 ~ 0*1; lx7 ~ 0*1

dx1 ~ 0*1; dx2 ~0*1; dx3 ~ 0*1; dx4 ~ 0*1
dx5 ~ 0*1; dx6 ~ 0*1

lx1 ~~ 0*lx1; lx2 ~~ 0*lx2; lx3 ~~ 0*lx3; lx4 ~~0*lx4
lx5 ~~ 0*lx5; lx6 ~~ 0*lx6; lx7 ~~ 0*lx7

dx1 ~~ 0*dx1; dx2 ~~ 0*dx2; dx3 ~~ 0*dx3
dx4 ~~ 0*dx4; dx5 ~~ 0*dx5;dx6 ~~ 0*dx6;

#auto-proportions
#dx1 ~ start(0.5)*beta*lx1; dx2 ~ start(0.5)*beta*lx2; dx3 ~ start(0.5)*beta*lx3;
#dx4 ~ start(0.5)*beta*lx4; dx5 ~ start(0.5)*beta*lx5; dx6 ~ start(0.5)*beta*lx6;

#residuals equal
#read1 ~~ start(180)*resid*read1;read2 ~~ resid*read2;read3 ~~ resid*read3;
#read4 ~~ resid*read4;read5 ~~ resid*read5;read6 ~~ resid*read6;
#read7 ~~ resid*read7;
read1 ~~ read1;read2 ~~ read2;read3 ~~ read3;
read4 ~~ read4;read5 ~~ read5;read6 ~~ read6;read7 ~~ read7;
"

fit.lds <- lavaan(lds_1, data=dat2)
#fit.bla = blavaan(lds_1,data=dat2)
summary(fit.lds)



#### plotting the lds #####
ecls <- dat2[,c(10:16)]
ntot <- nrow(ecls)    # total number of observations
ecls.sel <- ecls[sample(ntot, 50), ]



ecls.long <- reshape(ecls, varying = c("read1","read2","read3","read4","read5","read6","read7"),
                     v.names = "read",
                     times = c(0,.5,1,1.5,3.5,5.5,8.5), direction = "long")

ecls.long.sel <- reshape(ecls.sel, varying = c("read1","read2","read3","read4","read5","read6","read7"),
                         v.names = "read", times = c(0,.5,1,1.5,3.5,5.5,8.5), direction = "long")



library(ggplot2)
(pp = qplot(time,read,group=id,data=ecls.long.sel,alpha=I(1/2),
      geom = c("line","point"),xlab = "Grade of Testing", ylab = "Read[t]")) 


# additive change (linear model, alpha=1*year)
alpha = 1
y1.1 = 34.6
y2.1 = y1.1 + alpha*23.7
y3.1 = y2.1 + alpha*23.7
y4.1 = y3.1 + alpha*23.7
y5.1 = y4.1 + alpha*23.7
y6.1 = y5.1 + alpha*23.7
y7.1 = y6.1 + alpha*23.7

df1 <- data.frame(x=c(0,.5,1,1.5,3.5,5.5,8.5), y = c(y1.1,y2.1,y3.1,y4.1,y5.1,y6.1,y7.1),
                  id=c(999,999,999,999,999,999,999))
(p2 = pp + geom_line(data=df1,aes(x=x,y=y,id=id),colour="red",size=1))



alpha = 1
beta= .227

y1.2 = 36.5
y2.2 = y1.2 + alpha*5.5 + beta*y1.2
y3.2 = y2.2 + alpha*5.5 + beta*y2.2
y4.2 = y3.2 + alpha*5.5 + beta*y3.2
y5.2 = y4.2 + alpha*5.5 + beta*y4.2
y6.2 = y5.2 + alpha*5.5 + beta*y5.2
y7.2 = y6.2 + alpha*5.5 + beta*y6.2
#lines(c(1,2,4,6),c(y1.2,y2.2,y3.2,y4.2),lwd=3,col="red")

df2 <- data.frame(x=c(0,.5,1,1.5,3.5,5.5,8.5), y = c(y1.2,y2.2,y3.2,y4.2,y5.2,y6.2,y7.2),
                  id=c(999,999,999,999,999,999,999))
(p3 = p2 + geom_line(data=df2,aes(x=x,y=y,id=id),colour="green",size=1))



