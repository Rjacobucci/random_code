install.packages("blavaan", repos="http://faculty.missouri.edu/~merklee", type="source")
library(blavaan)

data("Demo.growth")
colnames(Demo.growth)[1:4] = c("X1","X2","X3","X4")

#LCS specification (constant change)
constantX.syntax <-'

#X

#latent variables
lX1 =~ 1*X1; lX2 =~ 1*X2; lX3 =~ 1*X3; lX4 =~ 1*X4;

#autoregressions
lX2 ~ 1*lX1; lX3 ~ 1*lX2; lX4 ~ 1*lX3;

#change - delta; d
dX1 =~ 1*lX2; dX2 =~ 1*lX3; dX3 =~ 1*lX4;

#intercept and slope
intX =~ 1*lX1;
slopeX =~ 1*dX1 + 1*dX2 + 1*dX3;

#residuals equal
X1 ~~ residX*X1; X2 ~~ residX*X2; X3 ~~ residX*X3; X4 ~~ residX*X4;


#manifest means @0
X1 ~ 0*1; X2 ~0*1; X3 ~ 0*1; X4 ~ 0*1;

#auto-proportions
dX1 ~ 0*lX1; dX2 ~ 0*lX2; dX3 ~ 0*lX3;

#slope and intercept means
slopeX ~ 1;
intX ~ 1;

#Latent variances and covariance
slopeX ~~ slopeX;
intX ~~ intX;
slopeX ~~ intX;

#means and vars @0
lX1 ~ 0*1; lX2 ~0*1; lX3 ~ 0*1; lX4 ~ 0*1;
dX1 ~ 0*1; dX2 ~0*1; dX3 ~ 0*1;

lX1 ~~ 0*lX1; lX2 ~~ 0*lX2; lX3 ~~ 0*lX3; lX4 ~~ 0*lX4;
dX1 ~~ 0*dX1; dX2 ~~ 0*dX2; dX3 ~~ 0*dX3;

'
constantX.run <- lavaan(constantX.syntax, data = Demo.growth)
summary(constantX.run)
fitmeasures(constantX.run)

bsem.out <- blavaan(constantX.syntax, data = Demo.growth)
summary(bsem.out)
fitmeasures(bsem.out)
