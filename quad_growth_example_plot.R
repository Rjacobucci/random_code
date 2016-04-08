
loads <- matrix(c(0,0,
                  1,1,
                  2,4,
                  3,9,
                  4,16,
                  5,25,
                  6,49),7,2,byrow=TRUE)

inst.vel = 60
accel = -5

expected.growth <- matrix(
 rep(t(c(0,0)), each=7)+
  rep(t(c(inst.vel,accel)), each=7)*loads, nrow=2, byrow=T)

expected.growth.comb <- matrix(
 rep(t(c(0)), each=7)+
  rep(t(c(inst.vel,accel)), each=7)*loads, nrow=2, byrow=T)
expected.growth.comb = colSums(expected.growth.comb)



plot(c(0,6), c(-120,120), xlab="Grade", ylab="Reading Score", type="n")
lines(0:6, expected.growth[1,], col="black", type="b", lw=3,lty=3)
lines(0:6, expected.growth[2,], col="red", type="b", lw=3,lty=1)
lines(0:6, expected.growth.comb, col="green", type="b", lw=3,lty=5)
