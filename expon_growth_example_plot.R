
time = c(0,1,2,3,4,5,6)
first.deriv = function(mu,time){
 1 - exp(-mu * time)
}

first.deriv(0,4)

second.deriv = function(mu1,mu2,time){
 mu2*time * exp(-mu1*time)
}

first.deriv(0,4)
second.deriv(0,0,4)

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
