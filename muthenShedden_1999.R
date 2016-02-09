

## chapter 7
# example 7.17 Mixture CFA Modeling
library(boot)
library(mvtnorm)




n=500
K = 2
theta = diag(c(0.5,0.5,0.5,0.5,0.5))
lambda_y = matrix(c(0.8,0.8,0.8,0.8,0.8),5,1)
epsilon = rnorm(n) # error for y_i
A = matrix(c(-1,1),1,2)
c.fill = rbinom(n,1,0.5)
c <- cbind(c.fill,1-c.fill)

u = cbind(rbinom(n,1,0.5),rbinom(n,1,0.5),rbinom(n,1,0.5))
x = matrix(rnorm(n*3),n,3) # covariates
Gamma_eta = t(c(-1,0,1)) # regression coefficients from x to eta
Zeta = rnorm(n)
psi = 1
eta = rep(0,n)
 
for(i in 1:n){
   A_add = ifelse(c[i]==1,A[,1],A[,2])
   eta[i] = A_add + Gamma_eta %*% x[i,] + Zeta[i]
}
 
y = t(lambda_y %*% eta + epsilon)

# binary outcome variables
# 3 r variables

lambda_u = matrix(c(0.1,0.5,
                    0.5,0.1,
                    0.1,0.5),3,2)

tau = matrix(0,n,3)
for(i in 1:n){
  ind = c[i]
  tau[i,] = inv.logit(lambda_u[,ind] * c[i])
}

### not sure what tau_ijk is

pi = matrix(0,n,2) # class probabilities   Eq 6
a = cbind(rnorm(n),rnorm(n)) # class probability intercepts
Gamma_c =matrix(c(-1,1,
                   1,1,
                  -1,1),3,2)

for(i in 1:n){
  for(cc in 1:2){
   pi[i,cc] = inv.logit(a[i,cc] + Gamma_c[,cc] %*% x[i,])
  }
}



# make sure sums to 1
for(i in 1:n){
 pi[i,] = pi[i,]/sum(pi[i,])
}





count = 0
iters = 50

while(count < iters){

count = count + 1
print(count)
#---------------------------------
# --------------------------------
#      E-Step
# -------------------------------
# --------------------------------


# eq 9
yu_x = matrix(0,500,2)
for(i in 1:n){
 for(j in 1:3){
  for(k in 1:K){
   yu_x[i,k] = yu_x[i,k] + pi[i,k] * dmvnorm(as.vector(lambda_y %*% (A[k] + Gamma_eta %*% x[i, ])), rep(0,5),
                                             lambda_y %*% psi %*% t(lambda_y) + theta) * tau[i,j]
  }
 }
}


# eq 14 
# computing posterior probability, p_ik

p_ik = 1/yu_x

#for(i in 1:n){
# for(k in 1:K){
#  pi[i,K] * dmvnorm(as.vector(lambda_y %*%(A[k] + Gamma_eta %*% x[i,])),rep(0,5), 
#                    lambda_y %*% psi %*% t(lambda_y) + theta)
# }
#}



S_cc = matrix(0,2,2)
for(i in 1:n){
 S_cc = S_cc + 1/n * diag(p_ik[i,])
}


S_cx = matrix(0,2,3)
for(i in 1:n){
 S_cx = S_cx + 1/n * p_ik[i,] %*% t(x[i,])
}

S_cy = matrix(0,2,5)
for(i in 1:n){
 S_cy = S_cy + 1/n * p_ik[i,] %*% t(y[i,])
}


S_xx = cov(x)
S_yy = cov(y)
S_xy = cov(x,y)
S_yx = cov(y,x)

pi_avg = c(0,0)
for(i in 1:n){
  pi_avg = pi_avg + 1/n * p_ik[i,]
}


V = solve(t(lambda_y) %*% solve(theta) %*% lambda_y + solve(psi))

# not sure if S_xx calculation is correct
B = V %*% solve(psi) %*% (Gamma_eta %*% S_xx %*% t(Gamma_eta) + Gamma_eta %*% t(S_cx) %*% t(A) + 
           A %*% S_cx %*% t(Gamma_eta) + A %*% S_cc %*% t(A)) %*% solve(psi) %*% V

C = V %*% t(lambda_y) %*% solve(theta) %*% S_yy %*% solve(theta) %*% lambda_y %*% V

D = V %*% solve(psi) %*% Gamma_eta %*% S_xy %*% solve(theta) %*% lambda_y + V %*% solve(psi) %*% A %*% S_cy %*% solve(theta) %*% lambda_y

S_eta = V + B + C + D

## not sure what S_yp is -- try t(S_cy)
S_yeta = S_yy %*% solve(theta) %*% lambda_y + (S_yx %*% t(Gamma_eta) + t(S_cy) %*% t(A)) %*% solve(psi) %*% V
## guessing S_px is S_cx
S_etax = V %*% (solve(psi) %*% (Gamma_eta %*% S_xx + A %*% S_cx) + t(lambda_y) %*% solve(theta) %*% S_yx)

S_etac = V %*% (solve(psi) %*% (Gamma_eta %*% t(S_cx) + A %*% S_cc) + t(lambda_y) %*% solve(theta) %*% t(S_cy))



# -------------------------------
# -------------------------------
#           M-Step
# -------------------------------
# -------------------------------


#(t(S_etax) %*% S_etac) %*% solve(t(S_cx))

psi = S_eta + Gamma_eta %*% S_xx %*% t(Gamma_eta) + A %*% S_cc %*% t(A) - S_etax %*% t(Gamma_eta) - 
      Gamma_eta %*% t(S_etax) - S_etac %*% t(A) - A %*% t(S_etac) + Gamma_eta %*% t(S_cx) %*% t(A) +
      A %*% S_cx %*% t(Gamma_eta)
psi = 1 ### just keep set it at this for this example

lambda_y = S_yeta %*% solve(S_eta)

theta = diag(diag(S_yy - S_yeta %*% t(lambda_y) - lambda_y %*% t(S_yeta) + lambda_y %*% S_eta %*% t(lambda_y)))


# equation 29 ?????????

A = S_etac %*% (1/cov(c))
Gamma_eta = S_etax %*% (1/cov(x))

logu_c = 0
#12
for(i in 1:n){
 for(j in 1:ncol(tau)){
  for(k in 1:K){
   if(c[i,k] == 1)
    logu_c =  logu_c +  log(tau[i,j]) + (1-u[n,j])*log(1-tau[i,j])
  }
 }
}

logc_x = 0
for(i in 1:n){
 for(k in 1:K){
  add = ifelse(c[i,k] == 1,log(pi[i,k]),0)
  logc_x = logc_x + add
 }
}



logeta_cx = 0
for(i in 1:n){
 #A_c = ifelse(c[i,] == 1,A[1],A[2])
 logeta_cx = logeta_cx + dnorm(A %*% c[i,] + Gamma_eta %*% x[i,],0, psi,log=T)
}


logy_eta = 0
for(i in 1:n){
 logy_eta = logy_eta + dmvnorm(as.vector(lambda_y * eta[i]),rep(0,nrow(theta)), theta,log=TRUE)
}


logLik_c = logc_x + logeta_cx + logy_eta + logu_c
print(logLik_c)

}
