library(poLCA)
data(carcinoma)

d = 4 # number of latent classes
Y = t(carcinoma)



# starting values

# pi
pi = c(.25,.25,.25,.25)
# theta jk, j is for each class, k to b is number of responses (7)
theta = matrix(0.5,d,nrow(Y))
# mu


pi = c(.444,.3544,.0281,.1735)
theta = matrix(
c(1,.98,.8588,.5867,1,.4771,1,
0,.139,0,0,.0593,0,0,
1,.3958,0,0,0,0,0,
.5025,1,0,.06,.78,0,.66),
,nrow(Y),d)

theta[theta==0] <- .1
theta[theta==1] <- .9
###### n denotes iteration

#### c^y is number of people with response vector the same

# p. 620 

# log lik
cy = 1
loglik = 0 
sum1 = 0
 for(j in 1:d){
   for(k in 1:ncol(theta)){
    yk1 = sum(Y[k,] ==1)  ##### check later, might be backwards
    yk2 = sum(Y[k,] ==2)
    val1 = theta[j,k]^yk2 * (1 - theta[j,k])^yk2
    if(k == 1){
     val2 = val1
    }else{
     val2 = val2*val1
    }
   }
   val5 = pi[j] * val2
   sum1 = sum1 + val5
 }
log(sum1)
 
 loglik = loglik 
  
}