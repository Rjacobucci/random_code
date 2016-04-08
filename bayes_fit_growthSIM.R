library(blavaan)


# three simulated dataset
# 1. no growth
# 2. linear growth
# 3. quadratic growth

# vary sample size simulated

mod.noGrowth <-"
i =~ 1*t1 + 1*t2 + 1*t3 + 1*t4
s =~ 0*t1 + 0*t2 + 0*t3 + 0*t4
i~5*1
s~0*1
s~~1*s
i~~1*i
t1~~1*t1
t2~~1*t2
t3~~1*t3
t4~~1*t4
i~~0*s
"

mod.linearGrowth <-"
i =~ 1*t1 + 1*t2 + 1*t3 + 1*t4
s =~ 0*t1 + 1*t2 + 2*t3 + 3*t4
i~5*1
s~0*1
s~~1*s
i~~1*i
t1~~1*t1
t2~~1*t2
t3~~1*t3
t4~~1*t4
i~~0*s
"

mod.quadGrowth <-"
i =~ 1*t1 + 1*t2 + 1*t3 + 1*t4
s =~ 0*t1 + 1*t2 + 2*t3 + 3*t4
s2 =~ 0*t1 + 1*t2 + 4*t3 + 9*t4
i~5*1
s~0*1
s2~0*1
s~~1*s
i~~1*i
s2~~1*s2
t1~~1*t1
t2~~1*t2
t3~~1*t3
t4~~1*t4
i~~0*s
i~~0*s2
s~~0*s2
"

 

count= 0
samps = c(50,100,200,500,5000)
iters=20
mods = list(mod.noGrowth,mod.linearGrowth,mod.quadGrowth)

fit.ret = data.frame(matrix(NA,iters*length(samps)*length(mods),22))
colnames(fit.ret) = c("samp.size","mod","bic1","dic1","waic1","looic1","jags_dic1",
                       "bic2","dic2","waic2","looic2","jags_dic2",
                       "bic3","dic3","waic3","looic3","jags_dic3",
                       "bic4","dic4","waic4","looic4","jags_dic4")

for(i in 1:iters){
 for(j in 1:length(samps)){
  for(k in 1:length(mods)){
   count = count + 1

fit.ret[count,"samp.size"] = samps[[j]]
fit.ret[count,"mod"] = k

#dat1 <- simulateData(mod.noGrowth,sample.nobs=100,model.type="lavaan")
#dat2 <- simulateData(mod.linearGrowth,sample.nobs=100,model.type="lavaan")
#dat3 <- simulateData(mod.quadGrowth,sample.nobs=100,model.type="lavaan")
dat <- simulateData(mods[[k]],sample.nobs=samps[[j]],model.type="lavaan")



mod1 <- ' i =~ 1*t1 + 1*t2 + 1*t3 + 1*t4
          s =~ 0*t1 + 0*t2 + 0*t3 + 0*t4 '
fit1 <- bgrowth(mod1, data=dat,sample=100000,burnin=10000, adapt=2000)

fits1= fitmeasures(fit1)
fit.ret[count,"bic1"] = fits1["bic"]
fit.ret[count,"dic1"] = fits1["dic"]
fit.ret[count,"waic1"] = fits1["waic"]
fit.ret[count,"looic1"] = fits1["looic"]
fit.ret[count,"jags_dic1"] = extract(fit1@external$runjags,"dic")



mod2 <- ' i =~ 1*t1 + 1*t2 + 1*t3 + 1*t4
s =~ 0*t1 + 1*t2 + 2*t3 + 3*t4 '
fit2 <- bgrowth(mod2, data=dat,sample=100000,burnin=10000, adapt=2000)

fits2 = fitmeasures(fit2)
fit.ret[count,"bic2"] = fits2["bic"]
fit.ret[count,"dic2"] = fits2["dic"]
fit.ret[count,"waic2"] = fits2["waic"]
fit.ret[count,"looic2"] = fits2["looic"]
fit.ret[count,"jags_dic2"] = extract(fit2@external$runjags,"dic")



mod3 <- ' i =~ 1*t1 + 1*t2 + 1*t3 + 1*t4
          s =~ 0*t1 + 1*t2 + 2*t3 + 3*t4
          s2 =~ 0*t1 + 1*t2 + 4*t3 + 9*t4'
fit3 <- bgrowth(mod3, data=dat,sample=100000,burnin=10000, adapt=2000)

fits3 = fitmeasures(fit3)
fit.ret[count,"bic3"] = fits3["bic"]
fit.ret[count,"dic3"] = fits3["dic"]
fit.ret[count,"waic3"] = fits3["waic"]
fit.ret[count,"looic3"] = fits3["looic"]
fit.ret[count,"jags_dic3"] = extract(fit3@external$runjags,"dic")



mod4 <- ' i =~ 1*t1 + 1*t2 + 1*t3 + 1*t4
          s =~ 0*t1 + l1*t2 + l2*t3 + 1*t4'
fit4 <- bgrowth(mod4, data=dat,sample=100000,burnin=10000, adapt=2000)

fits4 = fitmeasures(fit4)
fit.ret[count,"bic4"] = fits4["bic"]
fit.ret[count,"dic4"] = fits4["dic"]
fit.ret[count,"waic4"] = fits4["waic"]
fit.ret[count,"looic4"] = fits4["looic"]
fit.ret[count,"jags_dic4"] = extract(fit4@external$runjags,"dic")


  }
 }
}