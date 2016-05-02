library(blavaan)
library(snowfall)

# three simulated dataset
# 1. no growth
# 2. linear growth
# 3. quadratic growth

# vary sample size simulated






#iters=50
grid = expand.grid(samps,1:3)


sim_bayes_fit <- function(N,mod){
mod.num = mod
samp.num = N
 
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

mods = list(mod.noGrowth,mod.linearGrowth,mod.quadGrowth)

fit.ret = data.frame(matrix(NA,1,34))
colnames(fit.ret) = c("samp.size","mod",
                       "conv1","logl1","bic1","dic1","waic1","looic1","margloglik1","jags_dic1",
                       "conv2","logl2","bic2","dic2","waic2","looic2","margloglik2","jags_dic2",
                       "conv3","logl3","bic3","dic3","waic3","looic3","margloglik3","jags_dic3",
                       "conv4","logl4","bic4","dic4","waic4","looic4","margloglik4","jags_dic4")



fit.ret["samp.size"] = samp.num
fit.ret["mod"] = mod.num

#dat1 <- simulateData(mod.noGrowth,sample.nobs=100,model.type="lavaan")
#dat2 <- simulateData(mod.linearGrowth,sample.nobs=100,model.type="lavaan")
#dat3 <- simulateData(mod.quadGrowth,sample.nobs=100,model.type="lavaan")
dat <- simulateData(mods[[mod.num]],sample.nobs=samp.num,model.type="lavaan")



mod1 <- ' i =~ 1*t1 + 1*t2 + 1*t3 + 1*t4
          s =~ 0*t1 + 0*t2 + 0*t3 + 0*t4 '
fit1 = try(bgrowth(mod1, data=dat,sample=100000,burnin=10000, adapt=2000,
                   jagcontrol=list(silent.jags=TRUE)))
if(inherits(fit1, "try-error")){
 fit.ret["conv1"] = -9999
 fit.ret["logl1"] = -9999
 fit.ret["bic1"] = -9999
 fit.ret["dic1"] = -9999
 fit.ret["waic1"] = -9999
 fit.ret["looic1"] = -9999
 fit.ret["margloglik1"] = -9999
 fit.ret["jags_dic1"] = -9999
}else{
 fit.ret["conv1"] = all(fit1@external$runjags$psrf$psrf[,1] < 1.2)
 fits1 = fitmeasures(fit1)
 fit.ret["logl1"] = fits1["logl"]
 fit.ret["bic1"] = fits1["bic"]
 fit.ret["dic1"] = fits1["dic"]
 fit.ret["waic1"] = fits1["waic"]
 fit.ret["looic1"] = fits1["looic"]
 fit.ret["margloglik1"] = fits1["margloglik"]
 fit.ret["jags_dic1"] = extract(fit1@external$runjags,"dic")
}


mod2 <- ' i =~ 1*t1 + 1*t2 + 1*t3 + 1*t4
s =~ 0*t1 + 1*t2 + 2*t3 + 3*t4 '
fit2 = try(bgrowth(mod2, data=dat,sample=100000,burnin=10000, adapt=2000,
                   jagcontrol=list(silent.jags=TRUE)))
if(inherits(fit2, "try-error")){
 fit.ret["conv2"] = -9999
 fit.ret["logl2"] = -9999
 fit.ret["bic2"] = -9999
 fit.ret["dic2"] = -9999
 fit.ret["waic2"] = -9999
 fit.ret["looic2"] = -9999
 fit.ret["margloglik2"] = -9999
 fit.ret["jags_dic2"] = -9999
}else{
 fit.ret["conv2"] = all(fit2@external$runjags$psrf$psrf[,1] < 1.2)
 fits2 = fitmeasures(fit2)
 fit.ret["logl2"] = fits2["logl"]
 fit.ret["bic2"] = fits2["bic"]
 fit.ret["dic2"] = fits2["dic"]
 fit.ret["waic2"] = fits2["waic"]
 fit.ret["looic2"] = fits2["looic"]
 fit.ret["margloglik2"] = fits2["margloglik"]
 fit.ret["jags_dic2"] = extract(fit2@external$runjags,"dic")
}


mod3 <- ' i =~ 1*t1 + 1*t2 + 1*t3 + 1*t4
          s =~ 0*t1 + 1*t2 + 2*t3 + 3*t4
          s2 =~ 0*t1 + 1*t2 + 4*t3 + 9*t4'

fit3 = try(bgrowth(mod3, data=dat,sample=100000,burnin=10000, adapt=2000,
                   jagcontrol=list(silent.jags=TRUE)))
if(inherits(fit3, "try-error")){
 fit.ret["conv3"] = -9999
 fit.ret["logl3"] = -9999
 fit.ret["bic3"] = -9999
 fit.ret["dic3"] = -9999
 fit.ret["waic3"] = -9999
 fit.ret["looic3"] = -9999
 fit.ret["margloglik3"] = -9999
 fit.ret["jags_dic3"] = -9999
}else{
 fit.ret["conv3"] = all(fit3@external$runjags$psrf$psrf[,1] < 1.2)
 fits3 = fitmeasures(fit3)
 fit.ret["logl3"] = fits3["logl"]
 fit.ret["bic3"] = fits3["bic"]
 fit.ret["dic3"] = fits3["dic"]
 fit.ret["waic3"] = fits3["waic"]
 fit.ret["looic3"] = fits3["looic"]
 fit.ret["margloglik3"] = fits3["margloglik"]
 fit.ret["jags_dic3"] = extract(fit3@external$runjags,"dic")
}


mod4 <- ' i =~ 1*t1 + 1*t2 + 1*t3 + 1*t4
          s =~ 0*t1 + l1*t2 + l2*t3 + 1*t4'
fit4 = try(bgrowth(mod4, data=dat,sample=100000,burnin=10000, adapt=2000,
                   jagcontrol=list(silent.jags=TRUE)))
if(inherits(fit4, "try-error")){
 fit.ret["conv4"] = -9999
 fit.ret["logl4"] = -9999
 fit.ret["bic4"] = -9999
 fit.ret["dic4"] = -9999
 fit.ret["waic4"] = -9999
 fit.ret["looic4"] = -9999
 fit.ret["margloglik4"] = -9999
 fit.ret["jags_dic4"] = -9999
}else{
 fit.ret["conv4"] = all(fit4@external$runjags$psrf$psrf[,1] < 1.2)
 fits4 = fitmeasures(fit4)
 fit.ret["logl4"] = fits4["logl"]
 fit.ret["bic4"] = fits4["bic"]
 fit.ret["dic4"] = fits4["dic"]
 fit.ret["waic4"] = fits4["waic"]
 fit.ret["looic4"] = fits4["looic"]
 fit.ret["margloglik4"] = fits4["margloglik"]
 fit.ret["jags_dic4"] = extract(fit4@external$runjags,"dic")
}
}

samps = c(50,200,500,5000)
grid = expand.grid(samps,1)

sfStop()
sfInit(parallel=TRUE,cpus=4)
sfExport("sim_bayes_fit")
sfLibrary(blavaan); sfLibrary(lavaan)
sfLibrary(runjags)

out = system.time(sfApply(grid,1,function(x,y) sim_bayes_fit(x[1],x[2])))
