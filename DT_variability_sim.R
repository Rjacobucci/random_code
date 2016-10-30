library(evtree)

library(MASS);library(rpart)
data(Boston)
rpart.Boston <- rpart(medv ~., data=Boston)
#summary(rpart.Boston)
plot(rpart.Boston);text(rpart.Boston)


ev.out <- evtree(medv ~., data=Boston)
plot(ev.out)

mat <- matrix(NA,506,100)
mod.mat <- rep(NA,100)
for(i in 1:100){
  set.seed(i)
  rpart.Boston <- rpart(medv ~., data=Boston)
  val <- min(rpart.Boston$cptable[,"xerror"])
  ind <- which(rpart.Boston$cptable[,"xerror"] == val)[1]
  mod.mat[i] <- rpart.Boston$cptable[ind,"nsplit"]
  mat[,i] <- predict(rpart.Boston)
}



mat <- matrix(NA,506,100)
mod.mat <- rep(NA,100)
for(i in 1:100){
  set.seed(i)
  
  ids <- sample(1:nrow(Boston),nrow(Boston),replace=TRUE)
  data <- Boston[ids,]
  
  rpart.Boston <- rpart(medv ~., data=data)
  val <- min(rpart.Boston$cptable[,"xerror"])
  ind <- which(rpart.Boston$cptable[,"xerror"] == val)[1]
  mod.mat[i] <- rpart.Boston$cptable[ind,"nsplit"]
  mat[,i] <- predict(rpart.Boston)
}





library(tree)

Boston2 <- Boston[1:50,]

mat <- matrix(NA,50,100)
mod.mat <- rep(NA,100)
for(i in 1:100){
  set.seed(i)
  
  ids <- sample(1:nrow(Boston2),nrow(Boston2),replace=TRUE)
  data <- Boston2[ids,]
  
  tree.Boston <- tree(medv ~., data=data)
  mod.mat[i] <- summary(tree.Boston)$size
  mat[ids,i] <- predict(tree.Boston)
}

library(psych)
corFiml(mat)
cor(mat,use="pairwise.complete.obs")
mean.pred <- rowMeans(mat,na.rm=T)

tree.Boston <- tree(medv ~., data=Boston)
cor(mean.pred,predict(tree.Boston))





library(ISLR)
data(Caravan)
library(tree)

mat <- matrix(NA,nrow(Caravan),100)
mod.mat <- rep(NA,100)
for(i in 1:100){
  set.seed(i)
  
  ids <- sample(1:nrow(Caravan),nrow(Caravan),replace=TRUE)
  data <- Caravan[ids,]
  
  tree.Caravan <- tree(MOSTYPE ~., data=data,control=tree.control(nobs=nrow(Caravan),mindev=0.001))
  mod.mat[i] <- summary(tree.Caravan)$size
  mat[ids,i] <- predict(tree.Caravan)
}

library(psych)
#corFiml(mat)
cor(mat,use="pairwise.complete.obs");mod.mat
mean.pred <- rowMeans(mat,na.rm=T)

tree.Caravan <- tree(medv ~., data=Caravan)
cor(mean.pred,predict(tree.Caravan))
