context("test-native-pure.R")

makeIn<-function(){
 #Generate iris with nonsense data, with rows and cols in random order
 is<-iris[sample(nrow(iris),nrow(iris)),]
 X<-cbind(is[,-5],apply(is[,rep(1:4,each=3)],2,sample))
 names(X)<-c(names(is)[-5],sprintf("Nonsense%d",1:(ncol(X)-4)))
 X[,sample(ncol(X),ncol(X))]->X
 Y<-is$Species

 #Make it discrete
 X<-data.frame(apply(X,2,cut,10))

 #Value of k for all tests
 K<-7
 list(X=X,Y=Y,k=K)
}

fuzz<-function(algo,seed){
 set.seed(seed)
 makeIn()->input
 test_that(sprintf("Native %s works like pure %s, seed %d",algo,algo,seed),{
  do.call(sprintf("pure%s",algo),input)->pure
  do.call(algo,input)->native
  expect_equal(pure$scores,native$scores)
  expect_equal(pure$selection,native$selection)
 })
}

for(seed in 1:5)
 for(algo in c("MIM","JMIM","NJMIM","JMI","DISR","CMI","CMIM","MRMR"))
  fuzz(algo,seed)


