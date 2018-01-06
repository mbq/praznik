context("test-native-pure.R")

data.frame(apply(iris[,-5],2,cut,10))->X
X$const<-factor(rep(1,150))
X$tri<-factor(rep(1:3,50))
Y<-iris$Species
list(X=X,Y=Y,k=4)->input

for(algo in c("MIM","JMIM","NJMIM","JMI","DISR","CMIM","MRMR")){
 test_that(sprintf("Native %s works like pure %s",algo,algo),{
  do.call(sprintf("pure%s",algo),input)->pure
  do.call(algo,input)->native
  expect_equal(pure$score,setNames(native$score,NULL))
  expect_equal(pure$selection,names(native$selection))
 })
}

