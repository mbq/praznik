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

test_that("mi works like pure mi",{
 expect_equal(
  apply(X,2,mutinfo,Y),
  .Call(C_mi,X,Y)
 )
})

test_that("cmi works like pure cmi",{
 Z<-factor((1:150)%%7)
 expect_equal(
  apply(X,2,condmutinfo,Y,Z),
  cmi(X,Y,Z)
 )
})
