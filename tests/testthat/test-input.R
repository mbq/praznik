context("test-input.R")

test_that("crazy decision works",{
 data.frame(A=rep(letters[1:2],each=20))->X
 rep(c(TRUE,FALSE),each=20)->Y
 expect_equal(MIM(X,Y,1)$selection,"A")
})

test_that("crazy attributes work",{
 badfactor<-factor(c(rep(c("z","a"),each=7),letters))[1:14]
 data.frame(
  bool=rep(c(TRUE,FALSE),each=7),
  int=as.integer(rep(c(-37,33),each=7)),
  badfactor=badfactor
 )->X
 factor(rep(letters[1:2],each=7))->Y
 expect_equal(MIM(X,Y,3)$selection,c("bool","int","badfactor"))
})

test_that("X must be a data.frame",{
 expect_error(MIM(list(1:3),NULL,NULL),"X must be a data.frame")
})

test_that("X must be only reals, booleans, integers or factors",{
 Y<-c(TRUE,TRUE,FALSE,FALSE,FALSE)
 li<-data.frame(A=1:5)
 li$A<-list(1,1:2,1:3,1:4,1:5)
 badX<-list(
  char=data.frame(A=letters[1:5],stringsAsFactors=FALSE),
  realna=data.frame(A=c((1:4)*5.5,NA)),
  realinf=data.frame(A=c(1:4,Inf)),
  img=data.frame(A=1:5+3i),
  li=li
 )
 for(X in badX)
  expect_error(MIM(X,Y,1))
})

test_that("NAs and other quirks are caught",{
 Y<-iris$Species; Y[3]<-NA
 X<-iris[,-5]; X[12,3]<-NA
 expect_error(MIM(X,Y,1))
 X[12,3]<-Inf
 expect_error(MIM(X,Y,1))
 X[12,3]<-NaN
 expect_error(MIM(X,Y,1))
 X<-iris[,"Species",drop=FALSE]; X[17,1]<-NA
 expect_error(MIM(X,Y,1))
})
