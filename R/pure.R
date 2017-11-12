#' @importFrom stats setNames
#' @importFrom utils tail

mergef<-function(x,y)
 factor(as.numeric(y)*length(levels(x))+as.numeric(x))

mutinfo<-function(x,y)
 .Call(C_getMi,factor(x),factor(y))

nmutinfo<-function(x,y)
 .Call(C_getNmi,factor(x),factor(y))

cmutinfo<-function(a,b,c){
 a<-factor(a)
 b<-factor(b)
 c<-factor(c)
 bc<-mergef(b,c)
 mutinfo(a,bc)-mutinfo(a,c)
}

pureMIM<-function(X,Y,k=3){
 apply(X,2,mutinfo,Y)->mim
 sort(mim,decreasing=TRUE)[1:k]->ans
 list(
  selection=names(ans),
  scores=setNames(ans,NULL)
 )
}

pureCMIM<-function(X,Y,k=3){
 X<-data.frame(X)
 ascores<-apply(X,2,mutinfo,Y)
 selection<-names(which.max(ascores))
 fscores<-max(ascores)
 scores<-rep(Inf,ncol(X))
 if(k>1) for(e in 1:(k-1)){
  factor(X[,tail(selection,1)])->w
  scores[colnames(X)!=tail(selection,1)]->scores
  X[,colnames(X)!=tail(selection,1),drop=FALSE]->X
  newScores<-apply(X,2,function(xx) cmutinfo(xx,Y,w))
  scores<-pmin(
   newScores,  
   scores
  )
  selection<-c(selection,names(which.max(scores)))
  fscores<-c(fscores,max(scores))
 }
 list(
  selection=selection,
  scores=setNames(fscores,NULL)
 )
}

pureJMIM<-function(X,Y,k=3){
 X<-data.frame(X)
 ascores<-apply(X,2,mutinfo,Y)
 selection<-names(which.max(ascores))
 fscores<-max(ascores)
 scores<-rep(Inf,ncol(X))
 if(k>1) for(e in 1:(k-1)){
  factor(X[,tail(selection,1)])->x
  scores[colnames(X)!=tail(selection,1)]->scores
  X[,colnames(X)!=tail(selection,1),drop=FALSE]->X
  newScores<-apply(X,2,function(xx) mutinfo(mergef(x,factor(xx)),Y))
  scores<-pmin(
   newScores,
   scores
  )
  if(max(scores)==0) break

  selection<-c(selection,names(which.max(scores)))
  fscores<-c(fscores,max(scores))
 }
 list(
  selection=selection,
  scores=setNames(fscores,NULL)
 )
}

pureNJMIM<-function(X,Y,k=3){
 X<-data.frame(X)
 ascores<-apply(X,2,mutinfo,Y)
 selection<-names(which.max(ascores))
 fscores<-max(ascores)
 scores<-rep(Inf,ncol(X))
 if(k>1) for(e in 1:(k-1)){
  factor(X[,tail(selection,1)])->x
  scores[colnames(X)!=tail(selection,1)]->scores
  X[,colnames(X)!=tail(selection,1),drop=FALSE]->X
  newScores<-apply(X,2,function(xx) nmutinfo(mergef(x,factor(xx)),Y))
  scores<-pmin(
   newScores,
   scores
  )
  if(max(scores)==0) break

  selection<-c(selection,names(which.max(scores)))
  fscores<-c(fscores,max(scores))
 }
 list(
  selection=selection,
  scores=setNames(fscores,NULL)
 )
}

pureJMI<-function(X,Y,k=3){
 X<-data.frame(X)
 ascores<-apply(X,2,mutinfo,Y)
 selection<-names(which.max(ascores))
 fscores<-max(ascores)
 scores<-rep(0,ncol(X))
 if(k>1) for(e in 1:(k-1)){
  factor(X[,tail(selection,1)])->x
  scores[colnames(X)!=tail(selection,1)]->scores
  X[,colnames(X)!=tail(selection,1),drop=FALSE]->X
  scores+apply(X,2,function(xx) mutinfo(mergef(x,factor(xx)),Y))->scores
  if(max(scores)==0) break
  selection<-c(selection,names(which.max(scores)))
  fscores<-c(fscores,max(scores))
 }
 list(
  selection=selection,
  scores=setNames(fscores,NULL)
 )
}

pureMRMR<-function(X,Y,k=3){
 X<-data.frame(X)
 rel<-apply(X,2,mutinfo,Y)
 red<-rep(0,ncol(X))
 selection<-names(which.max(rel))
 fscores<-max(rel)
 if(k>1) for(e in 1:(k-1)){
  factor(X[,tail(selection,1)])->x
  rel[colnames(X)!=tail(selection,1)]->rel
  red[colnames(X)!=tail(selection,1)]->red
  X[,colnames(X)!=tail(selection,1),drop=FALSE]->X

  apply(X,2,function(xx) mutinfo(x,factor(xx)))->nred
  red<-red+nred;
  scores<-rel-red/e

  selection<-c(selection,names(which.max(scores)))
  fscores<-c(fscores,max(scores))
 }
 list(
  selection=selection,
  scores=fscores
 )
}

pureDISR<-function(X,Y,k=3){
 X<-data.frame(X)
 ascores<-apply(X,2,mutinfo,Y)
 selection<-names(which.max(ascores))
 fscores<-max(ascores)
 rep(0,ncol(X))->scores
 if(k>1) for(e in 1:(k-1)){
  factor(X[,tail(selection,1)])->x
  scores[colnames(X)!=tail(selection,1)]->scores
  X[,colnames(X)!=tail(selection,1),drop=FALSE]->X
  
  scores+apply(X,2,function(xx) nmutinfo(mergef(x,factor(xx)),Y))->scores
  if(max(scores)==0) break

  selection<-c(selection,names(which.max(scores)))
  fscores<-c(fscores,max(scores))
 }
 list(
  selection=selection,
  scores=fscores
 )
}

