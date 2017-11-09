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

pureCMI<-function(X,Y,k=3){
 X<-data.frame(X)
 S<-factor(rep(1,nrow(X)))
 selection<-c()
 ascores<-c()
 for(e in 1:k){
  apply(X,2,cmutinfo,Y,S)->scores
  if(max(scores)==0) break
  sel<-names(which.max(scores))
  selection<-c(selection,sel)
  ascores<-c(ascores,max(scores))
  S<-mergef(S,factor(X[,sel]))
  X[,colnames(X)!=sel,drop=FALSE]->X
 }
 list(
  selection=selection,
  scores=setNames(ascores,NULL)
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

#Verifiable with BetaGamma(Gamma=0,Beta=beta)
pureMIFS<-function(X,Y,k=3,beta=1){
 if(beta==0) return(pureMIM(X,Y,k))
 X<-data.frame(X)
 ascores<-apply(X,2,mutinfo,Y)
 selection<-names(which.max(ascores))
 fscores<-max(ascores)
 if(k>1) for(e in 1:(k-1)){
  factor(X[,tail(selection,1)])->x
  ascores[colnames(X)!=tail(selection,1)]->ascores
  X[,colnames(X)!=tail(selection,1),drop=FALSE]->X

  apply(X,2,function(xx) mutinfo(x,factor(xx)))->scores
  ascores<-ascores-beta*scores

  selection<-c(selection,names(which.max(ascores)))
  fscores<-c(fscores,max(ascores))
 }
 list(
  selection=selection,
  scores=fscores
 )
}

pureBetaGamma<-function(X,Y,k=3,beta=1,gamma=1){
 if(gamma==0) return(pureMIFS(X,Y,k,beta))
 X<-data.frame(X)
 ascores<-apply(X,2,mutinfo,Y)
 selection<-names(which.max(ascores))
 fscores<-max(ascores)
 if(k>1) for(e in 1:(k-1)){
  factor(X[,tail(selection,1)])->x
  ascores[colnames(X)!=tail(selection,1)]->ascores
  X[,colnames(X)!=tail(selection,1),drop=FALSE]->X

  apply(X,2,function(xx) mutinfo(x,factor(xx)))->crossMi
  apply(X,2,function(xx) cmutinfo(x,factor(xx),Y))->crossCmi
  ascores<-ascores-beta*crossMi+gamma*crossCmi

  selection<-c(selection,names(which.max(ascores)))
  fscores<-c(fscores,max(ascores))
 }
 list(
  selection=selection,
  scores=fscores
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

puremRMR_D<-function(X,Y,k=3){
 X<-data.frame(X)
 jscores<-apply(X,2,mutinfo,Y)
 bscores<-rep(0,ncol(X))
 selection<-names(which.max(jscores))
 fscores<-max(jscores)
 if(k>1) for(e in 1:(k-1)){
  factor(X[,tail(selection,1)])->x
  jscores[colnames(X)!=tail(selection,1)]->jscores
  bscores[colnames(X)!=tail(selection,1)]->bscores
  X[,colnames(X)!=tail(selection,1),drop=FALSE]->X

  apply(X,2,function(xx) mutinfo(x,factor(xx)))->scores
  bscores<-bscores+scores
  ascores<-jscores-bscores/e
  
  selection<-c(selection,names(which.max(ascores)))
  fscores<-c(fscores,max(ascores))
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

#Note: CondMI can stop at zero; others will either return nothing (for 0 in the initial MI scan) or k. Methods with negative scores will be even able to jump through zero.
