SEXP C_MIM(SEXP X,SEXP Y,SEXP K){
 int n,k,m,ny,*y,*nx,**x;
 struct ht *ht;
 prepareInput(X,Y,K,&ht,&n,&m,&k,&y,&ny,&x,&nx);

 int *cX=(int*)R_alloc(sizeof(int),n);
 int *cY=(int*)R_alloc(sizeof(int),n);
 
 double *score;
 int *idx;
 SEXP Ans; PROTECT(Ans=makeAns(k,&score,&idx));

 for(int e=0;e<k;e++){
  score[e]=-INFINITY; idx[e]=0;
 }

 for(int e=0;e<m;e++){
  fillHt(ht,n,ny,y,nx[e],x[e],NULL,e?NULL:cY,cX,0);
  double nmi=miHt(ht,cY,cX);

  //Insertion, since we need stable sort in principle
  if(score[k-1]>nmi) continue;
  int ee=k-2;
  for(;ee>=0 && score[ee]<nmi;ee--){
   score[ee+1]=score[ee]; idx[ee+1]=idx[ee];
  }
  score[ee+1]=nmi; idx[ee+1]=e+1;
 }

 UNPROTECT(1);
 return(Ans);
}

