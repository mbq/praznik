SEXP C_MI(SEXP X,SEXP Y,SEXP K){
 int n,k,m,ny,*y,*nx,**x;
 struct ht *ht;
 prepareInput(X,Y,K,&ht,&n,&m,&k,&y,&ny,&x,&nx);

 int *cX=(int*)R_alloc(sizeof(int),n);
 int *cY=(int*)R_alloc(sizeof(int),n);
 
 SEXP Ans; PROTECT(Ans=allocVector(REALSXP,m));
 double *mi=REAL(Ans);

 for(int e=0;e<m;e++){
  fillHt(ht,n,ny,y,nx[e],x[e],NULL,e?NULL:cY,cX,0);
  mi[e]=miHt(ht,cY,cX);
 }

 UNPROTECT(1);
 return(Ans);
}

