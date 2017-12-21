SEXP C_MI(SEXP X,SEXP Y,SEXP K){
 int nt=omp_get_max_threads();
 int n,k,m,ny,*y,*nx,**x;
 struct ht *hta[nt];
 prepareInput(X,Y,K,hta,&n,&m,&k,&y,&ny,&x,&nx,1);

 int *cXc=(int*)R_alloc(sizeof(int),n*nt);
 int *cYc=(int*)R_alloc(sizeof(int),n*nt);
 
 SEXP Ans; PROTECT(Ans=allocVector(REALSXP,m));
 double *mi=REAL(Ans);
 
 #pragma omp parallel
 {
  int tn=omp_get_thread_num(),*cX=cXc+(tn*n),*cY=cYc+(tn*n),dy=0;
  struct ht *ht=hta[tn];
  #pragma omp for
  for(int e=0;e<m;e++){
   fillHt(ht,n,ny,y,nx[e],x[e],NULL,dy?NULL:cY,cX,0); dy=1;
   mi[e]=miHt(ht,cY,cX);
  }
 }

 UNPROTECT(1);
 return(Ans);
}

