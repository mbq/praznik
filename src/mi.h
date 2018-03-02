SEXP C_mi(SEXP X,SEXP Y){
 int nt=omp_get_max_threads();
 int n,m,ny,*y,*nx,**x;
 struct ht *hta[nt];
 prepareInput(X,Y,R_NilValue,hta,&n,&m,NULL,&y,&ny,&x,&nx,nt);
 int *cXc=(int*)R_alloc(sizeof(int),n*nt);
 int *cYc=(int*)R_alloc(sizeof(int),n*nt);
 SEXP Ans=PROTECT(allocVector(REALSXP,m));
 double *score=REAL(Ans);

 #pragma omp parallel
 {
  int tn=omp_get_thread_num(),*cX=cXc+(tn*n),*cY=cYc+(tn*n),dy=0;
  struct ht *ht=hta[tn];
  #pragma omp for
  for(int e=0;e<m;e++){
   fillHt(ht,n,ny,y,nx[e],x[e],NULL,dy?NULL:cY,cX,0); dy=1;
   score[e]=miHt(ht,cY,cX);
  }
 }
 //Copy attribute names
 setAttrib(Ans,R_NamesSymbol,getAttrib(X,R_NamesSymbol));
 
 UNPROTECT(1);
 return(Ans);
}
