SEXP C_mi(SEXP X,SEXP Y,SEXP Threads){
 int n,m,ny,*y,*nx,**x,nt;
 struct ht **hta;
 prepareInput(X,Y,R_NilValue,Threads,&hta,&n,&m,NULL,&y,&ny,&x,&nx,&nt);
 int *cXc=(int*)R_alloc(sizeof(int),n*nt);
 int *cYc=(int*)R_alloc(sizeof(int),n*nt);
 SEXP Ans=PROTECT(allocVector(REALSXP,m));
 double *score=REAL(Ans);

 #pragma omp parallel num_threads(nt)
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

SEXP C_ig(SEXP X,SEXP Y,SEXP Threads){
 int n,m,ny,*y,*nx,**x,nt;
 struct ht **hta;
 prepareInput(X,Y,R_NilValue,Threads,&hta,&n,&m,NULL,&y,&ny,&x,&nx,&nt);
 int *cXc=(int*)R_alloc(sizeof(int),n*nt);
 int *cYc=(int*)R_alloc(sizeof(int),n*nt);
 SEXP Ans=PROTECT(allocVector(REALSXP,m));
 double *score=REAL(Ans);

 //fillHtOne(*hta,n,y,NULL,0);
 double offset=0.;//oneMinusGiHt(*hta);

 #pragma omp parallel num_threads(nt)
 {
  int tn=omp_get_thread_num(),*cX=cXc+(tn*n),*cY=cYc+(tn*n),dy=0;
  struct ht *ht=hta[tn];
  #pragma omp for
  for(int e=0;e<m;e++){
   fillHt(ht,n,ny,y,nx[e],x[e],NULL,dy?NULL:cY,cX,0); dy=1;
   score[e]=igHt(ht,cY,offset);
  }
 }
 //Copy attribute names
 setAttrib(Ans,R_NamesSymbol,getAttrib(X,R_NamesSymbol));
 
 UNPROTECT(1);
 return(Ans);
}
