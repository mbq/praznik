SEXP C_im(SEXP X,SEXP Y,SEXP Threads){
 int n,m,ny,*y,*nx,**x,nt;
 struct ht **hta;
 prepareInput(X,Y,R_NilValue,Threads,&hta,&n,&m,NULL,&y,&ny,&x,&nx,&nt);
 int *cXc=(int*)R_alloc(sizeof(int),n*nt);
 SEXP Ans=PROTECT(allocVector(REALSXP,m));
 double *score=REAL(Ans);
 double offset=imOffRaw(*hta,n,y,ny);

 #pragma omp parallel num_threads(nt)
 {
  int tn=omp_get_thread_num(),*cX=cXc+(tn*n);
  struct ht *ht=hta[tn];
  #pragma omp for
  for(int e=0;e<m;e++){
   fillHt(ht,n,nx[e],x[e],ny,y,NULL,cX,NULL,0);
   score[e]=imHt(ht,cX)-offset;
  }
 }
 //Copy attribute names
 setAttrib(Ans,R_NamesSymbol,getAttrib(X,R_NamesSymbol));
 
 UNPROTECT(1);
 return(Ans);
}
