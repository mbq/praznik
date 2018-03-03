SEXP C_cmi_jmi(SEXP X,SEXP Y,SEXP Z,SEXP Mode){
 int nt=omp_get_max_threads();
 int n,m,ny,*y,nz,*z,*nx,**x;
 struct ht *hta[nt];
 prepareInput(X,Y,R_NilValue,hta,&n,&m,NULL,&y,&ny,&x,&nx,nt);
 z=convertSEXP(*hta,n,Z,&nz);
 int *cXZc=(int*)R_alloc(sizeof(int),n*nt),
  *cY=(int*)R_alloc(sizeof(int),n),
  *xzc=(int*)R_alloc(sizeof(int),n*nt);

 double miYZ;
 {
  //Calculate I(Y;Z), which is a constant factor for CMI
  int *cZ=cXZc;
  fillHt(*hta,n,ny,y,nz,z,NULL,cY,cZ,0);
  miYZ=miHt(*hta,cY,cZ);
 }
 
 if(length(Mode)!=1) error("Invalid mode");
 int mode=INTEGER(Mode)[0];
 
 SEXP Ans=PROTECT(allocVector(REALSXP,m));
 double *score=REAL(Ans);

 #pragma omp parallel
 {
  int tn=omp_get_thread_num(),*cXZ=cXZc+(tn*n),*xz=xzc+(tn*n),nzx;
  struct ht *ht=hta[tn];
  #pragma omp for
  for(int e=0;e<m;e++){
   //Mix X and Z
   int nxz=fillHt(ht,n,nz,z,nx[e],x[e],xz,NULL,NULL,1);
   fillHt(ht,n,ny,y,nxz,xz,NULL,NULL,cXZ,0);
   // I(X;Y|Z)=I(Y;X,Z)-I(Y;Z)
   score[e]=miHt(ht,cY,cXZ)-miYZ;
  }
 }
 //Copy attribute names
 setAttrib(Ans,R_NamesSymbol,getAttrib(X,R_NamesSymbol));
 
 UNPROTECT(1);
 return(Ans);
}
