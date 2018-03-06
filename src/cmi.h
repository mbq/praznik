enum cmi_jmi_mode {cjmCMI=791,cjmJMI=792,cjmNJMI=793};

SEXP C_cmi_jmi(SEXP X,SEXP Y,SEXP Z,SEXP Mode,SEXP Threads){
 if(length(Mode)!=1) error("Invalid mode");
 int mode=INTEGER(Mode)[0];
 if(mode!=cjmCMI && mode!=cjmJMI && mode!=cjmNJMI)
  error("Invalid mode"); 

 int n,m,ny,*y,nz,*z,*nx,**x,nt;
 struct ht **hta;
 prepareInput(X,Y,R_NilValue,Threads,&hta,&n,&m,NULL,&y,&ny,&x,&nx,&nt);

 if(length(Z)!=n) error("Z vector size mismatch");
 z=convertSEXP(*hta,n,Z,&nz);

 int *cXZc=(int*)R_alloc(sizeof(int),n*nt),
  *cY=(int*)R_alloc(sizeof(int),n),
  *xzc=(int*)R_alloc(sizeof(int),n*nt);

 double scoreOff;
 {
  //Calculate I(Y;Z), which is a constant factor for CMI
  int *cZ=cXZc;
  //This is always needed for cY
  fillHt(*hta,n,ny,y,nz,z,NULL,cY,cZ,0);
  scoreOff=(mode==cjmCMI)?-miHt(*hta,cY,cZ):0.;
 }
 
 SEXP Ans=PROTECT(allocVector(REALSXP,m));
 double *score=REAL(Ans);

 #pragma omp parallel num_threads(nt)
 {
  int tn=omp_get_thread_num(),*cXZ=cXZc+(tn*n),*xz=xzc+(tn*n);
  struct ht *ht=hta[tn];
  #pragma omp for
  for(int e=0;e<m;e++){
   //Mix X and Z
   int nxz=fillHt(ht,n,nz,z,nx[e],x[e],xz,NULL,NULL,1);
   fillHt(ht,n,ny,y,nxz,xz,NULL,NULL,cXZ,0);
   // I(X;Y|Z)=I(Y;X,Z)-I(Y;Z)
   if(mode!=cjmNJMI){
    //CMI or JMI=I(Y;X,Z)
    score[e]=miHt(ht,cY,cXZ)+scoreOff;
   }else{
    //NJMI (i.e. DISR-like score)=I(Y;X,Z)/H(X,Y,Z)
    score[e]=nmiHt(ht,cY,cXZ);
   }
  }
 }
 //Copy attribute names
 setAttrib(Ans,R_NamesSymbol,getAttrib(X,R_NamesSymbol));
 
 UNPROTECT(1);
 return(Ans);
}
