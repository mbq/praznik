SEXP C_MIM(SEXP X,SEXP Y,SEXP K){
 int nt=omp_get_max_threads();
 int n,k,m,ny,*y,*nx,**x;
 struct ht *hta[nt];
 prepareInput(X,Y,K,hta,&n,&m,&k,&y,&ny,&x,&nx,nt);
 int *cXc=(int*)R_alloc(sizeof(int),n*nt);
 int *cYc=(int*)R_alloc(sizeof(int),n*nt);
 double *mi=(double*)R_alloc(sizeof(double),m);

 double *score;
 int *idx;
 SEXP Ans; PROTECT(Ans=makeAns(k,&score,&idx));

 for(int e=0;e<k;e++){
  score[e]=-INFINITY; idx[e]=0;
 }

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

 //Insertion sort since we need stable sort in principle
 for(int e=0;e<m;e++) if(score[k-1]<=mi[e]){
  int ee;
  for(ee=k-2;ee>=0 && score[ee]<mi[e];ee--){
   score[ee+1]=score[ee]; idx[ee+1]=idx[ee];
  }
  score[ee+1]=mi[e]; idx[ee+1]=e+1;
 }

 UNPROTECT(1);
 return(Ans);
}

