SEXP C_JMIM(SEXP X,SEXP Y,SEXP K){
 int nt=omp_get_max_threads();
 int n,k,m,ny,*y,*nx,**x;
 struct ht *hta[nt];
 prepareInput(X,Y,K,hta,&n,&m,&k,&y,&ny,&x,&nx,nt);

 double bs=0.; int *cY,*ctmp,bi=0;
 initialMiScan(hta,n,m,y,ny,x,nx,&cY,&ctmp,NULL,&bs,&bi,nt);
 if(bs==0) return(makeAns(0,NULL,NULL));

 //Save selected X as W and discard from further consideration
 int **w=(int**)R_alloc(sizeof(int*),k),
  *nw=(int*)R_alloc(sizeof(int),k),
  *lk=(int*)R_alloc(sizeof(int),m);
 w[0]=x[bi]; nw[0]=nx[bi]; x[bi]=NULL;
 for(int e=0;e<m;e++) lk[e]=0;

 //Yet put it as a first selected attribute
 double *score; int *idx;
 SEXP Ans; PROTECT(Ans=makeAns(k,&score,&idx));
 score[0]=bs; idx[0]=bi+1;
 
 //Time for an actual algorithm
 double *ms=(double*)R_alloc(sizeof(double),m);
 for(int e=0;e<m;e++) ms[e]=INFINITY;
 int *wxc=(int*)R_alloc(sizeof(int),n*nt),*cWXc=ctmp,ke=k;
 bs=-INFINITY;

 #pragma omp parallel
 for(int e=1;e<ke;e++){
  double tbs=-INFINITY;
  int tbi=-1,tn=omp_get_thread_num();
  struct ht *ht=hta[tn];
  int *wx=wxc+(tn*n),*cWX=cWXc+(tn*n);
  #pragma omp for schedule(dynamic)
  for(int ee=0;ee<m;ee++){
   //Ignore attributes already selected
   if(!x[ee] || tbs>ms[ee]) continue;

   //Push forward towards e
   for(;lk[ee]<e;lk[ee]++){
    int ew=lk[ee];
    //Mix x[ee] with lx making wx
    int nwx=fillHt(ht,n,nx[ee],x[ee],nw[ew],w[ew],wx,NULL,NULL,1);
    //Make MI of mix and Y and increase its accumulated score
    fillHt(ht,n,ny,y,nwx,wx,NULL,NULL,cWX,0); //cY stuff is red.
    double ns=miHt(ht,cY,cWX);
    ms[ee]=(ns<ms[ee])?ns:ms[ee];
    //Maybe it is already eliminated?
    if(tbs>ms[ee]) break;
   }
   //Check again
   if(ms[ee]>tbs){
    tbs=ms[ee]; tbi=ee;
   }
  }
  #pragma omp critical
  if(tbs>bs){
   bs=tbs;
   bi=tbi;
  }
  #pragma omp barrier 
  #pragma omp single
  if(bs>0.){
   w[e]=x[bi]; nw[e]=nx[bi]; x[bi]=NULL; 
   score[e]=bs; idx[e]=bi+1;
   bs=-INFINITY;
  }else ke=e;
 }

 Ans=finishAns(ke,Ans,X);
 UNPROTECT(1);
 return(Ans);
}

