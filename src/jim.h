SEXP C_JIM(SEXP X,SEXP Y,SEXP K,SEXP Threads){
 int n,k,m,ny,*y,*nx,**x,nt;
 struct ht **hta;
 prepareInput(X,Y,K,Threads,&hta,&n,&m,&k,&y,&ny,&x,&nx,&nt);

 double bs=0.,off; int *ctmp,bi=0;
 initialImScan(hta,n,m,y,ny,x,nx,&ctmp,NULL,&bs,&bi,nt,&off);
 if(bs==0) return(makeAns(0,NULL,NULL));

 //Save selected X as W and discard from further consideration
 int* w=x[bi],nw=nx[bi]; x[bi]=NULL;

 //Yet put it as a first selected attribute
 double *score; int *idx;
 SEXP Ans; PROTECT(Ans=makeAns(k,&score,&idx));
 score[0]=bs; idx[0]=bi+1;
 
 //Time for an actual algorithm
 double *as=(double*)R_alloc(sizeof(double),m); //Accumulated score
 for(int e=0;e<m;e++) as[e]=0.;
 int *wxc=(int*)R_alloc(sizeof(int),n*nt),*cWXc=ctmp;
 bs=0.;

 #pragma omp parallel num_threads(nt)
 for(int e=1;e<k;e++){
  double tbs=0.;
  int tbi=-1,tn=omp_get_thread_num();
  struct ht *ht=hta[tn];
  int *wx=wxc+(tn*n),*cWX=cWXc+(tn*n);
  #pragma omp for
  for(int ee=0;ee<m;ee++){
   //Ignore attributes already selected
   if(!x[ee]) continue;

   //Mix x[ee] with lx making wx
   int nwx=fillHt(ht,n,nx[ee],x[ee],nw,w,wx,NULL,NULL,1);

   //Make IM of mix and Y and increase its accumulated score
   fillHt(ht,n,nwx,wx,ny,y,NULL,cWX,NULL,0);
   as[ee]+=imHt(ht,cWX)-off;

   if(as[ee]>tbs){
    tbs=as[ee]; tbi=ee;
   }
  }
  #pragma omp critical
  if((tbs>bs) || (tbs==bs && tbi<bi)){
   bs=tbs;
   bi=tbi;
  }
  
  #pragma omp barrier 
  #pragma omp single
  {
   w=x[bi]; nw=nx[bi]; x[bi]=NULL; 
   score[e]=bs; idx[e]=bi+1;
   bs=0.;
  }
 }

 Ans=finishAns(k,Ans,X);

 UNPROTECT(1);
 return(Ans);
}

