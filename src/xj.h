//DISR, JMI, JMIM and NJMIM are basically the same method, hence a common function

enum xjm {xjmAccumulate,xjmMinimum};
SEXP static inline xjcore(SEXP X,SEXP Y,SEXP K,enum xjm mode,double sthHt(struct ht*,int*,int*)){
 int nt=omp_get_max_threads();
 int n,k,m,ny,*y,*nx,**x;
 struct ht *hta[nt];
 prepareInput(X,Y,K,hta,&n,&m,&k,&y,&ny,&x,&nx,nt);

 double bs=0.; int *cY,*ctmp,bi=0;
 initialMiScan(hta,n,m,y,ny,x,nx,&cY,&ctmp,NULL,&bs,&bi,nt);
 if(bs==0) return(makeAns(0,NULL,NULL));

 //Save selected X as W and discard from further consideration
 int* w=x[bi],nw=nx[bi]; x[bi]=NULL;

 //Yet put it as a first selected attribute
 double *score; int *idx;
 SEXP Ans; PROTECT(Ans=makeAns(k,&score,&idx));
 score[0]=bs; idx[0]=bi+1;
 
 //Time for an actual algorithm
 double *as=(double*)R_alloc(sizeof(double),m); //Accumulated score
 for(int e=0;e<m;e++) as[e]=(mode==xjmAccumulate)?0.:INFINITY;
 int *wxc=(int*)R_alloc(sizeof(int),n*nt),*cWXc=ctmp;
 for(int e=1;e<k;e++){
  bs=0.;  
  #pragma omp parallel
  {
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

    //Make MI of mix and Y and increase its accumulated score
    fillHt(ht,n,ny,y,nwx,wx,NULL,NULL,cWX,0); //cY stuff is red.
    if(mode==xjmAccumulate){
     as[ee]+=sthHt(ht,cY,cWX);
    }else{
     double ns=sthHt(ht,cY,cWX);
     as[ee]=(ns<as[ee])?ns:as[ee];
    }

    if(as[ee]>tbs){
     tbs=as[ee]; tbi=ee;
    }
   }
   //Find max among all threads
   #pragma omp critical
   if(tbs>bs){
    bs=tbs;
    bi=tbi;
   }
  }
  //Commit, in a single thread now
  w=x[bi]; nw=nx[bi]; x[bi]=NULL; 
  score[e]=bs; idx[e]=bi+1;
 }

 UNPROTECT(1);
 return(Ans);
}

SEXP C_JMI(SEXP X,SEXP Y,SEXP K){
 return(xjcore(X,Y,K,xjmAccumulate,miHt));
}

SEXP C_DISR(SEXP X,SEXP Y,SEXP K){
 return(xjcore(X,Y,K,xjmAccumulate,nmiHt));
}

SEXP C_JMIM(SEXP X,SEXP Y,SEXP K){
 return(xjcore(X,Y,K,xjmMinimum,miHt));
}

SEXP C_NJMIM(SEXP X,SEXP Y,SEXP K){
 return(xjcore(X,Y,K,xjmMinimum,nmiHt));
}

