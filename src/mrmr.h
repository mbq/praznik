//This is exactly the same as IF, hence the same implementation
SEXP C_MRMR(SEXP X,SEXP Y,SEXP K){
 int n,k,m,ny,*y,*nx,**x;
 struct ht *ht;
 prepareInput(X,Y,K,&ht,&n,&m,&k,&y,&ny,&x,&nx);

 double bs=0.,*rels=(double*)R_alloc(sizeof(double),m);
 int bi=0,*ctmp,*ctmp2;
 initialMiScan(ht,n,m,y,ny,x,nx,&ctmp,&ctmp2,rels,&bs,&bi);
 if(bs==0) return(makeAns(0,NULL,NULL));
 
 //Save selected X as W and discard from further consideration
 int* w=x[bi],nw=nx[bi]; x[bi]=NULL;

 //Yet put it as a first selected attribute
 double *score; int *idx;
 SEXP Ans; PROTECT(Ans=makeAns(k,&score,&idx));
 score[0]=bs; idx[0]=bi+1;
 
 //Time for an actual MRMR
 double *reds=(double*)R_alloc(sizeof(double),m); //Redundancy
 for(int e=0;e<m;e++) reds[e]=0.;

 for(int e=1;e<k;e++){
  bs=-INFINITY;
  int *cW=ctmp,*cX=ctmp2,dw=0;
  for(int ee=0;ee<m;ee++){
   //Ignore attributes already selected
   if(!x[ee]) continue;

   //WY x[ee] with w making the redundancy part
   fillHt(ht,n,nx[ee],x[ee],nw,w,NULL,cX,dw?NULL:cW,0); dw=1;
   reds[ee]+=miHt(ht,cX,cW);
   double sc=rels[ee]-reds[ee]/(double)e;

   if(sc>bs){
    bs=sc; bi=ee;
   }
  }
  w=x[bi]; nw=nx[bi]; x[bi]=NULL;
  score[e]=bs; idx[e]=bi+1;
 }

 UNPROTECT(1);
 return(Ans);
}

