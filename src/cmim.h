//This is exactly the same as IF, hence the same implementation
//TODO: Use Fleuret's trick to speed it up
SEXP C_CMIM(SEXP X,SEXP Y,SEXP K){
 int n,k,m,ny,*y,*nx,**x;
 struct ht *ht;
 prepareInput(X,Y,K,&ht,&n,&m,&k,&y,&ny,&x,&nx);

 double bs=0.; int bi=0,*cY,*ctmp;
 initialMiScan(ht,n,m,y,ny,x,nx,&cY,&ctmp,NULL,&bs,&bi);
 if(bs==0) return(makeAns(0,NULL,NULL));
 
 //Save selected X as W and discard from further consideration
 int* w=x[bi],nw=nx[bi]; x[bi]=NULL;

 //Yet put it as a first selected attribute
 double *score; int *idx;
 SEXP Ans; PROTECT(Ans=makeAns(k,&score,&idx));
 score[0]=bs; idx[0]=bi+1;
 
 //Time for an actual CMIM
 double *ms=(double*)R_alloc(sizeof(double),m); //Minimal score
 for(int e=0;e<m;e++) ms[e]=INFINITY;
 int *wy=(int*)R_alloc(sizeof(int),2*n),*cX=wy+n,*cWY=cY;

 for(int e=1;e<k;e++){
  bs=-INFINITY;
  int nwy=fillHt(ht,n,nw,w,ny,y,wy,NULL,NULL,1);
  for(int ee=0;ee<nwy;ee++) cWY[ee]=ht->cnt[ee].c;
  int *cW=ctmp,dw=0;
  for(int ee=0;ee<m;ee++){
   //Ignore attributes already selected
   if(!x[ee]) continue;

   //WY x[ee] with wy making first part of cmi
   fillHt(ht,n,nx[ee],x[ee],nwy,wy,NULL,cX,NULL,1);
   double cmi=miHt(ht,cX,cWY);
   //WY x[ee] with w to make second part of CMI
   fillHt(ht,n,nw,w,nx[ee],x[ee],NULL,dw?NULL:cW,NULL,0); dw=1;
   cmi-=miHt(ht,cW,cX);
   ms[ee]=(ms[ee]<cmi)?ms[ee]:cmi;

   if(ms[ee]>bs){
    bs=ms[ee]; bi=ee;
   }
  }
  w=x[bi]; nw=nx[bi]; x[bi]=NULL;
  score[e]=bs; idx[e]=bi+1;
 }

 UNPROTECT(1);
 return(Ans);
}

