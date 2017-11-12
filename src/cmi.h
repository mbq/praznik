SEXP C_CMI(SEXP X,SEXP Y,SEXP K){
 int n,k,m,ny,*y,*nx,**x;
 struct ht *ht;
 prepareInput(X,Y,K,&ht,&n,&m,&k,&y,&ny,&x,&nx);

 int *cX=(int*)R_alloc(sizeof(int),3*n),*cYZ=cX+n,*cZ=cYZ+n;

 //We need to remember last, currently calculated and best 
 int nyz,nz,ncyz,ncz,nbyz=0,nbz=0,
  *yz=(int*)R_alloc(sizeof(int),n*6),
  *z=yz+n,*cyz=z+n,*cz=cyz+n,*byz=cz+n,*bz=byz+n;
 //For optimisation, last CMI record is also useful (I(X;Y)>=I(X;Y|Z) for any Z
 double *cmi=(double*)R_alloc(sizeof(double),m);
 for(int e=0;e<m;e++) cmi[e]=INFINITY;

 //For answer (we don't know the length a priori, CMI can become 0 at any time)
 double *_score=(double*)R_alloc(sizeof(double),k);
 int *_idx=(int*)R_alloc(sizeof(int),k),_k=0; 

 //We need to copy yz since it will be modified and y can be a result of INTEGER
 for(int e=0;e<n;e++){
  z[e]=1; yz[e]=y[e];
 }
 nz=1; nyz=ny;
 
 for(int e=0;e<k;e++){
  double bs=0.; int bi=0;
  int dyz=0,dz=0;
  for(int ee=0;ee<m;ee++){
   if(!x[ee]) continue;
   if(cmi[ee]<bs) continue; //Updated CMI can be only lower
   ncyz=fillHt(ht,n,nx[ee],x[ee],nyz,yz,cyz,cX,dyz?NULL:cYZ,1); dyz=1;
   cmi[ee]=miHt(ht,cX,cYZ);
   //Shortcut: miX_Z>0, hence final CMI won't get higher, and Ht takes a bit of time
   if(cmi[ee]<=bs) continue;
   ncz=fillHt(ht,n,nx[ee],x[ee],nz,z,cz,NULL,dz?NULL:cZ,1); dz=1;
   cmi[ee]-=miHt(ht,cX,cZ);

   if(cmi[ee]>bs){
    bs=cmi[ee]; bi=ee;
    nbyz=ncyz; nbz=ncz;
    ISWAP(byz,cyz);
    ISWAP(bz,cz);
   }
  }
  if(bs<EPS) break; //Won't be anymore hits
  _score[_k]=bs; _idx[_k]=bi+1; _k++;
  nyz=nbyz; nz=nbz;
  x[bi]=NULL;
  ISWAP(yz,byz);
  ISWAP(z,bz);
 }

 double *score; int *idx;
 SEXP Ans; PROTECT(Ans=makeAns(_k,&score,&idx));
 for(int e=0;e<_k;e++){
  score[e]=_score[e];
  idx[e]=_idx[e];
 }

 UNPROTECT(1);
 return(Ans);
}

