void prepareInput(SEXP X,SEXP Y,SEXP K,struct ht **ht,int *n,int *m,int *k,int **y,int *ny,int ***x,int **nx){
 //X must be a factor data frame, Y must be a factor
 *n=length(Y);
 *k=INTEGER(K)[0];
 *m=length(X);

 if(*k>*m) error("Parameter k must be smaller than the number of attributes");
 if(*k<1) error("Parameter k must be positive");
 //TODO: Resolve this:
 //if(n[0]>=(1<<33)) error("Only at most 2^32 (4.2 billion) objects allowed");
 //TODO: Also eat matrices? --> then fix it
 if(*n!=length(VECTOR_ELT(X,0))) error("X and Y size mismatch");

 *ht=R_allocHt(*n);

 if(!isFactor(Y)) error("Y has to be a factor");
 *ny=length(getAttrib(Y,R_LevelsSymbol));
 if(*ny<*n){
  *y=INTEGER(Y);
 }else{
  *y=(int*)R_alloc(sizeof(int),*n);
  *ny=fillHt(*ht,*n,*ny,INTEGER(Y),0,NULL,*y,NULL,NULL,1);
 }

 *nx=(int*)R_alloc(sizeof(int),*m);
 *x=(int**)R_alloc(sizeof(int*),*m);
 for(int e=0;e<*m;e++){
  SEXP XX;
  PROTECT(XX=VECTOR_ELT(X,e));
  //TODO: Auto-cut numerics into 10-bin
  if(!isFactor(XX)) error("All X columns must be factor!");
  (*nx)[e]=length(getAttrib(XX,R_LevelsSymbol));
  if((*nx)[e]<*n){
   (*x)[e]=INTEGER(XX);
  }else{
   (*x)[e]=(int*)R_alloc(sizeof(int),*n);
   (*nx)[e]=fillHt(*ht,*n,(*nx)[e],INTEGER(XX),0,NULL,(*x)[e],NULL,NULL,1);
  }
  UNPROTECT(1);
 }
}

void static inline initialMiScan(struct ht* ht,int n,int m,int *y,int ny,int **x,int *nx,int **_cY,int **_cX,double *_mi,double *bs,int *bi){
 int *cX=(int*)R_alloc(sizeof(int),n); if(_cX) *_cX=cX;
 int *cY=(int*)R_alloc(sizeof(int),n); if(_cY) *_cY=cY;
 *bs=0.;
 for(int e=0;e<m;e++){
  fillHt(ht,n,ny,y,nx[e],x[e],NULL,e?NULL:cY,cX,0);
  double mi=miHt(ht,cY,cX);
  if(_mi) _mi[e]=mi;
  if(mi>*bs){
   *bs=mi;
   *bi=e;
  }
 }
}

SEXP makeAns(int k,double **score,int **idx){
 SEXP Ans; PROTECT(Ans=allocVector(VECSXP,2));
 SEXP Idx; PROTECT(Idx=allocVector(INTSXP,k));
 SEXP Score; PROTECT(Score=allocVector(REALSXP,k));
 SET_VECTOR_ELT(Ans,0,Idx);
 SET_VECTOR_ELT(Ans,1,Score);

 if(score) *score=REAL(Score);
 if(idx) *idx=INTEGER(Idx);
 UNPROTECT(3);
 return(Ans);
}

//Macros

#define ISWAP(x,y) do{int *__tmp=x;x=y;y=__tmp;}while(0)

