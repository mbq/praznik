int *convertSEXP(struct ht *ht,int n,SEXP in,int *nout){
 int lc=length(getAttrib(in,R_LevelsSymbol)),*out;
 if(isFactor(in) && lc<n){
  //Well-behaving factor; user shall give us such for optimal performance
  *nout=lc;
  out=INTEGER(in);
  for(int e=0;e<n;e++)
   if(out[e]==NA_INTEGER) error("NA values are not allowed");
  return(out);
 }
 if(isFactor(in)||isLogical(in)||isInteger(in)){
  //Integer-alike which needs collapsing into 1..n_levels
  int *out=(int*)R_alloc(sizeof(int),n);
  *nout=fillHtOne(ht,n,INTEGER(in),out,1);
  return(out);
 }
 if(isReal(in)){
  //Magically make discrete by scattering into 10-bins
  double *x=REAL(in),min=INFINITY,max=-INFINITY;
  for(int e=0;e<n;e++){
   if(!R_FINITE(x[e])) error("Non-finite numeric values are not allowed");
   min=min<x[e]?min:x[e];
   max=max>x[e]?max:x[e];
  }
  int *out=(int*)R_alloc(sizeof(int),n);
  if(n<6){
   *nout=2;
  }else if(n>30){
   *nout=10;
  }else{
   *nout=n/3;
  }
  for(int e=0;e<n;e++)
   out[e]=((int)((x[e]-min)/(max-min)*(double)(*nout)))%(*nout)+1;
  return(out);
 }
 //Other stuff
 return(NULL);
}

void prepareInput(SEXP X,SEXP Y,SEXP K,struct ht **ht,int *n,int *m,int *k,int **y,int *ny,int ***x,int **nx){
 if(!isFrame(X)) error("X must be a data.frame");
 *n=length(Y);
 *k=INTEGER(K)[0];
 *m=length(X);

 if(*k>*m) error("Parameter k must be at most the number of attributes");
 if(*k<1) error("Parameter k must be positive");
 if(n[0]>2147483648) error("Only at most 2^31 (2.1 billion) objects allowed");
 //TODO: Also eat matrices? --> then fix it
 if(*n!=length(VECTOR_ELT(X,0))) error("X and Y size mismatch");

 *ht=R_allocHt(*n);

  *y=convertSEXP(*ht,*n,Y,ny);
  if(!*y) error("Wrong Y type");
 
 *nx=(int*)R_alloc(sizeof(int),*m);
 *x=(int**)R_alloc(sizeof(int*),*m);
 for(int e=0;e<*m;e++){
  SEXP XX;
  PROTECT(XX=VECTOR_ELT(X,e));
  (*x)[e]=convertSEXP(*ht,*n,XX,(*nx)+e);
  if(!(*x)[e]) error("Wrong X[,%d] type",e+1);
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
#define EPS 1e-10
