SEXP C_engineTest(SEXP A,SEXP B){
 int N;
 N=length(A);
 if(N!=length(B)) error("A and B size mismatch!");

 struct ht *Q=R_allocHt(N);

 int nA=length(getAttrib(A,R_LevelsSymbol));
 int nB=length(getAttrib(B,R_LevelsSymbol));

 if(nA>N) error("A has  more levels then its length; fix that!");
 if(nB>N) error("B has more levels then its length; fix that!");
 int *a=INTEGER(A),*b=INTEGER(B);

 SEXP Ans; PROTECT(Ans=allocVector(VECSXP,5));

 SEXP CountA; PROTECT(CountA=allocVector(INTSXP,nA));
 SEXP CountB; PROTECT(CountB=allocVector(INTSXP,nB));
 SEXP AB; PROTECT(AB=allocVector(INTSXP,N));

 int *cA=INTEGER(CountA),*cB=INTEGER(CountB),*ab=INTEGER(AB);

 int nAB=fillHt(Q,N,nA,a,nB,b,ab,cA,cB,1);

 SEXP CountAB; PROTECT(CountAB=allocVector(INTSXP,nAB));
 int *cAB=INTEGER(CountAB);
 for(int e=0;e<nAB;e++) cAB[e]=Q->cnt[e].c;

 SEXP MI; PROTECT(MI=allocVector(REALSXP,1));
 REAL(MI)[0]=miHt(Q,cA,cB);

 SET_VECTOR_ELT(Ans,0,AB);
 SET_VECTOR_ELT(Ans,1,CountA);
 SET_VECTOR_ELT(Ans,2,CountB);
 SET_VECTOR_ELT(Ans,3,CountAB);
 SET_VECTOR_ELT(Ans,4,MI);
 UNPROTECT(6);
 return(Ans);
}

SEXP C_getMi(SEXP A,SEXP B){
 int N;
 N=length(A);
 if(N!=length(B)) error("A and B size mismatch!");

 struct ht *Q=R_allocHt(N);

 int nA=length(getAttrib(A,R_LevelsSymbol));
 int nB=length(getAttrib(B,R_LevelsSymbol));
 if((nA==0)||(nB==0)) error("A and B have to be factors!");
 if(nA>N) error("A has  more levels then its length; fix that!");
 if(nB>N) error("B has more levels then its length; fix that!");

 int *a=INTEGER(A),*b=INTEGER(B),*cA=(int*)R_alloc(sizeof(int),N*2),*cB=cA+N;
 fillHt(Q,N,nA,a,nB,b,NULL,cA,cB,1);

 SEXP Ans; PROTECT(Ans=allocVector(REALSXP,1));
 REAL(Ans)[0]=miHt(Q,cA,cB);
 
 UNPROTECT(1);
 return(Ans);
}

SEXP C_getNmi(SEXP A,SEXP B){
 int N;
 N=length(A);
 if(N!=length(B)) error("A and B size mismatch!");

 struct ht *Q=R_allocHt(N);

 int nA=length(getAttrib(A,R_LevelsSymbol));
 int nB=length(getAttrib(B,R_LevelsSymbol));
 if((nA==0)||(nB==0)) error("A and B have to be factors!");
 if(nA>N) error("A has  more levels then its length; fix that!");
 if(nB>N) error("B has more levels then its length; fix that!");

 int *a=INTEGER(A),*b=INTEGER(B),*cA=(int*)R_alloc(sizeof(int),N*2),*cB=cA+N;
 fillHt(Q,N,nA,a,nB,b,NULL,cA,cB,1);

 SEXP Ans; PROTECT(Ans=allocVector(REALSXP,1));
 REAL(Ans)[0]=nmiHt(Q,cA,cB);
 
 UNPROTECT(1);
 return(Ans);
}

SEXP C_setOmpThreads(SEXP threads){
 omp_set_num_threads(INTEGER(threads)[0]);
 return(R_NilValue);
}
