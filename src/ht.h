#include <stdint.h>

#define GET_A(x) ((x)>>32)
#define GET_B(x) ((x)&0x00000000ffffffff)
#define MAKE_AB(a,b) (((uint64_t)(a)<<32)|(b))

struct ht{
 struct hte **map;
 struct hte *cnt;
 int N;
 uint32_t nAB;
};

struct hte{
 uint64_t ab;
 struct hte *nxt;
 int c;
};

struct ht* R_allocHt(int N){
 struct ht *ans=(struct ht*)R_alloc(sizeof(struct ht),1);
 ans->N=N;
 ans->map=(struct hte**)R_alloc(sizeof(struct hte*),N);
 ans->cnt=(struct hte*)R_alloc(sizeof(struct hte),N);
 return(ans);
}

uint32_t static inline fillHt(struct ht *Q,int N,int nA,int *a,int nB,int *b,int *mix,int *cA,int *cB,int mixOff){
 //Investigate special cases when ht is not needed, i.e. nA/nB=1, nA/nB=N. nA/nB=0 /probably same as =1/ ~> probably not for this function
 if(!b){
  b=a;
  nB=nA;
 }
 //Zero HT
 uint32_t nAB=0;
 for(int e=0;e<N;e++) Q->map[e]=NULL;
 if(cA) for(int e=0;e<nA;e++) cA[e]=0;
 if(cB) for(int e=0;e<nB;e++) cB[e]=0;

 //Fill HT
 for(int e=0;e<N;e++){
  uint32_t _a=a[e]-1;
  uint32_t _b=b[e]-1;
  uint64_t _ab=MAKE_AB(_a,_b);
  uint32_t hab=(nA*nB<N)?(_a+_b*nA):((_a^_b)%N);

  struct hte **E=Q->map+hab;
  for(;E[0] && E[0]->ab!=_ab;E=&(E[0]->nxt));

  if(!E[0]){
   //Empty cell found
   Q->cnt[nAB].ab=_ab;
   Q->cnt[nAB].nxt=NULL;
   Q->cnt[nAB].c=0;
   E[0]=Q->cnt+nAB;
   nAB++;
  }

  //Count this combination
  E[0]->c++;
  if(cA) cA[_a]++;
  if(cB) cB[_b]++;

  if(mix) mix[e]=(E[0]-Q->cnt)+mixOff;
 }
 Q->nAB=nAB;

 return(nAB);
}

uint32_t static inline fillHtOne(struct ht *Q,int N,int *in,int *out,int mixOff){
 //Zero HT
 uint32_t nAB=0;
 for(int e=0;e<N;e++) Q->map[e]=NULL;
 
 //Fill HT
 for(int e=0;e<N;e++){
  if(in[e]==NA_INTEGER) error("NA values are not allowed");
  uint64_t _ab=(uint64_t)(in[e]);
  uint32_t hab=_ab%N;//TOOD: Better hash?

  struct hte **E=Q->map+hab;
  for(;E[0] && E[0]->ab!=_ab;E=&(E[0]->nxt));

  if(!E[0]){
   //Empty cell found
   Q->cnt[nAB].ab=_ab;
   Q->cnt[nAB].nxt=NULL;
   E[0]=Q->cnt+nAB;
   nAB++;
  }
  out[e]=(E[0]-Q->cnt)+mixOff;
 }
 return(nAB);
}


double miHt(struct ht *Q,int *cA,int *cB){
 double ans=0.;
 double N=Q->N;
 for(int e=0;e<Q->nAB;e++){
  if(!(Q->cnt[e].c)) continue;
  double cAB=Q->cnt[e].c;
  double _cA=cA[GET_A(Q->cnt[e].ab)];
  double _cB=cB[GET_B(Q->cnt[e].ab)];
  ans+=cAB*log((cAB*N)/(_cA*_cB));
 }
 return(ans/N);
}

double nmiHt(struct ht *Q,int *cA,int *cB){
 double I=0.;
 double H=0.;
 double N=Q->N;
 for(int e=0;e<Q->nAB;e++){
  if(!(Q->cnt[e].c)) continue;
  double cAB=Q->cnt[e].c;
  double _cA=cA[GET_A(Q->cnt[e].ab)];
  double _cB=cB[GET_B(Q->cnt[e].ab)];
  I+=cAB*log((cAB*N)/(_cA*_cB));
  H+=-cAB*log(cAB/N);
 }
 return(I/H);
}

