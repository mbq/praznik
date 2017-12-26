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
 //Investigate the mode of HT; is it really needed?
 
 if(nA*nB<N){
  //OK, HT not needed, just use it as a table; map is not needed, we just 
  int Neff=nA*nB;

  //Initiate
  uint32_t nAB=0;
  if(cA) for(int e=0;e<nA;e++) cA[e]=0;
  if(cB) for(int e=0;e<nB;e++) cB[e]=0;

  if(!mix){
   //Even deeper shortcut, no map at all
   Q->nAB=Neff;
   for(int e=0;e<Neff;e++) Q->cnt[e].c=0;
   for(int e=0;e<N;e++){
    uint32_t _a=a[e]-1;
    uint32_t _b=b[e]-1;
    uint64_t _ab=MAKE_AB(_a,_b);
    uint32_t hab=_a+_b*nA;
    Q->cnt[hab].c++;
    Q->cnt[hab].ab=_ab;
    if(cA) cA[_a]++;
    if(cB) cB[_b]++;
   }

   return(Neff);
  }
  
  //Zero the map
  for(int e=0;e<Neff;e++) Q->map[e]=NULL;

  //Fill HT
  for(int e=0;e<N;e++){
   uint32_t _a=a[e]-1;
   uint32_t _b=b[e]-1;
   uint64_t _ab=MAKE_AB(_a,_b);
   uint32_t hab=_a+_b*nA;

   struct hte **E=Q->map+hab;

   if(!*E){
    //Empty cell found
    Q->cnt[nAB].ab=_ab;
    Q->cnt[nAB].c=0;
    *E=Q->cnt+nAB;
    nAB++;
   }

   //Count this combination
   (*E)->c++;
   if(cA) cA[_a]++;
   if(cB) cB[_b]++;

   //TODO: Tautology here
   if(mix) mix[e]=(*E-Q->cnt)+mixOff;
  }
  Q->nAB=nAB;

  return(nAB);
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
  uint32_t hab=(_a^_b)%N; //TODO: Real hash?

  struct hte **E=Q->map+hab;
  for(;(*E)&&(*E)->ab!=_ab;E=&((*E)->nxt));

  if(!*E){
   //Empty cell found
   Q->cnt[nAB].ab=_ab;
   Q->cnt[nAB].nxt=NULL;
   Q->cnt[nAB].c=0;
   *E=Q->cnt+nAB;
   nAB++;
  }

  //Count this combination
  (*E)->c++;
  if(cA) cA[_a]++;
  if(cB) cB[_b]++;

  if(mix) mix[e]=(*E-Q->cnt)+mixOff;
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

