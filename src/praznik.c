#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <R_ext/Utils.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h> 
#include <omp.h>

//Hash table

#include "ht.h"

//Shared functions

#include "shared.h"

//Algorithms

#include "cmim.h"
#include "mi.h"
#include "mim.h"
#include "mrmr.h"
#include "xj.h"

//Auxiliary

#include "side.h"

//Registration

#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}
static const R_CallMethodDef R_CallDef[]={
 CALLDEF(C_engineTest,2),
 CALLDEF(C_getMi,2),
 CALLDEF(C_getNmi,2),
 CALLDEF(C_setOmpThreads,1),
 CALLDEF(C_MI,3),
 CALLDEF(C_CMIM,3),
 CALLDEF(C_JMI,3),
 CALLDEF(C_DISR,3),
 CALLDEF(C_JMIM,3),
 CALLDEF(C_NJMIM,3),
 CALLDEF(C_MIM,3),
 CALLDEF(C_MRMR,3),
 {NULL,NULL,0}
};

void attribute_visible R_init_praznik(DllInfo *dll){
 R_registerRoutines(dll,NULL,R_CallDef,NULL,NULL);
 R_useDynamicSymbols(dll,FALSE);
 R_forceSymbols(dll,TRUE);
}
