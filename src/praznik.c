#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <R_ext/Utils.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h> 

#ifdef _OPENMP
 #include <omp.h>
#else
 #define omp_get_thread_num() 0
 #define omp_get_max_threads() 1
 #define omp_set_num_threads(x)
#endif

//Hash table

#include "ht.h"

//Shared functions

#include "shared.h"

//Feature selection algorithms

#include "cmim.h"
#include "mim.h"
#include "mrmr.h"
#include "disr.h"
#include "jmi.h"
#include "jmim.h"
#include "njmim.h"

//Feature scoring algorithms

#include "mi.h"
#include "cmi.h"

//Auxiliary

#include "side.h"

//Registration

#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}
static const R_CallMethodDef R_CallDef[]={
 CALLDEF(C_engineTest,2),
 CALLDEF(C_getMi,2),
 CALLDEF(C_getNmi,2),
 CALLDEF(C_CMIM,4),
 CALLDEF(C_JMI,4),
 CALLDEF(C_DISR,4),
 CALLDEF(C_JMIM,4),
 CALLDEF(C_NJMIM,4),
 CALLDEF(C_MIM,4),
 CALLDEF(C_MRMR,4),
 CALLDEF(C_mi,3),
 CALLDEF(C_cmi_jmi,5),
 {NULL,NULL,0}
};

void attribute_visible R_init_praznik(DllInfo *dll){
 R_registerRoutines(dll,NULL,R_CallDef,NULL,NULL);
 R_useDynamicSymbols(dll,FALSE);
 R_forceSymbols(dll,TRUE);
}
