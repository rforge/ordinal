#include "distributions_c.h"
// #include "distributions_cpp.h"
#include <R_ext/Rdynload.h>


// static const R_ExternalMethodDef ExtEntries[] = {
//     EXTDEF(Mmatrix, 7),
//     {NULL, NULL, 0}
// };

void
#ifdef HAVE_VISIBILITY_ATTRIBUTE
__attribute__ ((visibility ("default")))
#endif
R_init_Rufus(DllInfo *dll)
{
    // R_registerRoutines(dll, NULL, CallEntries, NULL, ExtEntries);
    // R_useDynamicSymbols(dll, FALSE);

/* These are callable from other packages' C code: */

#define RREGDEF(name)  R_RegisterCCallable("Rufus", #name, (DL_FUNC) name)

    RREGDEF(d_pgumbel);
    RREGDEF(d_pAO);
    RREGDEF(d_plgamma);

    RREGDEF(d_dgumbel);
    RREGDEF(d_dAO);
    RREGDEF(d_dlgamma);

    RREGDEF(d_glogis);
    RREGDEF(d_gnorm);
    RREGDEF(d_gcauchy);
    RREGDEF(d_ggumbel);
    RREGDEF(d_gAO);
    RREGDEF(d_glgamma);
}

