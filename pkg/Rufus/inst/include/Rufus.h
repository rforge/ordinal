#ifndef RUFUS_H
#define RUFUS_H
/* That ifndef, etc. is an idiom to prevent the body of the header
 * being read more than once.
 */

#include <Rdefines.h> // Need these?
#include <Rconfig.h> // Need these?

// #include <R.h>
// #include <Rmath.h>


#ifdef	__cplusplus
extern "C" {
#endif
  /* That stanza allows the same header file to be used by C and C++
   * programs. There is a matching stanza at the end of this header
   * file.
   */

#ifdef HAVE_VISIBILITY_ATTRIBUTE
# define attribute_hidden __attribute__ ((visibility ("hidden")))
#else
# define attribute_hidden
#endif

  /* Additional scalar cumulative probability functions */
  double R_d_pgumbel  (double,double,double,int);
  // void   R_v_pgumbel  (double *, double *, double *, int *);
  double R_d_pAO      (double,double,int);
  double R_d_plgamma  (double,double,int);

  /* Additional scalar density functions */
  double R_d_dgumbel  (double,double,double,int);
  // void   R_v_dgumbel  (double *, double *, double *, int *);
  double R_d_dAO      (double,double,int);
  double R_d_dlgamma  (double,double,int);

  /* Scalar density gradients */
  // void   v_glogis   (double *);
  // void   v_gnorm    (double *, double *, double *);
  // void   v_ggumbel  (double *);
  // void   v_gcauchy  (double *);
  // void   v_gAO      (double *, double *);
  // void   v_glgamma  (double *, double *);

  double R_d_glogis   (double);
  double R_d_gnorm    (double, double, double);
  double R_d_ggumbel  (double);
  double R_d_gcauchy  (double);
  double R_d_gAO      (double,double);
  double R_d_glgamma  (double,double);

#ifdef	__cplusplus
}
#endif

#endif
