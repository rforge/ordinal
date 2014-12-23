#ifndef _ORDINAL_LINKS_H_
#define _ORDINAL_LINKS_H_
/* That ifndef, etc. is an idiom to prevent the body of the header
 * being read more than once.
 */

#include <R.h>
#include <Rmath.h>

#ifdef	__cplusplus
extern "C" {
#endif
  /* That stanza allows the same header file to be used by C and C++
   * programs. There is a matching stanza at the end of this header
   * file.
   */
  
  /* Additional scalar cumulative probability functions */
  double d_pgumbel  (double,double,double,int);
  double d_pgumbel2 (double,double,double,int);
  double d_pAO      (double,double,int);
  double d_plgamma  (double,double,int);
  
  /* Additional scalar density functions */
  double d_dgumbel  (double,double,double,int);
  double d_dgumbel2 (double,double,double,int);
  double d_dAO      (double,double,int);
  double d_dlgamma  (double,double,int);
  
  /* Scalar density gradients */
  double d_glogis   (double);
  double d_gnorm    (double);
  double d_gcauchy  (double);
  double d_ggumbel  (double);
  double d_ggumbel2 (double);
  double d_gAO      (double,double);
  double d_glgamma  (double,double);
  
#ifdef	__cplusplus
}
#endif

#endif
