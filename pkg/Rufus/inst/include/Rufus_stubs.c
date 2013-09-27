#include <Rconfig.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <Rversion.h> // probably don't need this one.
#include "Rufus.h"
// #include "distributions_c.h"
// #include "distributions_cpp.h" // probably don't need this one -- will it hurt?

#ifdef	__cplusplus
extern "C" {
// and  bool is defined
// #else
// # define bool Rboolean
#endif

#ifdef HAVE_VISIBILITY_ATTRIBUTE
# define attribute_hidden __attribute__ ((visibility ("hidden")))
#else
# define attribute_hidden
#endif

  // Exported CDFs p[dist] functions:

double attribute_hidden
R_d_pgumbel(double q, double loc, double scale, int lower_tail)
{
  static double(*fun)(double,double,double,int) = NULL;
  if(fun == NULL)
    fun = (double(*)(double,double,double,int))
      R_GetCCallable("Rufus", "d_pgumbel");
  return fun(q, loc, scale, lower_tail);
}

double attribute_hidden
R_d_pAO(double q, double lambda, int lower_tail)
{
  static double(*fun)(double,double,int) = NULL;
  if(fun == NULL)
    fun = (double(*)(double,double,int))
      R_GetCCallable("Rufus", "d_pAO");
  return fun(q, lambda, lower_tail);
}

double attribute_hidden
R_d_plgamma(double eta, double lambda, int lower_tail)
{
  static double(*fun)(double,double,int) = NULL;
  if(fun == NULL)
    fun = (double(*)(double,double,int))
      R_GetCCallable("Rufus", "d_plgamma");
  return fun(eta, lambda, lower_tail);
}

  // Exported PDFs d[dist] functions:

double attribute_hidden
R_d_dgumbel(double q, double loc, double scale, int give_log)
{
  static double(*fun)(double,double,double,int) = NULL;
  if(fun == NULL)
    fun = (double(*)(double,double,double,int))
      R_GetCCallable("Rufus", "d_dgumbel");
  return fun(q, loc, scale, give_log);
}

double attribute_hidden
R_d_dAO(double eta, double lambda, int give_log)
{
  static double(*fun)(double,double,int) = NULL;
  if(fun == NULL)
    fun = (double(*)(double,double,int))
      R_GetCCallable("Rufus", "d_dAO");
  return fun(eta, lambda, give_log);
}

double attribute_hidden
R_d_dlgamma(double x, double lambda, int give_log)
{
  static double(*fun)(double,double,int) = NULL;
  if(fun == NULL)
    fun = (double(*)(double,double,int))
      R_GetCCallable("Rufus", "d_dlgamma");
  return fun(x, lambda, give_log);
}

  // Exported gradients of PDFs: g[dist] functions:

double attribute_hidden
R_d_glogis(double x)
{
  static double(*fun)(double) = NULL;
  if(fun == NULL)
    fun = (double(*)(double))
      R_GetCCallable("Rufus", "d_glogis");
  return fun(x);
}

double attribute_hidden
R_d_gnorm(double x, double mean, double sd)
{
  static double(*fun)(double,double,double) = NULL;
  if(fun == NULL)
    fun = (double(*)(double,double,double))
      R_GetCCallable("Rufus", "d_gnorm");
  return fun(x, mean, sd);
}

double attribute_hidden
R_d_ggumbel(double x)
{
  static double(*fun)(double) = NULL;
  if(fun == NULL)
    fun = (double(*)(double))
      R_GetCCallable("Rufus", "d_ggumbel");
  return fun(x);
}

double attribute_hidden
R_d_gcauchy(double x)
{
  static double(*fun)(double) = NULL;
  if(fun == NULL)
    fun = (double(*)(double))
      R_GetCCallable("Rufus", "d_gcauchy");
  return fun(x);
}

double attribute_hidden
R_d_gAO(double x, double lambda)
{
  static double(*fun)(double,double) = NULL;
  if(fun == NULL)
    fun = (double(*)(double,double))
      R_GetCCallable("Rufus", "d_gAO");
  return fun(x, lambda);
}

double attribute_hidden
R_d_glgamma(double x, double lambda)
{
  static double(*fun)(double,double) = NULL;
  if(fun == NULL)
    fun = (double(*)(double,double))
      R_GetCCallable("Rufus", "d_glgamma");
  return fun(x, lambda);
}

  /*-------------------------------------------------------*/


#ifdef	__cplusplus
}
#endif

