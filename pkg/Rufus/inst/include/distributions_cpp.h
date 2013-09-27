#ifndef _RUFUS_DISTRIBUTIONS_CPP_H
#define _RUFUS_DISTRIBUTIONS_CPP_H
// #ifndef _r2misc_RCPP_HELLO_WORLD_H
// #define _r2misc_RCPP_HELLO_WORLD_H

#include <Rcpp.h>

using namespace Rcpp;

/*
 * note : RcppExport is an alias to `extern "C"` defined by Rcpp.
 *
 * It gives C calling convention to the rcpp_hello_world function so that
 * it can be called from .Call in R. Otherwise, the C++ compiler mangles the
 * name of the function and .Call can't find it.
 *
 * It is only useful to use RcppExport when the function is intended to be called
 * by .Call. See the thread http://thread.gmane.org/gmane.comp.lang.r.rcpp/649/focus=672
 * on Rcpp-devel for a misuse of RcppExport
 */
// RcppExport SEXP rcpp_hello_world() ;

/* Gumbel distribution functions: */

RcppExport SEXP pgumbel(SEXP q, SEXP location, SEXP scale, SEXP lower_tail, SEXP max);
RcppExport SEXP dgumbel(SEXP x, SEXP location, SEXP scale, SEXP give_log, SEXP max);
RcppExport SEXP ggumbel(SEXP x, SEXP max);
void            v_ggumbel(NumericVector);
NumericVector   ggumbel(NumericVector);


/* Additional gradient functions: */

RcppExport SEXP gnorm(SEXP, SEXP, SEXP);
RcppExport SEXP glogis(SEXP);
RcppExport SEXP gcauchy(SEXP);
RcppExport SEXP gAO(SEXP, SEXP);
RcppExport SEXP glgamma(SEXP, SEXP);

void            v_gnorm(NumericVector, double, double);
void            v_glogis(NumericVector);
void            v_gcauchy(NumericVector);
void            v_gAO(NumericVector, double);
void            v_glgamma(NumericVector, double);

NumericVector   gnorm(NumericVector, double, double);
NumericVector   glogis(NumericVector);
NumericVector   gAO(NumericVector, double);
NumericVector   glgamma(NumericVector, double);

#endif
