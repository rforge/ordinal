#include "distributions_c.h"

/* This file implements scalar distribution, density and gradient
   function */

/*-------------------------------------------------------*/
/* Scalar cumulative distribution functions (CDFs) */
/*-------------------------------------------------------*/

double d_pAO(double q, double lambda, int lower_tail)
{
  if(ISNAN(q) || ISNAN(lambda)) // true for NA and NaN
    return NA_REAL;
  if(q == R_PosInf)
    q = 1;
  else if(q == R_NegInf)
    q = 0;
  else {
    if(lambda < 1.0e-6)
      error("'lambda' has to be positive. lambda = %e was supplied\n",
	    lambda);
    q = 1 - R_pow(lambda * exp(q) + 1, -1/lambda);
  }
  return !lower_tail ? 1 - q : q;
}

double d_plgamma(double eta, double lambda, int lower_tail)
{
  double v;
  if(ISNAN(eta) || ISNAN(lambda)) // true for NA and NaN
    return NA_REAL;
  if(eta == R_PosInf)
    v = 1;
  else if(eta == R_NegInf)
    v = 0;
  else {
    v = R_pow_di(lambda, -2) * exp(lambda * eta);
    if(lambda < 1.0e-6)
      v = 1 - pgamma(v, R_pow_di(lambda, -2), /*scale = */ 1,
		     lower_tail, 0 /*give_log*/);
    if(lambda > -1.0e-6)
      v = pgamma(v, R_pow_di(lambda, -2), /*scale = */ 1,
		 lower_tail, 0/*give_log*/);
    if(lambda >= -1.0e-6 && lambda <= 1.0e-6)
      // pnorm(x, mu, sigma, lower_tail, give_log);
      v = pnorm(eta, 0., 1., 1, 0);
  }
  return !lower_tail ? 1 - v : v;
}

/*-------------------------------------------------------*/
/* Scalar probability density functions (PDFs) */
/*-------------------------------------------------------*/

double d_dAO(double eta, double lambda, int give_log)
{
  if(ISNAN(eta) || ISNAN(lambda)) // true for NA and NaN
    return NA_REAL;
  if(eta == R_PosInf || eta == R_NegInf)
    return 0;
  if(lambda < 1.0e-6)
    error("'lambda' has to be positive. lambda = %e was supplied\n",
	  lambda);
  eta -= (1 + 1 / lambda) * log(lambda * exp(eta) + 1);
  return give_log ? eta : exp(eta);
}

double d_dlgamma(double x, double lambda, int give_log)
{
  if(ISNAN(x) || ISNAN(lambda)) // true for NA and NaN
    return NA_REAL;
  if(x == R_PosInf || x == R_NegInf)
    return 0;
  if(lambda < 1.0e-5 && lambda > -1.0e-5) // lambda close to zero
    return dnorm(x, 0. , 1., give_log);

  double q_2 = R_pow_di(lambda, -2);
  x *= lambda;
  x = log(fabs(lambda)) + q_2 * log(q_2) -
    lgammafn(q_2) + q_2 * (x - exp(x));
  return !give_log ? exp(x) : x;
}

/*-------------------------------------------------------*/
/* Scalar gradients of probability density functions */
/*-------------------------------------------------------*/

void v_glogis(double *x)
{
  if(*x == INFINITY || *x == -INFINITY)
    *x = 0;
  else if(!ISNAN(*x)) {// NA and NaN handled correctly
    /* Store the sign of x, compute the gradient for the absolute value
       and restore the sign. This is needed to avoid exp(LARGE) to blow
       up and the function to return NaN.
    */
    int sign = *x > 0; //could use fsign() instead...
    *x = exp(-fabs(*x));
    *x = 2 * *x * *x * R_pow_di(1 + *x, -3) - *x *
      R_pow_di(1 + *x, -2);
    *x = sign ? *x : -*x;
  }
}

void v_gnorm(double *x, double *mean, double *sd)
{
  if(*x == INFINITY || *x == -INFINITY)
    *x = 0;
  else if(!ISNAN(*x)) // NA and NaN handled correctly
    *x = -(*x - *mean) / *sd * dnorm(*x, *mean, *sd, 0);
}

void v_gcauchy(double *x)
{
  if(*x == R_PosInf || *x == R_NegInf)
    *x = 0;
  else if(!ISNAN(*x))
    *x = -2 * *x / M_PI * R_pow_di(1 + *x * *x, -2);
}

void v_gAO(double *eta, double *lambda)
{
  if(*eta == R_PosInf || *eta == R_NegInf)
    *eta = 0.;
  else if(!ISNAN(*eta)) {
    double lex = *lambda * exp(*eta);
    if(lex == R_PosInf || lex == 0)
      *eta = 0.;
    else {
      double y = d_dAO(*eta, *lambda, 0/*give_log*/);
      *eta = y == 0. ? 0. : y * (1 - (1 + 1 / *lambda) * lex / (1 + lex));
    }
  }
}

void v_glgamma(double *x, double *lambda)
{
  if(*x == R_PosInf || *x == R_NegInf)
    *x = 0.;
  else if(!ISNAN(*x)) {
    if(*lambda < 1.0e-5 && *lambda > -1.0e-5) { // lambda close to zero
      *x = -*x * dnorm(*x, 0., 1., 0/*give_log*/);
    }
    else {
      double z = exp(*lambda * *x);
      if(z == R_PosInf || z == 0.) {
	*x = 0.;
      }
      else {
	double y = d_dlgamma(*x, *lambda, 0/*give_log*/);
	if(y <= 0.)
	  *x = 0.0;
	else
	  *x = y * (1 - exp(*lambda * *x)) / *lambda;
      }
    }
  }
}

double d_glogis(double x)
{
  double y = x;
  v_glogis(&y);
  return y;
}

double d_gnorm(double x, double mean, double sd)
{
  double y = x;
  v_gnorm(&y, &mean, &sd);
  return y;
}

double d_gcauchy(double x)
{
  double y = x;
  v_gcauchy(&y);
  return y;
}

double d_gAO(double x, double lambda)
{
  double y = x;
  v_gAO(&y, &lambda);
  return y;
}

double d_glgamma(double x, double lambda)
{
  double y = x;
  v_glgamma(&y, &lambda);
  return y;
}
