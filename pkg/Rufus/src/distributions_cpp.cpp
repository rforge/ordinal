#include "distributions_cpp.h"
#include "distributions_c.h"

using namespace Rcpp; // rather than repeating it in each function.

// SEXP pgumbel(SEXP q, SEXP location, SEXP scale, SEXP lower_tail, SEXP max)
// {
//   BEGIN_RCPP
//   NumericVector x = clone(q);
//
//   double loc_ = as<double>(location);
//   double scale_= as<double>(scale);
//   int lower_ = as<int>(lower_tail);
//   int max_ = as<int>(max);
//
//   //if(max_ <= 0)
//   if(!max_ )
//     x = -x;
//
//   for(int i = 0; i < x.size(); i++)
//     // v_pgumbel(&q_[i], &loc_, &scale, &lower_tail);
//     x[i] = d_pgumbel(x[i], loc_, scale_, lower_);
//
//   return wrap(x);
//   END_RCPP
// }
//
/* -------------------------------------------------------
   Vectorized void functions:
   ------------------------------------------------------- */

// Overloading the v_gnorm function also defined for double in links.c
void v_gnorm(NumericVector x, double mean, double sd)
// updating x in-place.
{
  for(int i = 0; i < x.size(); i++)
    // x[i] = d_gnorm(x[i], mean, sd);
    v_gnorm(&x[i], &mean, &sd);
}

void v_glogis(NumericVector x)
{
  for(int i = 0; i < x.size(); i++)
    v_glogis(&x[i]);
}

void v_gcauchy(NumericVector x)
{
  for(int i = 0; i < x.size(); i++)
    v_gcauchy(&x[i]);
}

void v_gAO(NumericVector x, double lambda)
{
  for(int i = 0; i < x.size(); i++)
    v_gAO(&x[i], &lambda);
}

void v_glgamma(NumericVector x, double lambda)
{
  for(int i = 0; i < x.size(); i++)
    v_glgamma(&x[i], &lambda);
}

/* -------------------------------------------------------
   Vectorized NumericVector functions:
   ------------------------------------------------------- */

NumericVector gnorm(NumericVector x, double mean, double sd)
{
  NumericVector y = clone(x);
  v_gnorm(y, mean, sd);
  // for(int i = 0; i < y.size(); i++)
  //   v_gnorm(&y[i], &mean, &sd);
  return y;
}

NumericVector glogis(NumericVector x)
{
  NumericVector y = clone(x);
  for(int i = 0; i < y.size(); i++)
    v_glogis(&y[i]);
  return y;
}

NumericVector gcauchy(NumericVector x)
{
  NumericVector y = clone(x);
  for(int i = 0; i < y.size(); i++)
    v_gcauchy(&y[i]);
  return y;
}

NumericVector gAO(NumericVector x, double lambda)
{
  NumericVector y = clone(x);
  v_gAO(y, lambda);
  return y;
}

NumericVector glgamma(NumericVector x, double lambda)
{
  NumericVector y = clone(x);
  v_glgamma(y, lambda);
  return y;
}

/* -------------------------------------------------------
   .Call interface functions:
   ------------------------------------------------------- */

// callable from R using .Call("gnorm", mean, sd)
SEXP gnorm(SEXP x, SEXP mean, SEXP sd)
{
  BEGIN_RCPP // Handle expections correctly
    // clone()ing x so as not to modify x n return:
    NumericVector z = clone(x);
  // Checking the right lenght of arguments:
  if(LENGTH(mean) >= 2)
    Rcpp::stop("mean has to be a numeric scalar");
  if(LENGTH(sd) >= 2)
    Rcpp::stop("length(sd) != 1");
  // coarsing SEXPs to double:
  double mean_ = as<double>(mean);
  double sd_ = as<double>(sd);

  // z = gnorm(z, mean_, sd_);
  // return wrap(z);
  // return wrap(gnorm(z, mean_, sd_));

  // Update z in-place with the values of the gradient:
  v_gnorm(z, mean_, sd_);

  // Updating z in-place with gradient values calling C-function
  // directly:
  // for(int i = 0; i < z.size(); i++)
  //   v_gnorm(&z[i], &mean_, &sd_);

  // wrap()ing result to the appropriate SEXP type on return to R:
  return wrap(z);
  END_RCPP
}

SEXP glogis(SEXP x)
{
  BEGIN_RCPP
  NumericVector z = clone(x);
  v_glogis(z);
  return wrap(z);
  END_RCPP
}

SEXP gcauchy(SEXP x)
{
  BEGIN_RCPP
  NumericVector z = clone(x);
  v_gcauchy(z);
  return wrap(z);
  END_RCPP
}

SEXP gAO(SEXP x, SEXP lambda)
{
  BEGIN_RCPP
  NumericVector z = clone(x);

  if(LENGTH(lambda) >= 2)
    Rcpp::stop("lambda has to be a numeric scalar");
  double lambda_ = as<double>(lambda);

  v_gAO(z, lambda_);
  return wrap(z);
  END_RCPP
}

SEXP glgamma(SEXP x, SEXP lambda)
{
  BEGIN_RCPP
  NumericVector z = clone(x);

  if(LENGTH(lambda) >= 2)
    Rcpp::stop("lambda has to be a numeric scalar");
  double lambda_ = as<double>(lambda);

  v_glgamma(z, lambda_);
  return wrap(z);
  END_RCPP
}

/* -------------------------------------------------------
   misc
   ------------------------------------------------------- */

SEXP gnorm2(SEXP x, SEXP mean, SEXP sd) {
  try {
    // NumericVector y(x);
    NumericVector z = clone(x);

    if(LENGTH(mean) >= 2)
      Rcpp::stop("mean has to be a numeric scalar");
    if(LENGTH(sd) >= 2)
      Rcpp::stop("length(sd) != 1");

    double l_mean = as<double>(mean);
    double l_sd = as<double>(sd);

    for(int i = 0; i < z.size(); i++) {
      z[i] = d_gnorm(z[i], l_mean, l_sd);
    }

    return wrap(z);

  } catch( std::exception &ex ) {
    forward_exception_to_r( ex );
  } catch(...) {
    ::Rf_error( "c++ exception (unknown reason)" );
  }
  return R_NilValue; // -Wall
}

// Using d_gnorm from links.h copied from the ordinal package:
SEXP gnorm3(SEXP x, SEXP mean, SEXP sd) {
  try {
    // NumericVector y(x);
    NumericVector z = clone(x);

    if(LENGTH(mean) >= 2)
      Rcpp::stop("mean has to be a numeric scalar");
    if(LENGTH(sd) >= 2)
      Rcpp::stop("length(sd) != 1");

    double l_mean = as<double>(mean);
    double l_sd = as<double>(sd);

    z = (z - l_mean) / l_sd;
    for(int i = 0; i < z.size(); i++) {
      z[i] = d_gnorm(z[i], 0., 1.) / l_sd;
    }

    return wrap(z);

  } catch( std::exception &ex ) {
    forward_exception_to_r( ex );
  } catch(...) {
    ::Rf_error( "c++ exception (unknown reason)" );
  }
  return R_NilValue; // -Wall
}

SEXP gnorm4(SEXP x, SEXP mean, SEXP sd) {
  try {
    NumericVector y(x);
    NumericVector z = clone(y);

    if(LENGTH(mean) >= 2)
      Rcpp::stop("mean has to be a numeric scalar");
    if(LENGTH(sd) >= 2)
      Rcpp::stop("length(sd) != 1");

    double l_mean = as<double>(mean);
    double l_sd = as<double>(sd);

    // Standardized variable:
    z = (y - l_mean) / l_sd;
    // Handle Inf and -Inf:
    z = ifelse((y == R_PosInf | y == R_NegInf), 0.,
    	       -z / l_sd * dnorm(z));
    // Handle missing values in y:
    if(is_true(any(is_na(y))))
      z = ifelse(is_na(y), NA_REAL, z);

    return wrap(z);

  } catch( std::exception &ex ) {
    forward_exception_to_r( ex );
  } catch(...) {
    ::Rf_error( "c++ exception (unknown reason)" );
  }
  return R_NilValue; // -Wall
}

// double d_glogis(double x)
// {
//   // Gradient of dlogis(x) wrt. x
//   if(ISNAN(x)) // true for NA and NaN
//     return NA_REAL;
//   if(x == R_PosInf || x == R_NegInf)
//     // if(x == INFINITE || x == -INFINITE) // seems to work as well.
//     return 0; // this special case needs to be handled separately
//
//   /* Store the sign of x, compute the gradient for the absolute value
//      and restore the sign. This is needed to avoid exp(LARGE) to blow
//      up and the function to return NaN.
//   */
//   int sign = x > 0; //could use fsign() instead...
//   x = exp(-fabs(x));
//   x = 2 * x * x * R_pow_di(1 + x, -3) - x *
//     R_pow_di(1 + x, -2);
//   return sign ? x : -x;
// }

// double d_gnorm2(double x, double mean, double sd)
// {
//   if(ISNAN(x)) // true for NA and NaN
//       return NA_REAL;
//   if(x == INFINITY || x == -INFINITY)
//     return 0;
//   else
//     return -(x - mean) / sd * R::dnorm(x, mean, sd, 0);
// }

// NumericVector nv_gnorm(NumericVector x, double mean, double sd)
// {
//   for(int i = 0; i < x.size(); i++) {
//     x[i] = d_gnorm(x[i], mean, sd);
//   }
//   return x;
// }

// SEXP gnorm(SEXP x, SEXP mean, SEXP sd) {
//   try {
//     NumericVector y(x);
//     NumericVector z = clone(y);
//
//     if(LENGTH(mean) >= 2)
//       Rcpp::stop("mean has to be a numeric scalar");
//     if(LENGTH(sd) >= 2)
//       Rcpp::stop("length(sd) != 1");
//
//     double l_mean = as<double>(mean);
//     double l_sd = as<double>(sd);
//
//     // Standardized variable:
//     z = (y - l_mean) / l_sd;
//     // Handle Inf and -Inf:
//     z = ifelse((y == R_PosInf | y == R_NegInf), 0.,
//     	       -z / l_sd * dnorm(z));
//     // Handle missing values in y:
//     if(is_true(any(is_na(y))))
//       z = ifelse(is_na(y), NA_REAL, z);
//
//     return wrap(z);
//
//   } catch( std::exception &ex ) {
//     forward_exception_to_r( ex );
//   } catch(...) {
//     ::Rf_error( "c++ exception (unknown reason)" );
//   }
//   return R_NilValue; // -Wall
// }


/* Development version of gnorm.
 */
// SEXP gnorm(SEXP x, SEXP mean, SEXP sd) {
//   try {
//     NumericVector y(x);
//     NumericVector z = clone(y);
//     // NumericVector l_mean(mean);
//     // if(l_mean.size() >= 2)
//     //   Rcpp::stop("mean has to be scalar");
//
//     if(LENGTH(mean) >= 2)
//       Rcpp::stop("mean has to be a numeric scalar");
//     if(LENGTH(sd) >= 2)
//       Rcpp::stop("length(sd) != 1");
//
//     double l_mean = as<double>(mean);
//     double l_sd = as<double>(sd);
//
//     z = (y - l_mean) / l_sd;
//     z = ifelse((y == R_PosInf | y == R_NegInf), 0.,
//     	       -z / l_sd * dnorm(z));
//     // ifelse((y == R_PosInf | y == R_NegInf), 0.,
//     //        -(y - l_mean) / (l_sd * l_sd) *
//     //        dnorm(y, l_mean, l_sd));
//     if(is_true(any(is_na(y))))
//       z = ifelse(is_na(y), NA_REAL, z);
//
// //     z = ifelse(is_na(y), NA_REAL,
// // 	       ifelse((y == R_PosInf | y == R_NegInf), 0.,
// // 		      -(y - l_mean) / (l_sd * l_sd) *
// // 		      dnorm(y, l_mean, l_sd)));
//
//     // for(int i = 0; i < y.size(); i++) {
//     //   if(NumericVector::is_na(y[i]))
//     // 	z[i] = NA_REAL;
//     // }
//
//     // for(int i = 0; i < y.size(); i++) {
//     //   if(NumericVector::is_na(y[i]))
//     // 	z[i] = NA_REAL;
//     //   else if(y[i] == R_PosInf || y[i] == R_NegInf)
//     // 	z[i] = 0.;
//     //   else
//     // 	z[i] = - (y[i] - l_mean) / pow(l_sd, 2.) *
//     // 	  R::dnorm(y[i], l_mean, l_sd, 0);
//     // }
//     return wrap(z);
//
//   } catch( std::exception &ex ) {
//     forward_exception_to_r( ex );
//   } catch(...) {
//     ::Rf_error( "c++ exception (unknown reason)" );
//   }
//   return R_NilValue; // -Wall
// }

// NumericVector gnorm5(NumericVector x, double mean = 0., double sd = 1.) {
//   // int n = x.size();
//   NumericVector y = clone(x);
//
//   for(int i = 0; i < x.size(); i++) {
//     if(NumericVector::is_na(x[i]))
//       y[i] = NA_REAL;
//     else if(x[i] == R_PosInf || x[i] == R_NegInf)
//       y[i] = 0.;
//     else
//       y[i] = - (x[i] - mean) / pow(sd, 2.) *
// 	R::dnorm(x[i], mean, sd, 0);
//   }
//
//   return y;
// }

// SEXP rcpp_hello_world(){
//
//     CharacterVector x = CharacterVector::create( "foo", "bar" )  ;
//     NumericVector y   = NumericVector::create( 0.0, 1.0 ) ;
//     List z            = List::create( x, y ) ;
//
//     return z ;
// }

