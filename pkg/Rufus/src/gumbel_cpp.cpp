#include "distributions_cpp.h"
#include "distributions_c.h"

using namespace Rcpp; // rather than repeating it in each function.


/*-------------------------------------------------------
  Cumulative distribution functions (CDFs)
  -------------------------------------------------------*/

SEXP pgumbel(SEXP q, SEXP location, SEXP scale, SEXP lower_tail, SEXP max)
{
  BEGIN_RCPP
  NumericVector x = clone(q);

  double loc_ = as<double>(location);
  double scale_= as<double>(scale);
  int lower_ = as<int>(lower_tail);
  int max_ = as<int>(max);

  if(!max_ ) x = -x;
  for(int i = 0; i < x.size(); i++)
    v_pgumbel(&x[i], &loc_, &scale_, &lower_);

  return wrap(x);
  END_RCPP
}

/*-------------------------------------------------------
  Probability mass functions (PDFs)
  -------------------------------------------------------*/

SEXP dgumbel(SEXP x, SEXP location, SEXP scale, SEXP give_log, SEXP max)
{
  BEGIN_RCPP
  NumericVector q = clone(x);

  double loc_ = as<double>(location);
  double scale_= as<double>(scale);
  int give_log_ = as<int>(give_log);
  int max_ = as<int>(max);

  if(!max_ ) q = -q;
  for(int i = 0; i < q.size(); i++)
    v_dgumbel(&q[i], &loc_, &scale_, &give_log_);

  return wrap(q);
  END_RCPP
}

/*-------------------------------------------------------
  Gradients of PDFs
  -------------------------------------------------------*/

SEXP ggumbel(SEXP x, SEXP max)
{
  BEGIN_RCPP
  NumericVector z = clone(x);
  int max_ = as<int>(max);
  if(!max_ ) z = -z;
  v_ggumbel(z);
  return wrap(z);
  END_RCPP
}


void v_ggumbel(NumericVector x)
{
  for(int i = 0; i < x.size(); i++)
    v_ggumbel(&x[i]);
}

NumericVector ggumbel(NumericVector x)
{
  NumericVector y = clone(x);
  v_ggumbel(y);
  return y;
}

