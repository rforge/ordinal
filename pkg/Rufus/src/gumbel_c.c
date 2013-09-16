#include "distributions_c.h"

void v_pgumbel(double *q, double *loc, double *scale, int *lower_tail)
// Consider implementing 'int give_log' to follow the convention from
// pnorm etc.
{
  if(!ISNAN(*q)) { // NaNs pass through unchanged
    if(*q == R_PosInf)
      *q = 1.;
    else if(*q == R_NegInf)
      *q = 0.;
    else {
      *q = (*q - *loc) / *scale;
      *q = exp( -exp( -*q));
    }
    if(!*lower_tail) *q = 1 - *q;
  }
}

double d_pgumbel(double q, double loc, double scale, int lower_tail)
{
  double x = q;
  v_pgumbel(&x, &loc, &scale, &lower_tail);
  return x;
}

void v_dgumbel(double *x, double *loc, double *scale, int *give_log)
{
  if(!ISNAN(*x)) {
    if(*x == R_PosInf || *x == R_NegInf)
      *x = 0;
    else {
      *x = (*x - *loc) / *scale;
      *x = -exp(-*x) - *x - log(*scale);
      if(!*give_log) *x = exp(*x);
    }
  }
}

double d_dgumbel(double q, double loc, double scale, int give_log)
{
  double x = q;
  v_dgumbel(&x, &loc, &scale, &give_log);
  return x;
}

void v_ggumbel(double *x)
{ // double y = *x
  if(*x == R_PosInf || *x == R_NegInf)
    *x = 0;
  else if(!ISNAN(*x)) {
    // double ey = exp(-y);
    *x = exp(-*x);
    if(*x == INFINITY)
      *x = 0;
    else {
      // double eq = exp(-exp(-y))
      double eq = exp(-*x);
      *x = -eq * *x + eq * *x * *x;
    }
  }
}

double d_ggumbel(double x)
{
  double y = x;
  v_ggumbel(&y);
  return y;
}
