#ifndef IREGNET_H
#define IREGNET_H

#include <Rcpp.h>
#include <cmath>

// used for all counts
#define ull unsigned long long

typedef enum {
  IREG_DIST_GAUSSIAN = 0,
  IREG_DIST_LOGISTIC,
  IREG_DIST_EXTREME_VALUE,
  IREG_DIST_LOG_GAUSSIAN,
  IREG_DIST_LOG_LOGISTIC,
  IREG_DIST_EXPONENTIAL,
  IREG_DIST_WEIBULL,
  IREG_DIST_UNKNOWN
} IREG_DIST;

/* Order set to agree with survival */
typedef enum {
  IREG_CENSOR_RIGHT = 0,
  IREG_CENSOR_NONE,
  IREG_CENSOR_LEFT,
  IREG_CENSOR_INTERVAL,
  IREG_CENSOR_INVALID
} IREG_CENSORING;

/* Functions from iregnet_fit.cpp */
typedef double (*transform_func) (double);

IREG_DIST
get_ireg_dist (Rcpp::String dist_str);

void
get_censoring_types (Rcpp::NumericMatrix &y, IREG_CENSORING *status);

/* Functions from distributions.cpp */
double
compute_grad_response(double *w, double *z, double *scale_update, const double *y_l, const double *y_r,
                      const double *eta, const double scale, const IREG_CENSORING *censoring_type,
                      const ull n_obs, IREG_DIST dist, double *mu, bool debug);

#endif  // IREGNET_H
