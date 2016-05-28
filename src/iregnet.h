#include <iostream>
#include <Rcpp.h>
#include <cmath>

// used for all counts
#define ull unsigned long long

typedef enum {
  IREG_DIST_GAUSSIAN = 0,
  IREG_DIST_LOGISTIC,
  IREG_DIST_EXTREME_VALUE,
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

typedef enum {
  IREG_DEBUG_NONE = 0,
  IREG_DEBUG_INPUT,     // 1
  IREG_DEBUG_N,
  IREG_DEBUG_CENSORING  // 3
} IREG_DEBUG;

/* Functions from distributions.cpp */
void compute_grad_response(double *w, double *z, double *y_l, double *y_r,
                           double *eta, double scale, IREG_CENSORING *censoring_type,
                           ull n_obs, IREG_DIST dist);
