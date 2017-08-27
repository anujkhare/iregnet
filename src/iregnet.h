#ifndef IREGNET_H
#define IREGNET_H

#include <RcppArmadillo.h>
#include <cmath>

using namespace arma;

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
get_censoring_types (mat &y, IREG_CENSORING *status);

/* Functions from distributions.cpp */
double
compute_grad_response(rowvec *w, rowvec *z, double *scale_update, const rowvec *y_l, const rowvec *y_r,
                      const rowvec &eta, const double scale, const IREG_CENSORING *censoring_type,
                      const ull n_obs, IREG_DIST dist, double *mu, bool debug, const bool estimate_scale, const rowvec &y_eta,
                      const rowvec &y_eta_square, const int *separator, rowvec *tempvar);

/* Functions from distributions_gaussian.cpp */
double
compute_grad_response_gaussian_none(rowvec *w, rowvec *z, double *scale_update, const rowvec *y_l, const rowvec *y_r,
                                    const rowvec &eta, const double scale, const IREG_CENSORING *censoring_type,
                                    const ull n_obs, IREG_DIST dist, double *mu, bool debug, const bool estimate_scale,
                                    const rowvec &y_eta, const rowvec &y_eta_square,const int *separator, rowvec *tempvar);

double
compute_grad_response_gaussian_right(rowvec *w, rowvec *z, double *scale_update, const rowvec *y_l, const rowvec *y_r,
                                    const rowvec &eta, const double scale, const IREG_CENSORING *censoring_type,
                                    const ull n_obs, IREG_DIST dist, double *mu, bool debug, const bool estimate_scale,
                                     const rowvec &y_eta, const rowvec &y_eta_square,const int *separator, rowvec *tempvar);

double
compute_grad_response_gaussian_left(rowvec *w, rowvec *z, double *scale_update, const rowvec *y_l, const rowvec *y_r,
                                     const rowvec &eta, const double scale, const IREG_CENSORING *censoring_type,
                                     const ull n_obs, IREG_DIST dist, double *mu, bool debug, const bool estimate_scale,
                                    const rowvec &y_eta, const rowvec &y_eta_square,const int *separator, rowvec *tempvar);

double
compute_grad_response_gaussian_interval(rowvec *w, rowvec *z, double *scale_update, const rowvec *y_l, const rowvec *y_r,
                                        const rowvec &eta, const double scale, const IREG_CENSORING *censoring_type,
                                        const ull n_obs, IREG_DIST dist, double *mu, bool debug, const bool estimate_scale,
                                        const rowvec &y_eta, const rowvec &y_eta_square,const int *separator, rowvec *tempvar);

/* Functions from distributions_logistic.cpp */
double
compute_grad_response_logistic_none(rowvec *w, rowvec *z, double *scale_update, const rowvec *y_l, const rowvec *y_r,
                                    const rowvec &eta, const double scale, const IREG_CENSORING *censoring_type,
                                    const ull n_obs, IREG_DIST dist, double *mu, bool debug, const bool estimate_scale,
                                    const rowvec &y_eta, const rowvec &y_eta_square,const int *separator, rowvec *tempvar);

double
compute_grad_response_logistic_right(rowvec *w, rowvec *z, double *scale_update, const rowvec *y_l, const rowvec *y_r,
                                     const rowvec &eta, const double scale, const IREG_CENSORING *censoring_type,
                                     const ull n_obs, IREG_DIST dist, double *mu, bool debug, const bool estimate_scale,
                                     const rowvec &y_eta, const rowvec &y_eta_square,const int *separator, rowvec *tempvar);

double
compute_grad_response_logistic_left(rowvec *w, rowvec *z, double *scale_update, const rowvec *y_l, const rowvec *y_r,
                                    const rowvec &eta, const double scale, const IREG_CENSORING *censoring_type,
                                    const ull n_obs, IREG_DIST dist, double *mu, bool debug, const bool estimate_scale,
                                    const rowvec &y_eta, const rowvec &y_eta_square,const int *separator, rowvec *tempvar);

double
compute_grad_response_logistic_interval(rowvec *w, rowvec *z, double *scale_update, const rowvec *y_l, const rowvec *y_r,
                                        const rowvec &eta, const double scale, const IREG_CENSORING *censoring_type,
                                        const ull n_obs, IREG_DIST dist, double *mu, bool debug, const bool estimate_scale,
                                        const rowvec &y_eta, const rowvec &y_eta_square,const int *separator, rowvec *tempvar);


#endif  // IREGNET_H
