/*
 * The log-likelihood and other parameters of the gaussian distributions will
 * be computed in this file.
 *
 * This file is copied (adapted) from the "SURVIVAL" package (survreg1.c) with
 * the kind permission of Dr. Terry Therneau.
 *
 */

/*
 ** < From Survival >
 **   j        ans[0]    ans[1]       ans[2]     ans[3]
 **   1                    f          f'/f        f''/ f
 **   2          F        1-F         f           f'
 **
 **  We do both F and 1-F to avoid the error in (1-F) for F near 1
 */
#include "R.h"
#include "iregnet.h"


#define SPI 2.506628274631001     /* sqrt(2*pi) */
#define ROOT_2 1.414213562373095
#define SMALL -200   /* what to use for log(f(x)) if f(x) is zero */
#define BIG_SIGMA_UPDATE 1 /* what to use for (log) scale_update if it should be very very large */


double
compute_only_none_censoring_type_logistic(rowvec *w, rowvec *z, double *scale_update, const double scale, const rowvec &eta,
                                 const rowvec &y_eta, const rowvec &y_eta_square, double &dsig_sum, double &ddsig_sum,
                                 const bool estimate_scale, const ull n_obs);

double
compute_grad_response_logistic_none(rowvec *w, rowvec *z, double *scale_update, const rowvec *y_l, const rowvec *y_r,
                                    const rowvec &eta, const double scale, const IREG_CENSORING *censoring_type,
                                    const ull n_obs, IREG_DIST dist, double *mu, bool debug, const bool estimate_scale,
                                    const rowvec &y_eta, const rowvec &y_eta_square,const int *separator, rowvec *tempvar) {
  double loglik = 0;
  double dsig_sum, ddsig_sum;
  dsig_sum = ddsig_sum = 0;

  loglik = compute_only_none_censoring_type_logistic(w, z, scale_update, scale, eta, y_eta, y_eta_square, dsig_sum, ddsig_sum,
                                            estimate_scale, n_obs);

  if (scale_update && estimate_scale == 1) {
    if (ddsig_sum != 0)
      *scale_update = -dsig_sum / ddsig_sum;
    else
      *scale_update = BIG_SIGMA_UPDATE;
  }

  return loglik;
}

double
compute_only_none_censoring_type_logistic(rowvec *w, rowvec *z, double *scale_update, const double scale, const rowvec &eta,
                                 const rowvec &y_eta, const rowvec &y_eta_square, double &dsig_sum, double &ddsig_sum,
                                 const bool estimate_scale, const ull n_obs)
{
  double scale_2 = scale * scale;
  double loglik = 0;
  double sum_log_scale = n_obs * log(scale);
  rowvec dg_vec;
  rowvec dsig_vec;
  rowvec ddsig_vec;

  rowvec temp_w_vec = rowvec(n_obs);
  rowvec sign_vec = rowvec(n_obs);
  rowvec sz_vec = scale * y_eta;

  for (int i = 0; i < n_obs; ++i) {
    if(y_eta(i)>0){
      temp_w_vec(i) = std::exp(-y_eta(i));
      sign_vec(i) = -1;
    } else {
      temp_w_vec(i) = std::exp(y_eta(i));
      sign_vec(i) = 1;
    }
  }

  rowvec temp_vec = 1+temp_w_vec;
  rowvec temp_square_vec = square(temp_vec);

  loglik = accu(log(temp_w_vec / temp_square_vec)) - sum_log_scale;

  (*w) = (-2 * temp_w_vec) / temp_square_vec / scale_2;

  dg_vec = ( ( sign_vec % (1-temp_w_vec) ) / temp_vec );
  dg_vec = -1 * ( dg_vec / scale );
  (*z) = eta - dg_vec /  (*w);

  if (scale_update) {
    dsig_vec = dg_vec % sz_vec;
    ddsig_vec = (*w) % square(sz_vec) - dsig_vec;
    dsig_vec -= 1;
    dsig_sum += accu(dsig_vec);
    ddsig_sum += accu(ddsig_vec);
  }

  return loglik;
}


