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
compute_grad_response_gaussian_none(rowvec *w, rowvec *z, double *scale_update, const rowvec *y_l, const rowvec *y_r,
                                    const rowvec *eta, const double scale, const IREG_CENSORING *censoring_type,
                                    const ull n_obs, IREG_DIST dist, double *mu, bool debug, const bool estimate_scale,
                                    rowvec *y_eta, rowvec *y_eta_square) {
  double scale_2 = scale * scale;
  double loglik;
  double dsig_sum, ddsig_sum;
  rowvec dsig_vec;
  rowvec ddsig_vec;
  double logspispecial = n_obs * log(SPI * scale);

  loglik = 0;
  dsig_sum = ddsig_sum = 0;

  loglik = accu((*y_eta_square) / -2) - logspispecial;
  (*z) = (*eta) + scale * (*y_eta);
  (*w).fill(-1 / scale_2);

  if (scale_update && estimate_scale == 1) {
    dsig_vec = (*y_eta_square);//m
    /*ddsig_vec_one = (square(y_eta_one_square) - y_eta_one_square) / (scale_2 * scale_2)
                          - dsig_vec_one % (dsig_vec_one + 1);*/
    ddsig_vec = -2 * (*y_eta_square);
    // dsg = sz * temp2 - dg*(dsig +1);
    dsig_vec -= 1;
    dsig_sum += accu(dsig_vec);
    ddsig_sum += accu(ddsig_vec);
    //dsig_vec = dsig_vec.array() - 1;
    if (ddsig_sum != 0)
      *scale_update = -dsig_sum / ddsig_sum;
    else
      *scale_update = BIG_SIGMA_UPDATE;
  }

  return loglik;
}
