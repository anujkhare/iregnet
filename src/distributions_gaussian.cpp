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
compute_only_none_censoring_type(rowvec *w, rowvec *z, double *scale_update, const double scale, const rowvec &eta,
                                 const rowvec &y_eta, const rowvec &y_eta_square, double &dsig_sum, double &ddsig_sum,
                                 const bool estimate_scale, const ull n_obs);

double
compute_left_or_right_censoring_type(rowvec *w, rowvec *z, double *scale_update, const double scale, const rowvec &eta,
                                  const rowvec &y_eta, const rowvec &y_eta_square, double &dsig_sum, double &ddsig_sum,
                                  const bool estimate_scale, const ull n_obs, const bool type);



double
compute_grad_response_gaussian_none(rowvec *w, rowvec *z, double *scale_update, const rowvec *y_l, const rowvec *y_r,
                                    const rowvec &eta, const double scale, const IREG_CENSORING *censoring_type,
                                    const ull n_obs, IREG_DIST dist, double *mu, bool debug, const bool estimate_scale,
                                    const rowvec &y_eta, const rowvec &y_eta_square,const int *separator, rowvec *tempvar) {
  double loglik;
  double dsig_sum, ddsig_sum;

  loglik = compute_only_none_censoring_type(w, z, scale_update, scale, eta, y_eta, y_eta_square, dsig_sum, ddsig_sum,
                                             estimate_scale, n_obs);

  return loglik;
}

double
compute_grad_response_gaussian_right(rowvec *w, rowvec *z, double *scale_update, const rowvec *y_l, const rowvec *y_r,
                                     const rowvec &eta, const double scale, const IREG_CENSORING *censoring_type,
                                     const ull n_obs, IREG_DIST dist, double *mu, bool debug, const bool estimate_scale,
                                     const rowvec &y_eta, const rowvec &y_eta_square,const int *separator, rowvec *tempvar)
{
  double loglik;
  double dsig_sum, ddsig_sum;

  ull none_censoring_type_number = separator[0];
  ull right_censoring_type_number = n_obs - separator[0];

  tempvar[0] = eta.subvec(0, separator[0] - 1); // eta_one
  tempvar[1] = eta.subvec(separator[0], n_obs - 1);// eta_two

  tempvar[2] = y_eta.subvec(0, separator[0] - 1); // y_eta_one
  tempvar[3] = y_eta_square.subvec(0, separator[0] - 1); // y_eta_one_square

  tempvar[4] = y_eta.subvec(separator[0], n_obs - 1); // y_eta_two
  tempvar[5] = y_eta_square.subvec(separator[0], n_obs - 1); // y_eta_two_square

  tempvar[6] = rowvec(none_censoring_type_number); // res_z_one
  tempvar[7] = rowvec(right_censoring_type_number); // res_z_two
  tempvar[8] = rowvec(none_censoring_type_number); // res_w_one
  tempvar[9] = rowvec(right_censoring_type_number); // res_w_two

  // compute for none
  loglik = compute_only_none_censoring_type(&tempvar[8], &tempvar[6], scale_update, scale, tempvar[0], tempvar[2], tempvar[3], dsig_sum, ddsig_sum,
                                            estimate_scale, none_censoring_type_number);

  // compute for right
  loglik += compute_left_or_right_censoring_type(&tempvar[9], &tempvar[7], scale_update, scale, tempvar[1], tempvar[4], tempvar[5], dsig_sum, ddsig_sum,
                                             estimate_scale, right_censoring_type_number, true);

  (*z) = join_rows(tempvar[6], tempvar[7]);
  (*w) = join_rows(tempvar[8], tempvar[9]);

  if (ddsig_sum != 0)
    *scale_update = -dsig_sum / ddsig_sum;
  else
    *scale_update = BIG_SIGMA_UPDATE;

  return loglik;
}


double
compute_grad_response_gaussian_left(rowvec *w, rowvec *z, double *scale_update, const rowvec *y_l, const rowvec *y_r,
                                     const rowvec &eta, const double scale, const IREG_CENSORING *censoring_type,
                                     const ull n_obs, IREG_DIST dist, double *mu, bool debug, const bool estimate_scale,
                                     const rowvec &y_eta, const rowvec &y_eta_square,const int *separator, rowvec *tempvar)
{
  double loglik;
  double dsig_sum, ddsig_sum;

  ull none_censoring_type_number = separator[0];
  ull left_censoring_type_number = n_obs - separator[0];

  tempvar[0] = eta.subvec(0, separator[0] - 1); // eta_one
  tempvar[1] = eta.subvec(separator[0], n_obs - 1);// eta_two

  tempvar[2] = y_eta.subvec(0, separator[0] - 1); // y_eta_one
  tempvar[3] = y_eta_square.subvec(0, separator[0] - 1); // y_eta_one_square

  tempvar[4] = y_eta.subvec(separator[0], n_obs - 1); // y_eta_two
  tempvar[5] = y_eta_square.subvec(separator[0], n_obs - 1); // y_eta_two_square

  tempvar[6] = rowvec(none_censoring_type_number); // res_z_one
  tempvar[7] = rowvec(left_censoring_type_number); // res_z_two
  tempvar[8] = rowvec(none_censoring_type_number); // res_w_one
  tempvar[9] = rowvec(left_censoring_type_number); // res_w_two

  // compute for none
  loglik = compute_only_none_censoring_type(&tempvar[8], &tempvar[6], scale_update, scale, tempvar[0], tempvar[2], tempvar[3], dsig_sum, ddsig_sum,
                                            estimate_scale, none_censoring_type_number);

  // compute for left
  loglik += compute_left_or_right_censoring_type(&tempvar[9], &tempvar[7], scale_update, scale, tempvar[1], tempvar[4], tempvar[5], dsig_sum, ddsig_sum,
                                                 estimate_scale, left_censoring_type_number, false);

  (*z) = join_rows(tempvar[6], tempvar[7]);
  (*w) = join_rows(tempvar[8], tempvar[9]);

  if (ddsig_sum != 0)
    *scale_update = -dsig_sum / ddsig_sum;
  else
    *scale_update = BIG_SIGMA_UPDATE;

  return loglik;
}

double
compute_only_none_censoring_type(rowvec *w, rowvec *z, double *scale_update, const double scale, const rowvec &eta,
                                 const rowvec &y_eta, const rowvec &y_eta_square, double &dsig_sum, double &ddsig_sum,
                                 const bool estimate_scale, const ull n_obs)
{
  double scale_2 = scale * scale;
  double loglik;
  double logspispecial = n_obs * log(SPI * scale);
  rowvec dsig_vec;
  rowvec ddsig_vec;

  loglik = 0;
  dsig_sum = ddsig_sum = 0;

  loglik = accu(y_eta_square / -2) - logspispecial;
  (*z) = eta + scale * y_eta;
  (*w).fill(-1 / scale_2);

  if (scale_update && estimate_scale == 1) {
    dsig_vec = y_eta_square;
    ddsig_vec = -2 * y_eta_square;
    dsig_vec -= 1;
    dsig_sum += accu(dsig_vec);
    ddsig_sum += accu(ddsig_vec);
    if (ddsig_sum != 0)
      *scale_update = -dsig_sum / ddsig_sum;
    else
      *scale_update = BIG_SIGMA_UPDATE;
  }

  return loglik;
}

double
compute_left_or_right_censoring_type(rowvec *w, rowvec *z, double *scale_update, const double scale, const rowvec &eta,
                                     const rowvec &y_eta, const rowvec &y_eta_square, double &dsig_sum, double &ddsig_sum,
                                     const bool estimate_scale, const ull n_obs, const bool type)
{

  double scale_2 = scale * scale;
  double loglik = 0;
  rowvec dsig_vec;
  rowvec ddsig_vec;

  rowvec temp_densities(n_obs);
  rowvec dg_vec(n_obs);
  rowvec ddg_vec(n_obs);
  rowvec f_vec(n_obs);

  f_vec = exp(-y_eta_square / 2) /SPI;

  temp_densities = y_eta;

  // if type is true, right censoring is computed
  // else left censoring is computed
  if(type){

    temp_densities.for_each( [](vec::elem_type& val) {

        if (val > 0) {
          val =  erfc(val /ROOT_2) /2;
        }
        else {
          val = (1 + erf(-val /ROOT_2))/2;
        }

    } );

    loglik += accu(y_eta_square / -2) - n_obs * log(SPI);

    dg_vec = (f_vec / temp_densities) / scale;
    ddg_vec = ((y_eta % f_vec) / temp_densities) / scale_2;
    ddg_vec -= square(dg_vec);

    if (scale_update) {
      dsig_vec = (f_vec % y_eta) / temp_densities;
      ddsig_vec = f_vec % y_eta;
      ddsig_vec = (temp_densities % (y_eta_square - 1) - ddsig_vec) % ddsig_vec;
      ddsig_vec = ddsig_vec / square(temp_densities);

      dsig_sum += accu(dsig_vec);
      ddsig_sum += accu(ddsig_vec);
    }

  } else {

    temp_densities.for_each( [](vec::elem_type& val) {

        if (val > 0) {
          val =  (1 + erf(val /ROOT_2))/2;;
        }
        else {
          val = erfc(-val /ROOT_2) /2;
        }

    } );

    loglik += accu(y_eta_square / -2) - n_obs * log(SPI);

    dg_vec = -(f_vec / temp_densities) / scale;
    ddg_vec = -((y_eta % f_vec) / temp_densities) / scale_2;
    ddg_vec -= square(dg_vec);

    if (scale_update) {
      dsig_vec = -(f_vec % y_eta) / temp_densities;
      ddsig_vec = f_vec % y_eta;
      ddsig_vec = (temp_densities % (1 - y_eta_square) - ddsig_vec) % ddsig_vec;
      ddsig_vec = ddsig_vec / square(temp_densities);

      dsig_sum += accu(dsig_vec);
      ddsig_sum += accu(ddsig_vec);
    }
  }

  (*z) = eta - dg_vec / ddg_vec;
  (*w) = ddg_vec;

  return loglik;
}


