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
compute_left_or_right_censoring_type_logistic(rowvec *w, rowvec *z, double *scale_update, const double scale, const rowvec &eta,
                                              const rowvec &y_eta, const rowvec &y_eta_square, double &dsig_sum, double &ddsig_sum,
                                              const bool estimate_scale, const ull n_obs, const bool type);

double
compute_only_interval_censoring_type_logistic(rowvec *w, rowvec *z, double *scale_update, const double scale, const rowvec &eta,
                                              const rowvec &y_eta, const rowvec &y_eta_square, const rowvec &y_eta_r,
                                              const rowvec &y_eta_square_r, double &dsig_sum, double &ddsig_sum,
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
compute_grad_response_logistic_right(rowvec *w, rowvec *z, double *scale_update, const rowvec *y_l, const rowvec *y_r,
                                     const rowvec &eta, const double scale, const IREG_CENSORING *censoring_type,
                                     const ull n_obs, IREG_DIST dist, double *mu, bool debug, const bool estimate_scale,
                                     const rowvec &y_eta, const rowvec &y_eta_square,const int *separator, rowvec *tempvar)
{
  double loglik = 0;
  double dsig_sum, ddsig_sum;
  dsig_sum = ddsig_sum = 0;

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
  loglik = compute_only_none_censoring_type_logistic(&tempvar[8], &tempvar[6], scale_update, scale, tempvar[0], tempvar[2],
                                            tempvar[3], dsig_sum, ddsig_sum, estimate_scale, none_censoring_type_number);

  // compute for right
  loglik += compute_left_or_right_censoring_type_logistic(&tempvar[9], &tempvar[7], scale_update, scale, tempvar[1], tempvar[4],
                                                 tempvar[5], dsig_sum, ddsig_sum, estimate_scale,
                                                 right_censoring_type_number, true);

  (*z) = join_rows(tempvar[6], tempvar[7]);
  (*w) = join_rows(tempvar[8], tempvar[9]);

  if (scale_update && estimate_scale == 1) {
    if (ddsig_sum != 0)
      *scale_update = -dsig_sum / ddsig_sum;
    else
      *scale_update = BIG_SIGMA_UPDATE;
  }

  return loglik;
}


double
compute_grad_response_logistic_left(rowvec *w, rowvec *z, double *scale_update, const rowvec *y_l, const rowvec *y_r,
                                    const rowvec &eta, const double scale, const IREG_CENSORING *censoring_type,
                                    const ull n_obs, IREG_DIST dist, double *mu, bool debug, const bool estimate_scale,
                                    const rowvec &y_eta, const rowvec &y_eta_square,const int *separator, rowvec *tempvar)
{
  double loglik = 0;
  double dsig_sum, ddsig_sum;
  dsig_sum = ddsig_sum = 0;

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
  loglik = compute_only_none_censoring_type_logistic(&tempvar[8], &tempvar[6], scale_update, scale, tempvar[0], tempvar[2],
                                            tempvar[3], dsig_sum, ddsig_sum, estimate_scale, none_censoring_type_number);

  // compute for left
  loglik += compute_left_or_right_censoring_type_logistic(&tempvar[9], &tempvar[7], scale_update, scale, tempvar[1], tempvar[4],
                                                 tempvar[5], dsig_sum, ddsig_sum, estimate_scale,
                                                 left_censoring_type_number, false);

  (*z) = join_rows(tempvar[6], tempvar[7]);
  (*w) = join_rows(tempvar[8], tempvar[9]);

  if (scale_update && estimate_scale == 1) {
    if (ddsig_sum != 0)
      *scale_update = -dsig_sum / ddsig_sum;
    else
      *scale_update = BIG_SIGMA_UPDATE;
  }

  return loglik;
}

double
compute_grad_response_logistic_interval(rowvec *w, rowvec *z, double *scale_update, const rowvec *y_l, const rowvec *y_r,
                                        const rowvec &eta, const double scale, const IREG_CENSORING *censoring_type,
                                        const ull n_obs, IREG_DIST dist, double *mu, bool debug, const bool estimate_scale,
                                        const rowvec &y_eta, const rowvec &y_eta_square,const int *separator, rowvec *tempvar)
{
  double loglik = 0;
  double dsig_sum, ddsig_sum;
  dsig_sum = ddsig_sum = 0;

  ull none_censoring_type_number = separator[0];
  ull left_censoring_type_number = separator[1];
  ull right_censoring_type_number = separator[2];
  ull interval_censoring_type_number = separator[3];

  ull separator_none_and_left = separator[0];
  ull separator_left_and_right = separator[0] + separator[1];
  ull separator_right_and_interval = separator[0] + separator[1] + separator[2];

  tempvar[0] = rowvec(none_censoring_type_number); // res_z_none
  tempvar[1] = rowvec(left_censoring_type_number); // res_z_left
  tempvar[2] = rowvec(right_censoring_type_number); // res_z_right
  tempvar[3] = rowvec(interval_censoring_type_number); // res_z_interval

  tempvar[4] = rowvec(none_censoring_type_number); // res_w_none
  tempvar[5] = rowvec(left_censoring_type_number); // res_w_left
  tempvar[6] = rowvec(right_censoring_type_number); // res_w_right
  tempvar[7] = rowvec(interval_censoring_type_number); // res_w_interval

  rowvec y_r_interval = (*y_r).subvec(separator_right_and_interval, n_obs - 1);
  rowvec y_eta_interval_r = (y_r_interval - eta.subvec(separator_right_and_interval, n_obs - 1)) / scale;
  rowvec y_eta_interval_square_r = square(y_eta_interval_r);

  (*z) = rowvec();
  (*w) = rowvec();


  // compute for none
  if(none_censoring_type_number > 0){

    loglik += compute_only_none_censoring_type_logistic(&tempvar[4], &tempvar[0], scale_update, scale,
                                               eta.subvec(0, separator_none_and_left - 1),
                                               y_eta.subvec(0, separator_none_and_left - 1),
                                               y_eta_square.subvec(0, separator_none_and_left - 1),
                                               dsig_sum, ddsig_sum, estimate_scale, none_censoring_type_number);

    (*z) = join_rows((*z), tempvar[0]);
    (*w) = join_rows((*w), tempvar[4]);
  }
  // compute for left
  if(left_censoring_type_number > 0){
    loglik += compute_left_or_right_censoring_type_logistic(&tempvar[5], &tempvar[1], scale_update, scale,
                                                   eta.subvec(separator_none_and_left, separator_left_and_right - 1),
                                                   y_eta.subvec(separator_none_and_left, separator_left_and_right - 1),
                                                   y_eta_square.subvec(separator_none_and_left, separator_left_and_right - 1),
                                                   dsig_sum, ddsig_sum, estimate_scale, left_censoring_type_number, false);
    (*z) = join_rows((*z), tempvar[1]);
    (*w) = join_rows((*w), tempvar[5]);
  }
  // compute for right
  if(right_censoring_type_number > 0){
    loglik += compute_left_or_right_censoring_type_logistic(&tempvar[6], &tempvar[2], scale_update, scale,
                                                   eta.subvec(separator_left_and_right, separator_right_and_interval - 1),
                                                   y_eta.subvec(separator_left_and_right, separator_right_and_interval - 1),
                                                   y_eta_square.subvec(separator_left_and_right, separator_right_and_interval - 1),
                                                   dsig_sum, ddsig_sum, estimate_scale, right_censoring_type_number, true);
    (*z) = join_rows((*z), tempvar[2]);
    (*w) = join_rows((*w), tempvar[6]);
  }

  // compute for interval
  loglik += compute_only_interval_censoring_type_logistic(&tempvar[7], &tempvar[3], scale_update, scale,
                                                 eta.subvec(separator_right_and_interval, n_obs - 1),
                                                 y_eta.subvec(separator_right_and_interval, n_obs - 1),
                                                 y_eta_square.subvec(separator_right_and_interval, n_obs - 1),
                                                 y_eta_interval_r, y_eta_interval_square_r,
                                                 dsig_sum, ddsig_sum, estimate_scale, interval_censoring_type_number);
  (*z) = join_rows((*z), tempvar[3]);
  (*w) = join_rows((*w), tempvar[7]);

  if (scale_update) {
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

double
compute_left_or_right_censoring_type_logistic(rowvec *w, rowvec *z, double *scale_update, const double scale, const rowvec &eta,
                                     const rowvec &y_eta, const rowvec &y_eta_square, double &dsig_sum, double &ddsig_sum,
                                     const bool estimate_scale, const ull n_obs, const bool type)
{
  double scale_2 = scale * scale;
  double loglik = 0;
  rowvec dsig_vec;
  rowvec ddsig_vec;
  rowvec dg_vec(n_obs);
  rowvec ddg_vec(n_obs);

  rowvec temp_w_vec = rowvec(n_obs);
  rowvec sign_vec = rowvec(n_obs);
  rowvec sz_vec = scale * y_eta;
  rowvec temp_1 = rowvec(n_obs);
  rowvec temp_2 = rowvec(n_obs);
  rowvec w_base = rowvec(n_obs);
  rowvec w_base_square = rowvec(n_obs);
  rowvec f_vec = rowvec(n_obs);

  // if type is true, right censoring is computed
  // else left censoring is computed

  if(type){

    for (int i = 0; i < n_obs; ++i) {
      if(y_eta(i)>0){
        w_base(i) = std::exp(-y_eta(i));
        f_vec(i) = w_base(i) / (1 + w_base(i));
        sign_vec(i) = -1;
      } else {
        w_base(i) = std::exp(y_eta(i));
        f_vec(i) = 1 / (1 + w_base(i));
        sign_vec(i) = 1;
      }
    }

    w_base_square = square(1 + w_base);

    temp_1 = -1 * ((w_base / w_base_square) / f_vec / scale);
    temp_2 = sign_vec % w_base % (1 - w_base);
    temp_2 = -1 * temp_2 / pow(1 + w_base, 3) / f_vec / scale_2;

    loglik = accu(log(f_vec));

    dg_vec = -1 * temp_1;
    ddg_vec = temp_2 - square(dg_vec);

    if (scale_update) {
      dsig_vec = dg_vec % (scale * y_eta);
      ddsig_vec = scale_2 * y_eta_square;
      ddsig_vec = ddsig_vec % temp_2;
      ddsig_vec = ddsig_vec - dsig_vec % (1 + dsig_vec);

      dsig_sum += accu(dsig_vec);
      ddsig_sum += accu(ddsig_vec);
    }

  } else {

    for (int i = 0; i < n_obs; ++i) {
      if(y_eta(i)>0){
        w_base(i) = std::exp(-y_eta(i));
        f_vec(i) = 1 / (1 + w_base(i));
        sign_vec(i) = -1;
      } else {
        w_base(i) = std::exp(y_eta(i));
        f_vec(i) = w_base(i) / (1 + w_base(i));
        sign_vec(i) = 1;
      }
    }

    w_base_square = square(1 + w_base);

    temp_1 = ((w_base / w_base_square) / f_vec / scale);
    temp_2 = sign_vec % w_base % (1 - w_base);
    temp_2 = temp_2 / pow(1 + w_base, 3) / f_vec / scale_2;

    loglik = accu(log(f_vec));

    dg_vec = -1 * temp_1;
    ddg_vec = temp_2 - square(dg_vec);

    if (scale_update) {
      dsig_vec = dg_vec % (scale * y_eta);
      ddsig_vec = scale_2 * y_eta_square;
      ddsig_vec = ddsig_vec % temp_2;
      ddsig_vec = ddsig_vec - dsig_vec % (1 + dsig_vec);

      dsig_sum += accu(dsig_vec);
      ddsig_sum += accu(ddsig_vec);
    }
  }

  (*z) = eta - dg_vec / ddg_vec;
  (*w) = ddg_vec;

  return loglik;
}

double
compute_only_interval_censoring_type_logistic(rowvec *w, rowvec *z, double *scale_update, const double scale, const rowvec &eta,
                                     const rowvec &y_eta, const rowvec &y_eta_square, const rowvec &y_eta_r,
                                     const rowvec &y_eta_square_r, double &dsig_sum, double &ddsig_sum,
                                     const bool estimate_scale, const ull n_obs)
{
  double scale_2 = scale * scale;
  double loglik = 0;
  rowvec dsig_vec;
  rowvec ddsig_vec;

  rowvec dg_vec(n_obs);
  rowvec ddg_vec(n_obs);
  rowvec f_vec_i_l(n_obs);
  rowvec f_vec_i_r(n_obs);

  rowvec w_base_l(n_obs);
  rowvec w_base_r(n_obs);
  rowvec sign_l(n_obs);
  rowvec sign_r(n_obs);
  rowvec temp_w_l(n_obs);
  rowvec temp_w_r(n_obs);

  rowvec temp(n_obs);
  rowvec temp2(n_obs);
  rowvec d_l(n_obs);
  rowvec d_r(n_obs);

  for (int i = 0; i < n_obs; ++i) {
    if(y_eta(i) > 0 && y_eta_r(i) > 0){

      w_base_l(i) = std::exp(-y_eta(i));
      w_base_r(i) = std::exp(-y_eta_r(i));
      temp_w_l(i) = 1 + w_base_l(i);
      temp_w_r(i) = 1 + w_base_r(i);

      temp(i) = w_base_l(i)/temp_w_l(i) - w_base_r(i)/temp_w_r(i);
      sign_l(i) = -1;
      sign_r(i) = -1;
    } else if (y_eta(i) > 0 && y_eta_r(i) <= 0){

      w_base_l(i) = std::exp(-y_eta(i));
      w_base_r(i) = std::exp(y_eta_r(i));
      temp_w_l(i) = 1 + w_base_l(i);
      temp_w_r(i) = 1 + w_base_r(i);

      temp(i) = w_base_l(i)/temp_w_l(i) - 1/temp_w_r(i);
      sign_l(i) = -1;
      sign_r(i) = 1;
    } else if (y_eta(i) <= 0 && y_eta_r(i) > 0){

      w_base_l(i) = std::exp(y_eta(i));
      w_base_r(i) = std::exp(-y_eta_r(i));
      temp_w_l(i) = 1 + w_base_l(i);
      temp_w_r(i) = 1 + w_base_r(i);

      temp(i) = 1/temp_w_r(i) - w_base_l(i)/temp_w_l(i);
      sign_l(i) = 1;
      sign_r(i) = -1;
    } else {

      w_base_l(i) = std::exp(y_eta(i));
      w_base_r(i) = std::exp(y_eta_r(i));
      temp_w_l(i) = 1 + w_base_l(i);
      temp_w_r(i) = 1 + w_base_r(i);

      temp(i) = w_base_r(i)/temp_w_r(i) - w_base_l(i)/temp_w_l(i);
      sign_l(i) = 1;
      sign_r(i) = 1;
    }
  }

  loglik += accu(log(temp));

  f_vec_i_l = w_base_l/square(temp_w_l);
  f_vec_i_r = w_base_r/square(temp_w_r);
  dg_vec = - (f_vec_i_r - f_vec_i_l) / (scale * temp);

  d_l = sign_l % f_vec_i_l % (1 - w_base_l) / temp_w_l;
  d_r = sign_r % f_vec_i_r % (1 - w_base_r) / temp_w_r;

  ddg_vec = (d_r - d_l) / (scale_2 * temp) - square(dg_vec);

  if(scale_update){
    dsig_vec = (y_eta % f_vec_i_l - y_eta_r % f_vec_i_r) / temp;
    ddsig_vec = ((y_eta_square_r % d_r + y_eta_square % d_l) / temp) - dsig_vec % (1 + dsig_vec);

    dsig_sum += accu(dsig_vec);
    ddsig_sum += accu(ddsig_vec);
  }

  if(any(ddg_vec == 0)){
    if(all(ddg_vec == 0))
      (*z) = eta;
    else{
      for (int i = 0; i < n_obs; ++i) {
        if(ddg_vec(i) == 0)
          (*z) = eta(i);
        else
          (*z) = eta(i) - dg_vec(i) / ddg_vec(i);
      }
    }

  } else{
    (*z) = eta - dg_vec / ddg_vec;
  }
  (*z) = eta - dg_vec / ddg_vec;
  (*w) = ddg_vec;

  return loglik;
}

