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
                                    const rowvec &eta, const double scale, const IREG_CENSORING *censoring_type,
                                    const ull n_obs, IREG_DIST dist, double *mu, bool debug, const bool estimate_scale,
                                    const rowvec &y_eta, const rowvec &y_eta_square,const int *separator, rowvec *tempvar) {
  double scale_2 = scale * scale;
  double loglik;
  double dsig_sum, ddsig_sum;
  rowvec dsig_vec;
  rowvec ddsig_vec;
  double logspispecial = n_obs * log(SPI * scale);

  loglik = 0;
  dsig_sum = ddsig_sum = 0;

  loglik = accu(y_eta_square / -2) - logspispecial;
  (*z) = eta + scale * y_eta;
  (*w).fill(-1 / scale_2);

  if (scale_update && estimate_scale == 1) {
    dsig_vec = y_eta_square;//m
    /*ddsig_vec_one = (square(y_eta_one_square) - y_eta_one_square) / (scale_2 * scale_2)
                          - dsig_vec_one % (dsig_vec_one + 1);*/
    ddsig_vec = -2 * y_eta_square;
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

double
compute_grad_response_gaussian_right(rowvec *w, rowvec *z, double *scale_update, const rowvec *y_l, const rowvec *y_r,
                                     const rowvec &eta, const double scale, const IREG_CENSORING *censoring_type,
                                     const ull n_obs, IREG_DIST dist, double *mu, bool debug, const bool estimate_scale,
                                     const rowvec &y_eta, const rowvec &y_eta_square,const int *separator, rowvec *tempvar)
{
  double dsig_sum = 0, ddsig_sum = 0;
  double loglik, scale_2 = scale * scale;
  loglik = 0;

  /* we have to separate eta[] and y_eta[] into two different type by their censoring type
   * so that we can compute them by quickly way in armadillo
   * there are some temp vars below in this step
   */

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
  tempvar[10] = rowvec(none_censoring_type_number); // dsig_vec_one
  tempvar[11] = rowvec(right_censoring_type_number); // dsig_vec_two
  tempvar[12] = rowvec(none_censoring_type_number); // ddsig_vec_one
  tempvar[13] = rowvec(right_censoring_type_number); // ddsig_vec_two
  tempvar[14] = rowvec(right_censoring_type_number); // temp_densities
  tempvar[15] = rowvec(right_censoring_type_number); // dg_vec
  tempvar[16] = rowvec(right_censoring_type_number); // ddg_vec
  tempvar[17] = rowvec(right_censoring_type_number); // f_vec


  //for none censoring type
  loglik += accu(tempvar[3] / -2) - none_censoring_type_number * log(SPI * scale);
  tempvar[6] = tempvar[0] + scale * tempvar[2];
  tempvar[8].fill(-1 / scale_2);

  if (scale_update) {
    tempvar[10] = tempvar[3];//m
    /*tempvar[12] = (square(tempvar[3]) - tempvar[3]) / (scale_2 * scale_2)
                          - tempvar[10] % (tempvar[10] + 1);*/
    tempvar[12] = -2 * tempvar[3];
    // dsg = sz * temp2 - dg*(dsig +1);

    tempvar[10] -= 1;
    dsig_sum += accu(tempvar[10]);
    ddsig_sum += accu(tempvar[12]);
    /*if(m == 0 && n ==1) {
      std::cout << "dsig_sum :" << dsig_sum << "\n";
      std::cout << "ddsig_sum : " << ddsig_sum << "\n";
    }*/

    //dsig_vec = dsig_vec.array() - 1;

  }

  tempvar[17] = exp(-tempvar[5] / 2) /SPI;

  tempvar[14] = tempvar[4];

  tempvar[14].for_each( [](vec::elem_type& val) {
      if (val > 0) {
        //          densities_l[0] = (1 + erf(y_eta(i)/ROOT_2))/2;
        val =  erfc(val /ROOT_2) /2;
      }
      else {
        val = (1 + erf(-val /ROOT_2))/2;
        //          densities_l[0] =  erfc(-y_eta(i)/ROOT_2) /2;
      }
  } );

  /*if (y_eta(i)>0) {
  //          densities_l[0] = (1 + erf(y_eta(i)/ROOT_2))/2;
      densities_l[1] =  erfc(y_eta(i)/ROOT_2) /2;
    }
    else {
      densities_l[1] = (1 + erf(-y_eta(i)/ROOT_2))/2;
  //          densities_l[0] =  erfc(-y_eta(i)/ROOT_2) /2;
    }*/
  /* densities_l[2] = f;
   densities_l[3] = -y_eta(i)*f;*/
  loglik += accu(tempvar[5] / -2) - right_censoring_type_number * log(SPI);

  tempvar[15] = (tempvar[17] / tempvar[14]) / scale;
  tempvar[16] = ((tempvar[4] % tempvar[17]) / tempvar[14]) / scale_2;
  tempvar[16] -= square(tempvar[15]);
/*
  temp  = -densities_l[2] / densities_l[1] / scale;
  temp2 = -densities_l[3] / densities_l[1] / scale_2;
  dg = -temp; // dg = -(1/sigma) * (f'(z) / f(z))
  ddg = temp2 - dg * dg;     // f'(z^l) / f()] / [1-F()] ...*/

  if (scale_update) {
    tempvar[11] = (tempvar[17] % tempvar[4]) / tempvar[14];
    //tempvar[13] = sz * sz* temp2 - dsig * (1 + dsig);
    tempvar[13] = tempvar[17] % tempvar[4];
    tempvar[13] = (tempvar[14] % (tempvar[5] - 1) - tempvar[13]) % tempvar[13];
    tempvar[13] = tempvar[13] / square(tempvar[14]);

    dsig_sum += accu(tempvar[11]);
    ddsig_sum += accu(tempvar[13]);
  }

  tempvar[7] = tempvar[1] - tempvar[15] / tempvar[16];
  tempvar[9] = tempvar[16];

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

  double dsig_sum = 0, ddsig_sum = 0;
  double loglik, scale_2 = scale * scale;
  loglik = 0;

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
  tempvar[10] = rowvec(none_censoring_type_number); // dsig_vec_one
  tempvar[11] = rowvec(left_censoring_type_number); // dsig_vec_two
  tempvar[12] = rowvec(none_censoring_type_number); // ddsig_vec_one
  tempvar[13] = rowvec(left_censoring_type_number); // ddsig_vec_two
  tempvar[14] = rowvec(left_censoring_type_number); // temp_densities
  tempvar[15] = rowvec(left_censoring_type_number); // dg_vec
  tempvar[16] = rowvec(left_censoring_type_number); // ddg_vec
  tempvar[17] = rowvec(left_censoring_type_number); // f_vec


  //for none
  loglik += accu(tempvar[3] / -2) - none_censoring_type_number * log(SPI) - none_censoring_type_number * log(scale);
  tempvar[6] = tempvar[0] + scale * tempvar[2];
  tempvar[8].fill(-1 / scale_2);

  if (scale_update) {
    tempvar[10] = tempvar[3];//m
    /*ddsig_vec_one = (square(y_eta_one_square) - y_eta_one_square) / (scale_2 * scale_2)
                          - dsig_vec_one % (dsig_vec_one + 1);*/
    tempvar[12] = -2 * tempvar[3];
    // dsg = sz * temp2 - dg*(dsig +1);

    tempvar[10] -= 1;
    dsig_sum += accu(tempvar[10]);
    ddsig_sum += accu(tempvar[12]);
  }

  //for right
  tempvar[17] = exp(-tempvar[5] / 2) /SPI;

  tempvar[14] = tempvar[4];

  tempvar[14].for_each( [](vec::elem_type& val) {
      if (val > 0) {
        val =  (1 + erf(val /ROOT_2))/2;;
      }
      else {
        val = erfc(-val /ROOT_2) /2;
      }
  } );

  loglik += accu(tempvar[5] / -2) - left_censoring_type_number * log(SPI);

  tempvar[15] = -(tempvar[17] / tempvar[14]) / scale;
  tempvar[16] = -((tempvar[4] % tempvar[17]) / tempvar[14]) / scale_2;
  tempvar[16] -= square(tempvar[15]);

  if (scale_update) {
    tempvar[11] = -(tempvar[17] % tempvar[4]) / tempvar[14];
    //ddsig_vec_two = sz * sz* temp2 - dsig * (1 + dsig);
    tempvar[13] = tempvar[17] % tempvar[4];
    tempvar[13] = (tempvar[14] % (1 - tempvar[5]) - tempvar[13]) % tempvar[13];
    tempvar[13] = tempvar[13] / square(tempvar[14]);

    dsig_sum += accu(tempvar[11]);
    ddsig_sum += accu(tempvar[13]);
  }

  tempvar[7] = tempvar[1] - tempvar[15] / tempvar[16];
  tempvar[9] = tempvar[16];

  (*z) = join_rows(tempvar[6], tempvar[7]);
  (*w) = join_rows(tempvar[8], tempvar[9]);

  if (ddsig_sum != 0)
    *scale_update = -dsig_sum / ddsig_sum;
  else
    *scale_update = BIG_SIGMA_UPDATE;

  return loglik;
}
