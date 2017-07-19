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
                                    rowvec &y_eta, rowvec &y_eta_square, int *separator, rowvec *tempvar) {
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
                                     rowvec &y_eta, rowvec &y_eta_square, int *separator, rowvec *tempvar)
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
                                    rowvec &y_eta, rowvec &y_eta_square, int *separator, rowvec *tempvar)
{

  double dsig_sum = 0, ddsig_sum = 0;
  double loglik, scale_2 = scale * scale;

  rowvec eta_one(n_obs);
  rowvec eta_two(n_obs);
  rowvec y_eta_one(n_obs);
  rowvec y_eta_one_square(n_obs);
  rowvec y_eta_two(n_obs);
  rowvec y_eta_two_square(n_obs);

  double count_one = 0;
  double count_two = 0;

  loglik = 0;

  for (int i = 0; i < n_obs; ++i) {
    if(censoring_type[i] == 1){
      y_eta_one(count_one) = y_eta(i);
      y_eta_one_square(count_one) = y_eta_square(i);
      eta_one(count_one) = eta(i);
      count_one ++;
    } else if(censoring_type[i] == 2){
      y_eta_two(count_two) = y_eta(i);
      y_eta_two_square(count_two) = y_eta_square(i);
      eta_two(count_two) = eta(i);
      count_two ++;
    }
  }

  y_eta_one.resize(count_one);
  y_eta_one_square.resize(count_one);
  eta_one.resize(count_one);

  y_eta_two.resize(count_two);
  y_eta_two_square.resize(count_two);
  eta_two.resize(count_two);

  rowvec res_z_one(count_one);
  rowvec res_z_two(count_two);

  rowvec res_w_one(count_one);
  rowvec res_w_two(count_two);

  rowvec dsig_vec_one(count_one);
  rowvec dsig_vec_two(count_two);

  rowvec ddsig_vec_one(count_one);
  rowvec ddsig_vec_two(count_two);


  //for none
  loglik += accu(y_eta_one_square / -2) - count_one * log(SPI) - count_one * log(scale);
  res_z_one = eta_one + scale * y_eta_one;
  res_w_one.fill(-1 / scale_2);

  if (scale_update) {
    dsig_vec_one = y_eta_one_square;//m
    /*ddsig_vec_one = (square(y_eta_one_square) - y_eta_one_square) / (scale_2 * scale_2)
                          - dsig_vec_one % (dsig_vec_one + 1);*/
    ddsig_vec_one = -2 * y_eta_one_square;
    // dsg = sz * temp2 - dg*(dsig +1);

    dsig_vec_one -= 1;
    dsig_sum += accu(dsig_vec_one);
    ddsig_sum += accu(ddsig_vec_one);
  }

  //for right
  rowvec temp_densities(count_two);
  rowvec dg_vec(count_two);
  rowvec ddg_vec(count_two);
  rowvec f_vec(count_two);

  f_vec = exp(-y_eta_two_square / 2) /SPI;

  temp_densities = y_eta_two;

  temp_densities.for_each( [](vec::elem_type& val) {
      if (val > 0) {
        val =  (1 + erf(val /ROOT_2))/2;;
      }
      else {
        val = erfc(-val /ROOT_2) /2;
      }
  } );

  loglik += accu(y_eta_two_square / -2) - count_two * log(SPI);

  dg_vec = -(f_vec / temp_densities) / scale;
  ddg_vec = -((y_eta_two % f_vec) / temp_densities) / scale_2;
  ddg_vec -= square(dg_vec);

  if (scale_update) {
    dsig_vec_two = -(f_vec % y_eta_two) / temp_densities;
    //ddsig_vec_two = sz * sz* temp2 - dsig * (1 + dsig);
    ddsig_vec_two = f_vec % y_eta_two;
    ddsig_vec_two = (temp_densities % (1 - y_eta_two_square) - ddsig_vec_two) % ddsig_vec_two;
    ddsig_vec_two = ddsig_vec_two / square(temp_densities);

    dsig_sum += accu(dsig_vec_two);
    ddsig_sum += accu(ddsig_vec_two);
  }

  res_z_two = eta_two - dg_vec / ddg_vec;
  res_w_two = ddg_vec;

  count_one = 0;
  count_two = 0;

  for (int i = 0; i < n_obs; ++i) {
    if(censoring_type[i] == 1){
      (*z)(i) = res_z_one(count_one);
      (*w)(i) = res_w_one(count_one);
      count_one ++;
    } else if(censoring_type[i] == 2){
      (*z)(i) = res_z_two(count_two);
      (*w)(i) = res_w_two(count_two);
      count_two ++;
    }
  }

  if (ddsig_sum != 0)
    *scale_update = -dsig_sum / ddsig_sum;
  else
    *scale_update = BIG_SIGMA_UPDATE;

  return loglik;
}
