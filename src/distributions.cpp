/*
 * The log-likelihood and other parameters of the distributions will
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

static void exvalue_d (double z, double ans[4], int j);
static void logistic_d      (double z, double ans[4], int j);
static void gauss_d      (double z, double ans[4], int j);

void (*sreg_gg)(double, double [4], int);

#define SPI 2.506628274631001     /* sqrt(2*pi) */
#define ROOT_2 1.414213562373095
#define SMALL -200   /* what to use for log(f(x)) if f(x) is zero */

/* Function to calculate the log likelihood of the data
 * Inputs:
 *
 * Outputs: log likelihood
 */

/* Function to calculate w and z given value of eta, and output values
 *
 * Inputs:
 *      mu: allocated vector for - grad of LL wrt eta; mu_i = del g / del eta_i
 *      w: allocated vector for - diagonal of hessian of LL wrt eta; w_i = del ^2 g / del eta_i ^ 2
 *      z: allocated vector for - working response; z_i = x_i'beta - mu_i / w_i
 *      y_l: left output (transformed to log scale)
 *      y_r: right output (transformed to log scale)
 *      eta: vector of linear predictors; eta_i = x_i' beta
 *      scale: value of the scale (sigma) of the fit
 *      censoring_type: vector denoting the censoring type of each observation
 *      n_obs: number of observations
 *
 * Returns: None
 *
 * Modifies:
 *      mu: grad of LL wrt eta
 *      w: vector with diagonal of hessian of LL wrt eta
 *      z: working response; z_i = x_i'beta - mu_i / w_i
 */
// TODO: if densities are close to 0! (survival)
void compute_grad_response(double *w, double *z, double *y_l, double *y_r,
                           double *eta, double scale, IREG_CENSORING *censoring_type,
                           ull n_obs, IREG_DIST dist)
{
  double normalized_y[2];     // z^l and z^u, where z^u_i = (y_i - eta_i) / scale
  double densities_l[4];      // F, 1-F, f, f', for the left observation y_l
  double densities_r[4];      // F, 1-F, f, f', for the right observation y_r
  double mu_i;                // grad of LL wrt eta; mu_i = del g / del eta_i
  double scale_2 = scale * scale, temp;

  switch(dist) {
    case IREG_DIST_EXTREME_VALUE:   sreg_gg = exvalue_d;  break;
    case IREG_DIST_LOGISTIC:        sreg_gg = logistic_d; break;
    case IREG_DIST_GAUSSIAN:        sreg_gg = gauss_d;    break;
  }

  /* We are skipping pointer checks to save computation time, assume valid ptrs are given */
  for (ull i = 0; i < n_obs; ++i) {

    switch(censoring_type[i]) {
      case IREG_CENSOR_INTERVAL:
        normalized_y[0] = (y_l[i] - eta[i]) / scale;
        normalized_y[1] = (y_r[i] - eta[i]) / scale;

        (*sreg_gg)(normalized_y[0], densities_l, 2);    // gives F, 1-F, f, f'
        (*sreg_gg)(normalized_y[1], densities_r, 2);    // gives F, 1-F, f, f'

        // temp = F(z^u) - F(z^l)
        if (normalized_y[0] > 0) temp = densities_l[1] - densities_r[1];  // stop roundoff in tails
        else                     temp = densities_r[0] - densities_l[0];

        mu_i = -(densities_r[2] - densities_l[2]) / temp / scale; // mu_i = -(1/sigma) * (f'(z) / f(z))
        w[i] = -(densities_r[3] - densities_l[3]) / temp / scale_2 - mu_i * mu_i;

        break;

      case IREG_CENSOR_NONE:
        normalized_y[0] = (y_l[i] - eta[i]) / scale;
        (*sreg_gg)(normalized_y[0], densities_l, 1);    // gives 0, f, f'/f, f''/f

        mu_i = -(densities_l[2]) / scale; // mu_i = -(1/sigma) * (f'(z) / f(z))
        w[i] = -(densities_l[3]) / scale_2 - mu_i;

        break;

      case IREG_CENSOR_LEFT:
        normalized_y[1] = (y_r[i] - eta[i]) / scale;
        (*sreg_gg)(normalized_y[1], densities_r, 2);    // gives F, 1-F, f, f'

        mu_i = -densities_r[2] / densities_r[0] / scale;
        w[i] = -densities_r[3] / densities_r[0] / scale_2 - mu_i * mu_i;

        break;

      case IREG_CENSOR_RIGHT:
        normalized_y[0] = (y_l[i] - eta[i]) / scale;
        (*sreg_gg)(normalized_y[0], densities_l, 2);    // gives F, 1-F, f, f'

        mu_i = densities_l[2] / densities_l[1] / scale; // mu_i = -(1/sigma) * (f'(z) / f(z))
        w[i] = densities_l[3] / densities_l[1] / scale_2 - mu_i * mu_i;     // f'(z^l) / f()] / [1-F()] ...

        break;

      default:
        break;
    }

    z[i] = eta[i] - mu_i / w[i];
  } // end for: n_obs
}

static void logistic_d (double z, double ans[4], int j)
{
  double w, temp;
  int    sign, ii;

  /*
  ** < From Survival >
  ** The symmetry of the logistic allows me to be careful, and never take
  **  exp(large number).  This routine should be very accurate.
  */

  if (z>0)  {
    w = std::exp(-z);
    sign = -1;
    ii=0;
  }
  else {
    w = std::exp(z);
    sign = 1;
    ii=1;
  }
  temp = 1+w;
  switch(j) {
  case 1:  ans[1] = w/(temp*temp);
    ans[2] = sign*(1-w)/temp;
    ans[3] = (w*w -4*w +1)/(temp*temp);
    break;
  case 2:  ans[1-ii] = w/temp;
    ans[ii]   = 1/temp;
    ans[2] = w/(temp*temp);
    ans[3] = sign*ans[2]*(1-w)/temp;
    break;
  }
}

static void gauss_d (double z, double ans[4], int j)
{
  double f;

  f = std::exp(-z*z/2) /SPI;
  switch(j) {
  case 1: ans[1] =f;
    ans[2] = -z;
    ans[3] = z*z -1;
    break;
  case 2: if (z>0) {
    ans[0] = (1 + erf(z/ROOT_2))/2;
    ans[1] =  erfc(z/ROOT_2) /2;
  }
  else {
    ans[1] = (1 + erf(-z/ROOT_2))/2;
    ans[0] =  erfc(-z/ROOT_2) /2;
  }
  ans[2] = f;
  ans[3] = -z*f;
  break;
  }
}

/*
 ** < From Survival >
 ** In the Gaussian and logistic cases, I could avoid numeric disaster by only
 **   evaluating exp(x) for x<0.  By symmetry, I got what I need for
 **   x >0.  The extreme value dist, howerver, is asymmetric, and I don't yet
 **   see the appropriate numeric tricks.
 ** Perhaps a Taylor series will could be used for large z.
 */

static void exvalue_d (double z, double ans[4], int j)
{
  double temp;
  double w;
  if (z < SMALL) w= std::exp(SMALL);
  else if (-z < SMALL) w = std::exp(-SMALL);  /* stop infinite answers */
else   w = std::exp(z);

temp = std::exp(-w);
switch(j) {
case 1:  ans[1] = w*temp;
  ans[2] = 1-w;
  ans[3] = w*(w-3) +1;
  break;
case 2:  ans[0] = 1-temp;
  ans[1] = temp;
  ans[2] = w*temp;
  ans[3] = w*temp*(1-w);
  break;
}
}

