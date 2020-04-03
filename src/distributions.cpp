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
#define BIG_SIGMA_UPDATE 1 /* what to use for (log) scale_update if it should be very very large */


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
 *      mu: grad of LL wrt eta and scale - size should be n_vars + 1
 *      w: vector with diagonal of hessian of LL wrt eta
 *      z: working response; z_i = x_i'beta - mu_i / w_i
 *      scale_update: the Newton update term for scale parameter
 */
double
compute_grad_response(double *w, double *z, double *scale_update, const double *y_l, const double *y_r,
                      const double *eta, const double scale, const IREG_CENSORING *censoring_type,
                      const ull n_obs, const IREG_DIST dist, double *mu, bool debug=false)
{
  double normalized_y[2];     // z^l and z^u, where z^u_i = (y_i - eta_i) / scale
  double densities_l[4];      // F, 1-F, f, f', for the left observation y_l
  double densities_r[4];      // F, 1-F, f, f', for the right observation y_r
  //double mu_i;                // grad of LL wrt eta; mu_i = del g / del eta_i
  double scale_2 = scale * scale, temp, temp2;
  double loglik, dg, ddg, response;
  double dsig, ddsig, dsg, sz;
  double dsig_sum, ddsig_sum;

  switch(dist) {
    case IREG_DIST_EXTREME_VALUE:   sreg_gg = exvalue_d;  break;
    case IREG_DIST_LOGISTIC:        sreg_gg = logistic_d; break;
    case IREG_DIST_GAUSSIAN:        sreg_gg = gauss_d;    break;
    // New code : fixes unsupported switch case errors
    case IREG_DIST_LOG_GAUSSIAN:{
      Rcpp::Function warning("warning");
      warning("compute_grad_response - Unsupported distribution providied: Log Guassian");
      return 0;
    }
    case IREG_DIST_LOG_LOGISTIC:{
      Rcpp::Function warning("warning");
      warning("compute_grad_response - Unsupported distribution providied: Log Logistic");
      return 0;
    }
    case IREG_DIST_EXPONENTIAL:{
      Rcpp::Function warning("warning");
      warning("compute_grad_response - Unsupported distribution providied: exponential");
      return 0;
    }
    case IREG_DIST_WEIBULL:{
      Rcpp::Function warning("warning");
      warning("compute_grad_response - Unsupported distribution providied: Weibull");
      return 0;
    }
    case IREG_DIST_UNKNOWN:{
      Rcpp::Function warning("warning");
      warning("compute_grad_response - Unkown distibution provided");
      return 0;
    }
  }

  loglik = 0;
  dsig_sum = ddsig_sum = 0;
  /* We are skipping pointer checks to save computation time, assume valid ptrs are given */
  for (ull i = 0; i < n_obs; ++i) {

    switch(censoring_type[i]) {
      case IREG_CENSOR_NONE:
        normalized_y[0] = (y_l[i] - eta[i]) / scale;
        sz = scale * normalized_y[0];
        (*sreg_gg)(normalized_y[0], densities_l, 1);    // gives 0, f, f'/f, f''/f

        if (densities_l[1] <= 0) {
          /* off the probability scale -- avoid log(0), and set the
          **  derivatives to gaussian limits (almost any deriv will
          **  do, since the function value triggers step-halving).
          */
          loglik += SMALL;

          dg = -normalized_y[0] / scale;
          ddg = -1 / scale;
          if (scale_update) {
            dsig = ddsig = 0;
          }

        } else {
          loglik += log(densities_l[1]) - log(scale);

          temp = densities_l[2] / scale;
          temp2 = (densities_l[3]) / scale_2;
          dg = -temp; // mu_i = -(1/sigma) * (f'(z) / f(z))
          ddg =  temp2 - dg * dg;

          if (scale_update) {
            dsig = -temp * sz;
            ddsig = sz * sz * temp2 - dsig * (1 + dsig);
            // dsg = sz * temp2 - dg*(dsig +1);
            dsig -= 1;
          }
        }

        break;

      case IREG_CENSOR_RIGHT:
        normalized_y[0] = (y_l[i] - eta[i]) / scale;
        sz = scale * normalized_y[0];
        (*sreg_gg)(normalized_y[0], densities_l, 2);    // gives F, 1-F, f, f'

        if (densities_l[1] <= 0) {
          loglik += SMALL;

          dg = normalized_y[0]/ scale;
          ddg = 0;
          if (scale_update) {
            dsig = ddsig = dsg = 0;
          }

        } else {
          loglik += log(densities_l[1]);

          temp  = -densities_l[2] / densities_l[1] / scale;
          temp2 = -densities_l[3] / densities_l[1] / scale_2;
          dg = -temp; // dg = -(1/sigma) * (f'(z) / f(z))
          ddg = temp2 - dg * dg;     // f'(z^l) / f()] / [1-F()] ...

          if (scale_update) {
            dsig = -temp * sz;
            ddsig = sz * sz* temp2 - dsig * (1 + dsig);
          }
          // dsg = sz * temp2 - dg * (dsig + 1);
        }
        break;

      case IREG_CENSOR_LEFT:
        normalized_y[1] = (y_l[i] - eta[i]) / scale;    // Note: because we store all vals in left column survival style
        sz = scale * normalized_y[1];
        (*sreg_gg)(normalized_y[1], densities_r, 2);    // gives F, 1-F, f, f'

        if (densities_r[0] <= 0) {
          loglik += SMALL;

          dg = -normalized_y[1] / scale;
          ddg = 0;
          if (scale_update) {
            dsig = ddsig = dsg = 0;
          }

        } else {
          loglik += log(densities_r[0]);

          temp = densities_r[2] / densities_r[0] / scale;
          temp2 = densities_r[3] / densities_r[0] / scale_2;
          dg = -densities_r[2] / densities_r[0] / scale;
          ddg = temp2 - dg * dg;

          if (scale_update) {
            dsig = -temp * sz;
            ddsig = sz * sz * temp2 - dsig * (1 + dsig);
          }
        }

        break;

      case IREG_CENSOR_INTERVAL:
        normalized_y[0] = (y_l[i] - eta[i]) / scale;
        normalized_y[1] = (y_r[i] - eta[i]) / scale;

        (*sreg_gg)(normalized_y[0], densities_l, 2);    // gives F, 1-F, f, f'
        (*sreg_gg)(normalized_y[1], densities_r, 2);    // gives F, 1-F, f, f'

        // temp = F(z^u) - F(z^l)
        if (normalized_y[0] > 0) temp = densities_l[1] - densities_r[1];  // stop roundoff in tails
        else                     temp = densities_r[0] - densities_l[0];

        if (temp <= 0) {
          // off the probability scale -- avoid log(0)
          loglik += SMALL;

          dg = 1;
          ddg = 0;
          if (scale_update) {
            dsig = ddsig = dsg = 0;
          }

        } else {
          loglik += log(temp);

          dg = -(densities_r[2] - densities_l[2]) / (temp * scale); // mu_i = -(1/sigma) * (f'(z) / f(z))
          ddg = (densities_r[3] - densities_l[3]) / (temp * scale_2) - dg * dg;

          if (scale_update) {
            dsig = (normalized_y[0] * densities_l[2] - normalized_y[1] * densities_r[2]) / temp;
            ddsig = ((normalized_y[1] * normalized_y[1] * densities_r[3] -
                     normalized_y[0] * normalized_y[0] * densities_l[3]) / temp) -
                    dsig * (1 + dsig);
            // dsg = ((normalized_y[1] * densities_r[3] -
            //         normalized_y[0] * densities_l[3]) / (temp * scale)) - dg * (dsig + 1);
          }
        }

        break;

      default:
        break;
    }
    if (dsig == 0 || ddg == 0)
      response = eta[i];
    else
      response = eta[i] - dg / ddg;

    if (mu) mu[i] = dg;
    if (w) w[i] = ddg;
    if (z) z[i] = response;
    if (scale_update) {
      dsig_sum += dsig;
      ddsig_sum += ddsig;
    }

    // if (debug) {
    //   std::cerr << "\t\t" << z[i] << "\t" << dg << "\t" << ddg << "\n";
    // }

  } // end for: n_obs

  if (scale_update) {
    if (ddsig_sum != 0)
      *scale_update = -dsig_sum / ddsig_sum;
    else
      *scale_update = BIG_SIGMA_UPDATE;
  }

  if (mu)
    mu[n_obs] = dsig_sum;

  return loglik;
}

static void
logistic_d (double z, double ans[4], int j)
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

static void
gauss_d (double z, double ans[4], int j)
{
  double f;

  f = exp(-z*z/2) /SPI;
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

static void
exvalue_d (double z, double ans[4], int j)
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


//' @title Compute the densities of given vector
//'
//' @description
//' C++ function to compute the density or distibution at each point of a
//' given vector. The supported distributions are \code{gaussian},
//' \code{logistic} and (least) \code{extreme value}.
//' \emph{This is used for testing, and may not be useful otherwise.}
//'
//' @param z Vector of points at which the densities are to be calculated.
//'
//' @param j If \code{1}, the returned vector contains
//' 0, f, f'/f and f''/f, where, f is the
//' density function of the given distribution.
//' If 2, returned vector contains F, 1-F, f and f', where F is the
//' distribution function.
//'
//' @param family The distribution of the data. \code{Guassian}, \code{logistic} and
//' (least) \code{extreme value} distributions are supported.
// [[Rcpp::export]]
Rcpp::NumericVector
compute_densities(Rcpp::NumericVector z, int j, Rcpp::String family)
{
  int n_obs = z.size();
  Rcpp::NumericMatrix ans(n_obs, 4);
  double temp[4];

  switch(get_ireg_dist(family)) {
    case IREG_DIST_EXTREME_VALUE:   sreg_gg = exvalue_d;  break;
    case IREG_DIST_LOGISTIC:        sreg_gg = logistic_d; break;
    case IREG_DIST_GAUSSIAN:        sreg_gg = gauss_d;    break;
    // New code : fixes unsupported switch case errors
    case IREG_DIST_LOG_GAUSSIAN:{
      Rcpp::Function warning("warning");
      warning("Unsupported distribution provided: Log Guassian");
      return ans;
    }
    case IREG_DIST_LOG_LOGISTIC:{
      Rcpp::Function warning("warning");
      warning("Unsupported distribution provided: Log Logistic");
      return ans;
    }
    case IREG_DIST_EXPONENTIAL:{
      Rcpp::Function warning("warning");
      warning("Unsupported distribution provided: exponential");
      return ans;
    }
    case IREG_DIST_WEIBULL:{
      Rcpp::Function warning("warning");
      warning("Unsupported distribution provided: Weibull");
      return ans;
    }
    case IREG_DIST_UNKNOWN:{
      Rcpp::Function warning("warning");
      warning("Unkown distibution provided");
      return ans;
    }
  }

  for (int ii = 0; ii < n_obs; ++ii) {
    sreg_gg(z[ii], temp, j);
    for (int jj = 0; jj < 4; ++jj) {
      ans(ii, jj) = temp[jj];
    }
  }
  return ans;
}

//' @title Compute the gradients of the log Likelihood wrt beta, and scale
//'
//' @description
//' C++ function to compute the gradients of log likelihood with respect to the
//' coefficients \code{beta} and \code{scale}.
//' The supported distributions are \code{gaussian},
//' \code{logistic} and (least) \code{extreme value}.
//' \emph{This is used for testing, and may not be useful otherwise.}
//'
//' @param X Design matrix.
//' @param y Output matrix in 2 column format.
//' @param eta Linear predictors.
//' @param scale Scale.
//' @param family The distribution of the data. \code{Guassian}, \code{logistic} and
//' (least) \code{extreme value} distributions are supported.
// [[Rcpp::export]]
Rcpp::List iregnet_compute_gradients(Rcpp::NumericMatrix X, Rcpp::NumericMatrix y,
                                     Rcpp::NumericVector eta, double scale,
                                     Rcpp::String family)
{
  int n_obs = y.nrow(), n_vars = X.ncol();
  Rcpp::NumericVector mu(n_obs+1), out_gradients(n_vars + 1);     // weights and 1 for scale
  double *gradients = REAL(out_gradients);

  /* Loggaussain = Gaussian with log(y), etc. */
  IREG_DIST dist = get_ireg_dist(family);

  IREG_CENSORING status[y.nrow()];
  get_censoring_types(y, status); // It will modify y as well!

  // Note that compute_grad_response returns the derivatives wrt eta (X^T beta),
  // so we need to adjust
  compute_grad_response(NULL, NULL, NULL, REAL(y), REAL(y) + n_obs, REAL(eta), scale,
                        status, n_obs, dist, REAL(mu));

  for (ull j = 0; j < n_vars; ++j) {
    gradients[j] = 0;
    for (ull i = 0; i < n_obs; ++i) {
      gradients[j] += X(i, j) * mu[i];
     }
     gradients[j] = gradients[j] / X.nrow();
  }

  gradients[n_vars] = mu[n_obs];    // wrt scale

  return Rcpp::List::create(Rcpp::Named("gradients") = out_gradients,
                            Rcpp::Named("mu") = mu);
}
