#include <iostream>
#include <Rcpp.h>
#include <cmath>

typedef enum {
  IREG_DIST_GAUSSIAN = 0,
  IREG_DIST_LOGISTIC,
  IREG_DIST_EXTREME_VALUE,
  IREG_DIST_UNKNOWN
} IREG_DIST;
