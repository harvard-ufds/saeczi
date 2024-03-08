#include <Rcpp.h>
#include "running_mse.h"

RunningMse::RunningMse():
  mse(0) {}

void RunningMse::reset() {
  mse = 0;
}

double RunningMse::get_mse() const {
  return mse;
}

void RunningMse::update(const Rcpp::NumericVector& newvalues) {
  
}
