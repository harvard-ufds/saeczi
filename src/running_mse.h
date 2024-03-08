#ifndef running_mse_H
#define running_mse_H

#include <Rcpp.h> 

class RunningMse {
  public:
    RunningMse();
    void update(const Rcpp::NumericVector& newvalues);
    void reset();
    double get_mse() const;
  private:
    double mse;
};

#endif