#ifndef UTILS_H
#define UTILS_H

#include "helper_functions.h"
double multiv_prob(double sig_x, double sig_y, double x_obs, double y_obs,
                   double mu_x, double mu_y);

LandmarkObs homogenous_trans(const LandmarkObs& obs, double theta, double x, double y);

#endif  // UTILS_H