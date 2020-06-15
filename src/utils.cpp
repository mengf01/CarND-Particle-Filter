#include "utils.h"
#include <cmath>

double multiv_prob(double sig_x, double sig_y, double x_obs, double y_obs,
                   double mu_x, double mu_y) {
  // calculate normalization term
  double gauss_norm;
  gauss_norm = 1 / (2 * M_PI * sig_x * sig_y);

  // calculate exponent
  double exponent;
  exponent = (pow(x_obs - mu_x, 2) / (2 * pow(sig_x, 2)))
               + (pow(y_obs - mu_y, 2) / (2 * pow(sig_y, 2)));
    
  // calculate weight using normalization terms and exponent
  double weight;
  weight = gauss_norm * exp(-exponent);
    
  return weight;
}

LandmarkObs homogenous_trans(const LandmarkObs& obs, double car_to_map_rot, double x, double y){
    // x and y are in car coordinate.
    LandmarkObs transformed_observation;
    // transform to map x coordinate
    transformed_observation.x = x + (cos(car_to_map_rot) * obs.x) - (sin(car_to_map_rot) * obs.y);
    // transform to map y coordinate
    transformed_observation.y = y + (sin(car_to_map_rot) * obs.x) + (cos(car_to_map_rot) * obs.y);
    transformed_observation.id = obs.id;
    return transformed_observation;
}
