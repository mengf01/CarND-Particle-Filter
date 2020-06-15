/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"
#include "utils.cpp"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "helper_functions.h"

using std::string;
using std::vector;
using std::normal_distribution;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  std::default_random_engine gen;
  normal_distribution<double> dist_x(x, std[0]);
  normal_distribution<double> dist_y(y, std[1]);
  normal_distribution<double> dist_theta(theta, std[2]);
  num_particles = 100;  // TODO: Set the number of particles

  for (int i=0; i<num_particles; i++){
    Particle part;
    part.id = i;
    part.weight = 1.0;
    part.x = dist_x(gen);
    part.y = dist_y(gen);
    part.theta = dist_theta(gen);  
    particles.push_back(part);
  }
  is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */
  std::default_random_engine gen;
  normal_distribution<double> dist_x(0, std_pos[0]);
  normal_distribution<double> dist_y(0, std_pos[1]);
  normal_distribution<double> dist_theta(0, std_pos[2]);
  for (auto& part : particles){
    double theta_new = part.theta + yaw_rate*delta_t;
    part.x += velocity/yaw_rate*(sin(theta_new) - sin(part.theta));
    part.y += velocity/yaw_rate*(cos(part.theta) - cos(theta_new));
    part.theta = theta_new;
    // add random noise
    part.x += dist_x(gen);
    part.y += dist_y(gen);
    part.theta += dist_theta(gen);  
  }

}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations) {
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */
  double inf = std::numeric_limits<double>::infinity();
  for (auto& obs : observations){
    int best_id = 0;
    double best_dis_sq = inf;
    for (const auto& p : predicted){
      double cur_dist_sq = (p.x - obs.x) * (p.x - obs.x) + (p.y - obs.y) * (p.y - obs.y) ;
      if (cur_dist_sq < best_dis_sq){
        best_id = p.id;
        best_dis_sq = cur_dist_sq;
      }
    }
    obs.id = best_id;
  }

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian 
   *   distribution. You can read more about this distribution here: 
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system. 
   *   Your particles are located according to the MAP'S coordinate system. 
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */

  // std::cout << "---------------------\n";
  // std::cout << observations.size() << "\n"; // small number such as 6
  // std::cout << map_landmarks.landmark_list.size() << "\n"; // always 42

  for (auto& part : particles){
    part.weight = 1;
    // frame transformation from vehicle coordinate to map coordinate
  	vector<LandmarkObs> transformed_observations;
    for (const auto& obs : observations){
      transformed_observations.emplace_back(homogenous_trans(obs, part.theta, part.x, part.y));
    }
    // prepare predicted based on lidar range
    vector<LandmarkObs> predicted_landmarks;
    for (const auto &landmark : map_landmarks.landmark_list) {
        if (fabs(landmark.x_f - part.x) <= sensor_range && fabs(landmark.y_f - part.y) <= sensor_range) {
            predicted_landmarks.emplace_back(LandmarkObs{landmark.id_i, landmark.x_f, landmark.y_f});
        }
    }
    // Associate the each trasnformed observation to their cloest landmarks
    dataAssociation(predicted_landmarks, transformed_observations);
    // update weight using the multivariate gussian distribution
    for (const auto &transformed_obs : transformed_observations) {
      for (const auto& landmark : map_landmarks.landmark_list){
        if (transformed_obs.id==landmark.id_i){
          part.weight *= multiv_prob(std_landmark[0], std_landmark[1], transformed_obs.x, transformed_obs.y, landmark.x_f, landmark.y_f); 
//           std::cout << "]]]]]]]]]]]]]]]]]]]]]]\n";
//           std::cout << transformed_obs.x << "\n";
//           std::cout << transformed_obs.y << "\n";
//           std::cout << landmark.x_f << "\n";
//           std::cout << landmark.y_f << "\n";
//           std::cout << part.weight << "\n";
          break;
        }
      }
    }
//     std::cout << "--------------\n";
//     std::cout << part.weight << "\n";
  }
}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */
  vector<Particle> parts_new;
  std::random_device rd;
  std::mt19937 gen(rd());
  // Setup the weights
  std::vector<double> weights_all;
  for(uint i=0; i<particles.size(); ++i) {
//     std::cout << "--------------\n";
//     std::cout << i << " "<< particles[i].weight << "\n";
    weights_all.push_back(particles[i].weight);
  }
  // Create the distribution with weights
  std::discrete_distribution<> d(weights_all.begin(), weights_all.end());
  std::vector<double> p = d.probabilities();
  // Resample particles
  for(uint i=0; i<particles.size(); ++i) {
    int index = d(gen);
//     std::cout << index << "\n";
    parts_new.push_back(particles[index]);
  }
  particles = parts_new;
}

void ParticleFilter::SetAssociations(Particle& particle, 
                                     const vector<int>& associations, 
                                     const vector<double>& sense_x, 
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association, 
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}