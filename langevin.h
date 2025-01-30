#ifndef LANGEVIN_H
#define LANGEVIN_H

#include <vector>
#include "input_params.h"

//#define VELOCITY_STAT
#define PDF_VELOCITY
//#define DEBUG
#define AUTOCORRELATION
//#define CONCENTRATION

void normalizePDF(std::vector<double> &pdf);

void numericalProcedure(std::vector<double> &concentration, 
std::vector<double> &velocityVariance, const input_params params, 
std::vector<double> &pdf_vel, std::vector<double> &w_autocorrelator,
std::vector<double> &phi_autocorrelator);

#endif /* LANGEVIN_H */