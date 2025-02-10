#ifndef LANGEVIN_H
#define LANGEVIN_H

#include <vector>
#include "input_params.h"

#define CONCENTRATION
//#define VELOCITY_STAT
#define PDF_VELOCITY
//#define DEBUG
//#define AUTOCORRELATION
//#define ANGLE_CORRELATION
//#define PHI_TRAJECTORY

void normalizePDF(std::vector<double> &pdf);

void numericalProcedure(std::vector<double> &concentration, 
std::vector<double> &velocityVariance, const input_params params, 
std::vector<double> &pdf_vel, std::vector<double> &w_autocorrelator,
std::vector<double> &phi_autocorrelator, std::vector<double> &phi_tau_corr,
std::vector<double> &phi_traj);

void numericalProcedure2(const input_params params, std::vector<double> &concentration);
double Sigma(double r);
double M(double r, double L);
double D(double r, double L);

void calcAutoCorr(std::vector<double> &W_BUFFER, std::vector<double> &autocorrelator);
void calcCorr(std::vector<double> &BUFFER, std::vector<double> &correlator);
void normCorr(std::vector<double> &W_BUFFER, std::vector<double> &autocorrelator, uint64_t counter);
void normalizePDF(std::vector<double> &concentration);

template <typename T>
void printVector(std::vector<T> &v);

#endif /* LANGEVIN_H */