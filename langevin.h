#ifndef LANGEVIN_H
#define LANGEVIN_H

#include <vector>
#include "input_params.h"

#define CONCENTRATION
#define PDF_VELOCITY
#define AUTOCORRELATION
#define DIFFUSION
#define FLUID_CORRELATION

void numericalProcedure(const input_params params, std::vector<double> &concentration, 
std::vector<std::vector<double>> & w_autocorr, std::vector<double> & diffusion);

double Sigma(double r);
double M(double r, double L);
double D(double r, double L);

void calcAutoCorr(std::vector<double> &W_BUFFER, std::vector<double> &autocorrelator);
void calcCorr(std::vector<double> &BUFFER, std::vector<double> &correlator);
void normCorr(std::vector<double> &W_BUFFER, std::vector<double> &autocorrelator, uint64_t counter);
void normalizePDF(std::vector<double> &pdf);
double Integrate(const std::vector<double> &func, double dt);

template <typename T>
void printVector(std::vector<T> &v);

#endif /* LANGEVIN_H */