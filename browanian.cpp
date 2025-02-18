#include <string>
#include <iostream>
#include <vector>
#include <math.h>
#include <random>
#include <fstream>
#include <sstream>
#include <iterator>
#include <chrono>
#include "dataio.h"
#include "input_params.h"
#include "langevin.h"

int main(int argc, char *argv[])
{
    auto t1 = std::chrono::high_resolution_clock::now();

    std::string paramsPath(argv[1]);
    input_params data;

    int64_t nBins = 100;
    int64_t nBinsCorr = 1000;

    std::vector<double> concentration(nBins, 0.0); 
    std::vector<double> velocityVariance(nBins, 0.0);
    std::vector<double> pdf_vel(nBins, 0.0);
    std::vector<double> diffusion(nBins, 0.0);
    std::vector<std::vector<double>> w_autocorrelator(nBins, std::vector<double>(nBinsCorr, 0.0));

    readParams(data, paramsPath);
    numericalProcedure(data, concentration, w_autocorrelator, diffusion);

#ifdef CONCENTRATION
    saveHist(concentration, argv[2]);
#endif // CONCENTRATION

#ifdef VELOCITY_STAT
    saveHist(velocityVariance, argv[3]);
#endif // VELOCITY_STAT

#ifdef PDF_VELOCITY
    saveHist(pdf_vel, argv[4]);
#endif // PDF_VELOCITY

#ifdef DIFFUSION
    saveHist(diffusion, argv[5]);
#endif // DIFFUSION

#ifdef AUTOCORRELATION
    saveHist(w_autocorrelator[3], argv[6]);
#endif // AUTOCORRELATION

    auto t2 = std::chrono::high_resolution_clock::now();
    auto ms_int = std::chrono::duration_cast<std::chrono::seconds>(t2 - t1);
    std::cout << ms_int.count() << "s\n";
    return 0;
}