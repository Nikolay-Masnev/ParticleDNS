#include "langevin.h"
#include <random>
#include <iostream>
#include <chrono>
#include <iomanip>

const double r_pdf = 40;
const double dr_pdf = 0.1;
const double Re = 2500;
const uint64_t BUFFER_SIZE = 1e7;

void normalizePDF(std::vector<double> &concentration)
{
    double sum = 0;
    const uint64_t size = uint64_t(concentration.size());

    for(uint64_t i = 0; i < size; ++i)
    {
        sum += concentration[i];
    }

    for(uint64_t i = 0; i < size; ++i)
    {
        concentration[i] /= sum;
    }
}

template <typename T>
void printVector(std::vector<T> &v)
{
    const uint64_t size = uint64_t(v.size());

    for(uint64_t i = 0; i < size; ++i)
        std::cout << v[i] << ' ';
    std::cout << '\n';
}

void calcAutoCorr(std::vector<double> &W_BUFFER, std::vector<double> &autocorrelator)
{
    const uint64_t autocorr_size = uint64_t(autocorrelator.size());
    const uint64_t buff_size =  uint64_t(W_BUFFER.size());

    for(uint64_t k = 0; k < autocorr_size; ++k)
    {
        for(uint64_t n = 0; n <buff_size - k; ++n)
        {
            autocorrelator[k] += W_BUFFER[n] * W_BUFFER[n + k];
        }
    }
}

void calcCorr(std::vector<double> &BUFFER, std::vector<double> &correlator)
{
    std::cout << "calcCorr\n";
    const uint64_t dt = 1;
    const uint64_t tau = 10;
    const uint64_t corr_size = uint64_t(correlator.size());
    const uint64_t buff_size =  uint64_t(BUFFER.size());

    for(uint64_t k = 0; k < corr_size; ++k)
    {
        for(uint64_t n = 0; n <buff_size - k * dt - tau * dt; ++n)
        {
            correlator[k] += BUFFER[n] - BUFFER[n + tau] ;
        }
    }
}

void normCorr(std::vector<double> &W_BUFFER, std::vector<double> &autocorrelator, uint64_t counter)
{
    const uint64_t autocorr_size = uint64_t(autocorrelator.size());
    const uint64_t buff_size =  uint64_t(W_BUFFER.size());

    for(uint64_t k = 0; k < autocorr_size; ++k)
    {
        autocorrelator[k] /=(counter * (buff_size - k));
    }
}

double Sigma(double r)
{
    return std::abs(10 * (0.2 * tanh(0.5 * r) - 0.1 * tanh(0.1 * r)));
}

double M(double r, double L)
{
    return 0.1 + pow(r/L, 2);
}

double D(double r, double L)
{
    return 10 * (0.1 + pow(r/L, 2));
}

double Integrate(const std::vector<double> &func, double dt)
{
    double sum = 0;
    uint64_t size = func.size();

    for(uint64_t i = 0; i < size; ++i)
    {
        sum += dt * func[i];
    }

    return sum;
}


void numericalProcedure(const input_params params, std::vector<double> &concentration, std::vector<std::vector<double>> &w_autocorrelator,
    std::vector<double> &diffusion)
{
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(-0.5, 0.5);

    double r = params.r_0;
    double w_r = 0;
    double w_phi = 0;
    double k_r1 = 0;
    double k_r2 = 0;
    double k_wr_1 = 0;
    double k_wr_2 = 0;
    double k_wphi_1 = 0;
    double k_wphi_2 = 0;

    double tau_invert = pow(params.BoxSize, 2) / (pow(params.a,2) * Re);
    double tau = 1/tau_invert;
    double dt = tau/10;
    double sqrt_dt = sqrt(dt);
    double sqrt12 = sqrt(12);
    double dt_tau_invert = dt * tau_invert;

    uint64_t steps = params.numSteps;
    double L = params.BoxSize;

#ifdef PDF_VELOCITY
    double w_phi_var = 0;
    double w_r_var = 0;
    double w_r_avg = 0;
    double w_phi_avg = 0;
    uint64_t w_count = 0;
#endif // PDF_VELOCITY

    int64_t ind = 0;
    int64_t numBins = concentration.size();
    double r_bin = L / numBins;

#ifdef FLUID_CORRELATION
    double W1_old, W2_old, W3_old, W4_old, W5_old, W6_old;
    W1_old = 0;
    W2_old = 0;
    W3_old = 0;
    W4_old = 0;
    W5_old = 0;
    W6_old = 0;
#endif // FLUID_CORRELATION

    double W1, W2, W3, W4, W5, W6;
    W1 = 0;
    W2 = 0;
    W3 = 0;
    W4 = 0;
    W5 = 0;
    W6 = 0;

    double tau_corr = tau * 0.1;
    double rho = exp(-dt/tau_corr);
    double sqrt_one_rho = sqrt(1 - rho * rho);
    
    double dr = 0;

#ifdef AUTOCORRELATION
    uint64_t autocorr_size = w_autocorrelator.size();
    std::vector<std::vector<double>> W_BUFFER(w_autocorrelator.size(), std::vector<double>(BUFFER_SIZE, 0.0));
    std::vector<uint64_t> tmp_ind(autocorr_size, 0);
    std::vector<uint64_t> autocorr_counter(autocorr_size, 0);
    double dl = L / autocorr_size;
    uint64_t autocorr_ind = 0;

    std::cout << dl << '\n';
#endif // AUTOCORRELATION

    for(uint64_t i = 0; i < steps; ++i)
    {
#ifdef FLUID_CORRELATION
        W1 = W1_old * rho + sqrt_one_rho * dis(gen);
        W2 = W2_old * rho + sqrt_one_rho * dis(gen);
        W3 = W3_old * rho + sqrt_one_rho * dis(gen);
        W4 = W4_old * rho + sqrt_one_rho * dis(gen);
        W5 = W5_old * rho + sqrt_one_rho * dis(gen);
        W6 = W6_old * rho + sqrt_one_rho * dis(gen);

        W1_old = W1;
        W2_old = W2;
        W3_old = W3;
        W4_old = W4;
        W5_old = W5;
        W6_old = W6;
#else // FLUID_CORRELATION
        W1 = dis(gen);
        W2 = dis(gen);
        W3 = dis(gen);
        W4 = dis(gen);
        W5 = dis(gen);
        W6 = dis(gen);
#endif // // FLUID_CORRELATION

        k_r1 = dt * w_r;

        k_wr_1 = - dt_tau_invert * w_r + sqrt_dt * sqrt12 * (W1 * sqrt(D(r,L)) + (W2 * w_r + W3 * w_phi) * sqrt(M(r, L)));

        k_wphi_1 = -dt_tau_invert * w_phi + sqrt_dt * sqrt12 * (W4 * sqrt(D(r,L)) + (W5 * w_r + W6 * w_phi) * sqrt(M(r,L))) - dt * w_r * Sigma(r);

        k_r2 = dt * (w_r + k_wr_1);

        k_wr_2 = - dt_tau_invert * (w_r + k_wr_1) + sqrt_dt * sqrt12 * (W1 * sqrt(D(r+k_r1,L))
         + (W2 * (w_r+k_wr_1) + W3 * (w_phi + k_wphi_1)) * sqrt(M(r + k_r1, L)));

        k_wphi_2 = -dt_tau_invert * (w_phi + k_wphi_1) 
        + sqrt_dt * sqrt12 * (W4 * sqrt(D(r+k_r1,L)) 
        + (W5 * (w_r + k_wr_1) + W6 * (w_phi + k_wphi_1)) * sqrt(M(r+k_r1,L))) 
        - dt * (w_r + k_wr_1) * Sigma(r + k_r1);

        dr = 0.5 * (k_r1 + k_r2);
        w_r += 0.5 * (k_wr_1 + k_wr_2);
        w_phi += 0.5 * (k_wphi_1 + k_wphi_2);

        if(r + dr > L)
        {
            r = 2 * L - r - dr;
            w_r *= -1;
        }
        else if (r + dr < 0)
        {
            r = - r - dr;
            w_r *= -1;
        }
        else
            r = r + dr;

#ifdef CONCENTRATION
        ind = std::min(int64_t(r / r_bin), numBins-1);
        concentration[ind]++;
#endif // CONCENTRATION

#ifdef PDF_VELOCITY
        if (r > r_pdf - dr_pdf && r < r_pdf + dr_pdf)
        {
            w_r_var += std::pow(w_r, 2);
            w_phi_var += std::pow(w_phi, 2);
            w_r_avg += w_r;
            w_phi_avg += w_phi;
            w_count++;
        }
#endif // PDF_VELOCITY

#ifdef AUTOCORRELATION
        autocorr_ind = uint64_t(r/dl);
        
        if(tmp_ind[autocorr_ind] >= BUFFER_SIZE - 1)
        {
                tmp_ind[autocorr_ind] = 0;
                autocorr_counter[autocorr_ind]++;
                calcAutoCorr(W_BUFFER[autocorr_ind], w_autocorrelator[autocorr_ind]);
                std::fill(W_BUFFER[autocorr_ind].begin(), W_BUFFER[autocorr_ind].end(), 0.0);
        }
        else
        {
                W_BUFFER[autocorr_ind][tmp_ind[autocorr_ind]] = w_r;
                tmp_ind[autocorr_ind]++;
        }
#endif // AUTOCORRELATION
    }

#ifdef CONCENTRATION
    normalizePDF(concentration);
#endif // CONCENTRATION

#ifdef PDF_VELOCITY
    std::cout << "dispersion of w_r: " << std::sqrt(w_r_var/w_count) << '\n'; 
    std::cout << "dispersion of w_phi: " << std::sqrt(w_phi_var/w_count) << '\n';
    std::cout << "dispersion of w: " << std::sqrt( (w_phi_var + w_r_var)/w_count) << '\n';
    std::cout << "average of w_r: " << w_r_avg/w_count << '\n';
    std::cout << "average of w_phi: " << w_phi_avg/w_count << '\n';
#endif // PDF_VELOCITY

#ifdef AUTOCORRELATION
    std::cout << "autocorr\n";
    for(uint64_t i = 0; i < uint64_t(w_autocorrelator.size()); ++i)
    {
        normCorr(W_BUFFER[i], w_autocorrelator[i], autocorr_counter[i]+1);
    }
#endif // AUTOCORRELATION

#ifdef DIFFUSION
    std::cout << "diffusion\n";
    for(uint64_t i = 0; i < uint64_t(w_autocorrelator.size()); ++i)
    {
        diffusion[i] =  Integrate(w_autocorrelator[i], dt);
    }
#endif // DIFFUSION

}
