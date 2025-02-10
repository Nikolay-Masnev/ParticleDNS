#include "langevin.h"
#include <random>
#include <iostream>
#include <chrono>
#include <iomanip>

const double r_pdf = 3;
const double dr_pdf = 0.1;
const double w_max = 1;
const double Re = 2500;
const uint64_t BUFFER_SIZE = 1e7;
const double maxTimeStep = 1e5;
const double twoPi = 2.0 * 3.141592865358979;

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
    //uint64_t dt = uint64_t(maxTimeStep / autocorrelator.size()); // p точек через каждые dt
    const uint64_t dt = 1;
    const uint64_t autocorr_size = uint64_t(autocorrelator.size());
    const uint64_t buff_size =  uint64_t(W_BUFFER.size());

    for(uint64_t k = 0; k < autocorr_size; ++k)
    {
        for(uint64_t n = 0; n <buff_size - k * dt; ++n)
        {
            autocorrelator[k] += W_BUFFER[n] * W_BUFFER[n + k * dt];
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
    //uint64_t dt = uint64_t(maxTimeStep / autocorrelator.size()); // p точек через каждые dt
    const uint64_t dt = 1;
    const uint64_t autocorr_size = uint64_t(autocorrelator.size());
    const uint64_t buff_size =  uint64_t(W_BUFFER.size());

    for(uint64_t k = 0; k < autocorr_size; ++k)
    {
        autocorrelator[k] /=(counter * (buff_size - k * dt));
    }
}

void numericalProcedure(std::vector<double> &concentration, 
std::vector<double> &velocityVariance, 
const input_params params, std::vector<double> &pdf_vel,
std::vector<double> &w_autocorrelator, std::vector<double> &phi_autocorrelator,
std::vector<double> &phi_tau_corr, std::vector<double> &phi_traj)
{
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(-0.5, 0.5);

    double r = params.r_0;
    double dr = 0;
    double phi = 0;
    double w = 0;
    double dphi = 0;
    double dw = 0;
    double kr1 = 0;
    double kr2 = 0;
    double kw1 = 0;
    double kw2 = 0;
    double kphi1 = 0;
    double kphi2 = 0;
    double sqrt12 = sqrt(12);
    double sqrt_dt = 0;
    double dt = 0;
    double r_bin = params.BoxSize / concentration.size();
    double W1 = 0;
    double W2 = 0;
    uint64_t numBins = concentration.size();
    uint64_t ind = 0;
    uint64_t ind_pdf = 0;

    double tau_invert = pow(params.BoxSize, 2) / (pow(params.a,2) * Re);
    double tau = 1/tau_invert;
    double D_sqrt = sqrt(10);
    double Sigma = 10;

    dt = tau/10;
    sqrt_dt = sqrt(dt);

    uint64_t size = uint64_t(velocityVariance.size());
    std::vector<uint64_t> counter(size, 0.0);

#ifdef DEBUG
    std::cout << r << ' ' << sqrt_dt << ' ' << r_bin << ' ' 
    << numBins << ' ' << ind << ' ' << tau_invert << ' ' << tau << ' ' << dt <<  '\n';
#endif

    std::vector<double> W_BUFFER(BUFFER_SIZE, 0.0);
    uint64_t tmp_ind = 0;

    std::vector<double> PHI_BUFFER(BUFFER_SIZE, 0.0);
    uint64_t phi_tmp_ind = 0;   

    std::vector<double> PHI_CORR_BUFFER(BUFFER_SIZE, 0.0);
    uint64_t phi_corr_tmp_ind = 0;  

    int64_t autocorr_counter = 0;
    int64_t phi_autocorr_counter = 0;
    int64_t phi_corr_counter = 0;

    double w_var = 0;
    double w_r_var = 0;
    int64_t w_count = 0;

    for(uint64_t i = 0; i < params.numSteps; ++i)
    {
        W1 = dis(gen);
        W2 = dis(gen);
        
        kr1 = dt * w * sin(phi);

        kw1 = - dt * tau_invert * w + D_sqrt * sqrt_dt * W1 
        * sqrt12 * sqrt(0.1 + pow(r/params.BoxSize, 2));

        kphi1 = Sigma * (0.2 * tanh(0.5 * r) 
        - 0.1 * tanh(0.1 * r)) * pow(sin(phi),2) 
        + sqrt_dt * W2 * sqrt12 * sqrt(0.1 + pow(r/params.BoxSize, 2));

        kr2 = dt * (w + kw1) * sin(phi + kphi1);

        kw2 = - dt * tau_invert * (w + kw1) 
        + D_sqrt * sqrt_dt * W1 * sqrt12 * sqrt(0.1 + pow((r + kr1)/params.BoxSize, 2));

        kphi2 = Sigma * (0.2 * tanh(0.5 * (r + kr1)) 
        - 0.1 * tanh(0.1 * (r + kr1))) * pow(sin(phi + kphi1),2) 
        + sqrt_dt * W2 * sqrt12 * sqrt(0.1 + pow((r + kr1)/params.BoxSize, 2));

        dr = 0.5 * (kr1 + kr2);
        dw = 0.5 * (kw1 + kw2);
        dphi = 0.5 * (kphi1 + kphi2);

        if(r + dr > params.BoxSize)
            r = 2 * params.BoxSize - r - dr;
        else if (r + dr < 0)
            r = - r - dr;
        else
            r = r + dr;

        w = std::abs(w + dw);
        phi = phi + dphi;

#ifdef CONCENTRATION
        ind = std::min(uint64_t(r / r_bin), numBins-1);
        concentration[ind]++;
#endif // CONCENTRATION

#ifdef PDF_VELOCITY
        if (r > r_pdf - dr_pdf && r < r_pdf + dr_pdf)
        {
            ind_pdf = std::min(int64_t(std::abs(w * sin(phi)) * pdf_vel.size()/w_max), int64_t(pdf_vel.size()-1));
            pdf_vel[ind_pdf]++;

            w_var += std::pow(w, 2);
            w_r_var += std::pow(w * sin(phi), 2);
            w_count++;
        }
#endif // PDF_VELOCITY

#ifdef VELOCITY_STAT
        velocityVariance[ind] += pow(w*sin(phi),2);
        counter[ind]++;
#endif // VELOCITY_STAT

#ifdef AUTOCORRELATION
        if (r > r_pdf - dr_pdf && r < r_pdf + dr_pdf)
        {   
            if(tmp_ind >= BUFFER_SIZE-1)
            {
                tmp_ind = 0;
                autocorr_counter++;
                calcAutoCorr(W_BUFFER, w_autocorrelator);
                std::fill(W_BUFFER.begin(), W_BUFFER.end(), 0.0);
            }
            else
            {
                W_BUFFER[tmp_ind] = w * sin(phi);
                tmp_ind++;
            }

            if(phi_tmp_ind >= BUFFER_SIZE-1)
            {
                phi_tmp_ind = 0;
                phi_autocorr_counter++;
                calcAutoCorr(PHI_BUFFER, phi_autocorrelator);
                std::fill(PHI_BUFFER.begin(), PHI_BUFFER.end(), 0.0);
            }
            else
            {
                PHI_BUFFER[phi_tmp_ind] = phi;
                phi_tmp_ind++;
            }
        }
#endif // AUTOCORRELATION

#ifdef ANGLE_CORRELATION
        if (r > r_pdf - dr_pdf && r < r_pdf + dr_pdf)
        {
            if(phi_corr_tmp_ind >= BUFFER_SIZE -1)
            {
                phi_corr_tmp_ind = 0;
                phi_corr_counter++;
                calcCorr(PHI_CORR_BUFFER, phi_tau_corr);
                std::fill(PHI_CORR_BUFFER.begin(), PHI_CORR_BUFFER.end(), 0.0);
            }
            else
            {
                PHI_CORR_BUFFER[phi_corr_tmp_ind] = phi;
                phi_corr_tmp_ind++;
            }
        }
#endif // ANGLE_CORRELATION

#ifdef PHI_TRAJECTORY
        if(i < phi_traj.size())
        {
            phi_traj[i] = phi;
        }
#endif // PHI_TRAJECTORY
    }

#ifdef CONCENTRATION
    normalizePDF(concentration);
#endif // CONCENTRATION

#ifdef VELOCITY_STAT
    for(uint64_t i = 0; i < uint64_t(size); ++i)
    {
        velocityVariance[i] /= (counter[i]+1);
        velocityVariance[i] = sqrt(velocityVariance[i]);
    }
#endif // VELOCITY_STAT

#ifdef PDF_VELOCITY
    normalizePDF(pdf_vel);
    std::cout << "dispersion of w: " << std::sqrt(w_var/w_count) << '\n'; 
    std::cout << "dispersion of w_r: " << std::sqrt(w_r_var/w_count) << '\n'; 
#endif // PDF_VELOCITY

#ifdef AUTOCORRELATION
    normCorr(W_BUFFER, w_autocorrelator, autocorr_counter);
    normCorr(PHI_BUFFER, phi_autocorrelator, phi_autocorr_counter);
    normCorr(PHI_CORR_BUFFER, phi_tau_corr, phi_corr_counter);
#endif // AUTOCORRELATION
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

void numericalProcedure2(const input_params params, std::vector<double> &concentration)
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

    double w_phi_var = 0;
    double w_r_var = 0;
    uint64_t w_count = 0;

    double W1, W2, W3, W4, W5, W6;

    double dr, dw_r, dw_phi;

    int64_t ind = 0;
    int64_t numBins = concentration.size();
    double r_bin = L / numBins;

    for(uint64_t i = 0; i < steps; ++i)
    {
        W1 = dis(gen);
        W2 = dis(gen);
        W3 = dis(gen);
        W4 = dis(gen);
        W5 = dis(gen);
        W6 = dis(gen);

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

        if(r + dr > L)
            r = 2 * L - r - dr;
        else if (r + dr < 0)
            r = - r - dr;
        else
            r = r + dr;

        w_r += 0.5 * (k_wr_1 + k_wr_2);
        w_phi += 0.5 * (k_wphi_1 + k_wphi_2);

#ifdef CONCENTRATION
        ind = std::min(int64_t(r / r_bin), numBins-1);
        concentration[ind]++;
#endif // CONCENTRATION

#ifdef PDF_VELOCITY
        if (r > r_pdf - dr_pdf && r < r_pdf + dr_pdf)
        {
            w_r_var += std::pow(w_r, 2);
            w_phi_var += std::pow(w_phi, 2);
            w_count++;
        }
#endif // PDF_VELOCITY
    }

#ifdef CONCENTRATION
    normalizePDF(concentration);
#endif // CONCENTRATION

#ifdef PDF_VELOCITY
    std::cout << "dispersion of w_r: " << std::sqrt(w_r_var/w_count) << '\n'; 
    std::cout << "dispersion of w_phi: " << std::sqrt(w_phi_var/w_count) << '\n';
    std::cout << "dispersion of w: " << std::sqrt( (w_phi_var + w_r_var)/w_count) << '\n';
#endif // PDF_VELOCITY
}