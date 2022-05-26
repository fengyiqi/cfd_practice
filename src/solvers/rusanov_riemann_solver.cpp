#include "rusanov_riemann_solver.h"
#include <cmath>
#include <algorithm>
#include <iostream>

RusanovRiemannSolver::RusanovRiemannSolver() : weno5_(WENO5()) {}

void RusanovRiemannSolver::ComputeCellFaceFlux(const double (&conservatives)[FI::CN()][GI::TCX()], double (&flux)[FI::CN()][GI::ICX() + 1]) {
    
    double reconstructed_conservatives_l[FI::CN()];
    double reconstructed_conservatives_r[FI::CN()];
    double rhoL, uL, eL, pL, hL;
    double rhoR, uR, eR, pR, hR;
    double uu, hh, aa, alpha, ps;
    double Ds[3] = {0.0, 1.0, 0.0};

    for (unsigned int i = 0; i < GI::ICX() + 1; i++){
        for (unsigned int cs = 0; cs < FI::CN(); cs++){
            reconstructed_conservatives_l[cs] = weno5_.Apply({
                conservatives[cs][reconstruction_start_ + i], 
                conservatives[cs][reconstruction_start_ + i + 1], 
                conservatives[cs][reconstruction_start_ + i + 2], 
                conservatives[cs][reconstruction_start_ + i + 3], 
                conservatives[cs][reconstruction_start_ + i + 4] 
                });
            reconstructed_conservatives_r[cs] = weno5_.Apply({
                conservatives[cs][reconstruction_start_ + i + 5], 
                conservatives[cs][reconstruction_start_ + i + 4], 
                conservatives[cs][reconstruction_start_ + i + 3], 
                conservatives[cs][reconstruction_start_ + i + 2], 
                conservatives[cs][reconstruction_start_ + i + 1] 
                });
        }

        rhoL = reconstructed_conservatives_l[ConservativePool::Mass];
        uL = reconstructed_conservatives_l[ConservativePool::MomentumX] / rhoL;
        eL = reconstructed_conservatives_l[ConservativePool::Energy] / rhoL;
        pL = (gamma_ - 1.0) * (eL * rhoL - 0.5 * rhoL * uL * uL);
        hL = eL + pL / rhoL;

        rhoR = reconstructed_conservatives_r[ConservativePool::Mass];
        uR = reconstructed_conservatives_r[ConservativePool::MomentumX] / rhoR;
        eR = reconstructed_conservatives_r[ConservativePool::Energy] / rhoR;
        pR = (gamma_ - 1.0) * (eR * rhoR - 0.5 * rhoR * uR * uR);
        hR = eR + pR / rhoR;

        alpha = 1 / (std::sqrt(std::abs(rhoL)) + std::sqrt(std::abs(rhoR)));

        uu = (std::sqrt(std::abs(rhoL))*uL + std::sqrt(std::abs(rhoR))*uR)*alpha;
        hh = (std::sqrt(std::abs(rhoL))*hL + std::sqrt(std::abs(rhoR))*hR)*alpha;
        aa = std::sqrt(std::abs((gamma_ - 1) * (hh - 0.5 * uu*uu)));
        ps = std::abs(aa + uu);

        for (unsigned int cs = 0; cs < FI::CN(); cs++){
            flux[cs][i] = 0.5 * (ConvertConservativesToFlux(reconstructed_conservatives_l, cs) + ConvertConservativesToFlux(reconstructed_conservatives_r, cs));
            flux[cs][i] -= 0.5 * ps * (reconstructed_conservatives_r[cs] - reconstructed_conservatives_l[cs]);
        }
    }
}

double RusanovRiemannSolver::ConvertConservativesToFlux(const double (&conservatives)[FI::CN()], const unsigned int& state) {
    if (state == ConservativePool::Mass) 
        return conservatives[ConservativePool::Mass];
    else if (state == ConservativePool::MomentumX) {
        const double pressure = (gamma_ - 1) * (conservatives[ConservativePool::Energy] - 0.5 * conservatives[ConservativePool::MomentumX] * conservatives[ConservativePool::MomentumX] / conservatives[ConservativePool::Mass]);
        return (conservatives[ConservativePool::MomentumX] * conservatives[ConservativePool::MomentumX]) / conservatives[ConservativePool::Mass] + pressure;
    }
    else if (state == ConservativePool::Energy) {
        const double pressure = (gamma_ - 1) * (conservatives[ConservativePool::Energy] - 0.5 * conservatives[ConservativePool::MomentumX] * conservatives[ConservativePool::MomentumX] / conservatives[ConservativePool::Mass]);
        return (conservatives[ConservativePool::MomentumX] * conservatives[ConservativePool::Energy] + pressure * conservatives[ConservativePool::MomentumX]) / conservatives[ConservativePool::Mass];
    }
    else
        return 0.0;
}