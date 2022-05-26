#include "hllc_riemann_solver.h"
#include <cmath>
#include <algorithm>
#include <iostream>

HllcRiemannSolver::HllcRiemannSolver() : weno5_(WENO5()) {}

void HllcRiemannSolver::ComputeCellFaceFlux(const double (&conservatives)[FI::CN()][GI::TCX()], double (&flux)[FI::CN()][GI::ICX() + 1]) {
    
    double reconstructed_conservatives_l[FI::CN()];
    double reconstructed_conservatives_r[FI::CN()];
    double rhoL, uL, eL, pL, aL;
    double rhoR, uR, eR, pR, aR;
    double SL, SR, SP, PLR;
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
        aL = std::sqrt(std::abs(gamma_ * pL / rhoL));

        rhoR = reconstructed_conservatives_r[ConservativePool::Mass];
        uR = reconstructed_conservatives_r[ConservativePool::MomentumX] / rhoR;
        eR = reconstructed_conservatives_r[ConservativePool::Energy] / rhoR;
        pR = (gamma_ - 1.0) * (eR * rhoR - 0.5 * rhoR * uR * uR);
        aR = std::sqrt(std::abs(gamma_ * pR / rhoR));
        


        SL = std::min(uL, uR) - std::max(aL, aR);
        SR = std::min(uL, uR) + std::max(aL, aR);
        SP = (pR - pL + rhoL * uL * (SL - uL) - rhoR * uR * (SR - uR)) / (rhoL * (SL - uL) - rhoR * (SR - uR));
        PLR = 0.5 * (pL + pR + rhoL * (SL - uL) * (SP - uL) + rhoR * (SR - uR) * (SP - uR));
        Ds[2] = SP;

        if (SL >= 0) {
            for (unsigned int cs = 0; cs < FI::CN(); cs++)
                flux[cs][i] = ConvertConservativesToFlux(reconstructed_conservatives_l, cs);
        }
        else if (SR <= 0) {
            for (unsigned int cs = 0; cs < FI::CN(); cs++)
                flux[cs][i] = ConvertConservativesToFlux(reconstructed_conservatives_r, cs);
        }
        else if (SP >= 0 && SL <= 0) {
            for (unsigned int cs = 0; cs < FI::CN(); cs++)
                flux[cs][i] = (SP * (SL * reconstructed_conservatives_l[cs] - ConvertConservativesToFlux(reconstructed_conservatives_l, cs)) + SL * PLR * Ds[cs]) / (SL - SP);
        }
        else if (SP <= 0 && SR >= 0) {
            for (unsigned int cs = 0; cs < FI::CN(); cs++)
                flux[cs][i] = (SP * (SR * reconstructed_conservatives_r[cs] - ConvertConservativesToFlux(reconstructed_conservatives_r, cs)) + SR * PLR * Ds[cs]) / (SR - SP);
        }
    }
}

double HllcRiemannSolver::ConvertConservativesToFlux(const double (&conservatives)[FI::CN()], const unsigned int& state) {
    if (state == ConservativePool::Mass) 
        return conservatives[ConservativePool::MomentumX];
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