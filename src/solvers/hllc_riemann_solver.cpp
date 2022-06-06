#include "hllc_riemann_solver.h"
#include <cmath>
#include <algorithm>
#include <iostream>
#include "utils.h"

// HllcRiemannSolver::HllcRiemannSolver() : weno5_(WENO5()) {}

void HllcRiemannSolver::ComputeCellFaceFlux() {
    
    double reconstructed_conservatives_l[FI::CN()];
    double reconstructed_conservatives_r[FI::CN()];
    double rhoL, uL, internal_eL, pL, aL;
    double rhoR, uR, internal_eR, pR, aR;
    double SL, SR, SP, PLR;
    double Ds[3] = {0.0, 1.0, 0.0};

    for (unsigned int i = 0; i < GI::ICX() + 1; i++){
        for (unsigned int cs = 0; cs < FI::CN(); cs++){
            reconstructed_conservatives_l[cs] = weno5_.Apply({
                block_.conservative_buffer_[cs][reconstruction_start_ + i], 
                block_.conservative_buffer_[cs][reconstruction_start_ + i + 1], 
                block_.conservative_buffer_[cs][reconstruction_start_ + i + 2], 
                block_.conservative_buffer_[cs][reconstruction_start_ + i + 3], 
                block_.conservative_buffer_[cs][reconstruction_start_ + i + 4] 
                });
            reconstructed_conservatives_r[cs] = weno5_.Apply({
                block_.conservative_buffer_[cs][reconstruction_start_ + i + 5], 
                block_.conservative_buffer_[cs][reconstruction_start_ + i + 4], 
                block_.conservative_buffer_[cs][reconstruction_start_ + i + 3], 
                block_.conservative_buffer_[cs][reconstruction_start_ + i + 2], 
                block_.conservative_buffer_[cs][reconstruction_start_ + i + 1] 
                });
        }

        rhoL = reconstructed_conservatives_l[ConservativePool::Mass];
        uL = reconstructed_conservatives_l[ConservativePool::MomentumX] / rhoL;
        internal_eL = reconstructed_conservatives_l[ConservativePool::Energy] / rhoL - 0.5 * uL * uL;
        pL = block_.eos_.ComputePressure(rhoL, internal_eL);
        aL = block_.eos_.ComputeSpeedOfSound(rhoL, pL);

        rhoR = reconstructed_conservatives_r[ConservativePool::Mass];
        uR = reconstructed_conservatives_r[ConservativePool::MomentumX] / rhoR;
        internal_eR = reconstructed_conservatives_r[ConservativePool::Energy] / rhoR - 0.5 * uR * uR;
        pR = block_.eos_.ComputePressure(rhoR, internal_eR);
        aR = block_.eos_.ComputeSpeedOfSound(rhoR, pR);
        


        SL = std::min(uL, uR) - std::max(aL, aR);
        SR = std::min(uL, uR) + std::max(aL, aR);
        SP = (pR - pL + rhoL * uL * (SL - uL) - rhoR * uR * (SR - uR)) / (rhoL * (SL - uL) - rhoR * (SR - uR));
        PLR = 0.5 * (pL + pR + rhoL * (SL - uL) * (SP - uL) + rhoR * (SR - uR) * (SP - uR));
        Ds[2] = SP;

        if (SL >= 0) {
            for (unsigned int cs = 0; cs < FI::CN(); cs++)
                block_.cell_face_flux_temp_[cs][i] = ConvertConservativesToFlux(reconstructed_conservatives_l, pL, cs);
        }
        else if (SR <= 0) {
            for (unsigned int cs = 0; cs < FI::CN(); cs++)
                block_.cell_face_flux_temp_[cs][i] = ConvertConservativesToFlux(reconstructed_conservatives_r, pR, cs);
        }
        else if (SP >= 0 && SL <= 0) {
            for (unsigned int cs = 0; cs < FI::CN(); cs++)
                block_.cell_face_flux_temp_[cs][i] = (SP * (SL * reconstructed_conservatives_l[cs] - ConvertConservativesToFlux(reconstructed_conservatives_l, pL, cs)) + SL * PLR * Ds[cs]) / (SL - SP);
        }
        else if (SP <= 0 && SR >= 0) {
            for (unsigned int cs = 0; cs < FI::CN(); cs++)
                block_.cell_face_flux_temp_[cs][i] = (SP * (SR * reconstructed_conservatives_r[cs] - ConvertConservativesToFlux(reconstructed_conservatives_r, pR, cs)) + SR * PLR * Ds[cs]) / (SR - SP);
        }
    }
}

