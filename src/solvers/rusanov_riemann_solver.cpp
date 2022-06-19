#include "rusanov_riemann_solver.h"
#include <cmath>
#include <algorithm>
#include <iostream>


void RusanovRiemannSolver::ComputeCellFaceFlux(const double (&conservatives)[FI::EN()][GI::TCX()], double (&flux)[FI::EN()][GI::ICX() + 1]) {
    
    double reconstructed_conservatives_l[FI::EN()];
    double reconstructed_conservatives_r[FI::EN()];
    double rhoL, uL, eL, pL, hL, internal_eL;
    double rhoR, uR, eR, pR, hR, internal_eR;
    double uu, hh, aa, alpha, ps;
    double Ds[3] = {0.0, 1.0, 0.0};
    double stencil_array[stencil_.GetStencilSize()];

    for (unsigned int i = 0; i < GI::ICX() + 1; i++){
        for (auto& e : FI::Equations){
            // copy appropriate number of values into stencil (from left to right)
            for (unsigned int j = 0; j < stencil_.GetStencilSize(); j++){
                stencil_array[j] = block_.conservative_buffer_[e][stencil_.GetStart() + i + j];
            }
            reconstructed_conservatives_l[e] = stencil_.Apply(stencil_array);
            // copy appropriate number of values into stencil (from right to left)
            for (unsigned int j = 0; j < stencil_.GetStencilSize(); j++){
                stencil_array[j] = block_.conservative_buffer_[e][stencil_.GetStart() + i + stencil_.GetStencilSize() - j];
            }  
            reconstructed_conservatives_r[e] = stencil_.Apply(stencil_array);
        }

        rhoL = reconstructed_conservatives_l[EIndex(FI::EquationEum::Mass)];
        uL = reconstructed_conservatives_l[EIndex(FI::EquationEum::MomentumX)] / rhoL;
        eL = reconstructed_conservatives_l[EIndex(FI::EquationEum::Energy)] / rhoL;
        internal_eL = reconstructed_conservatives_l[EIndex(FI::EquationEum::Energy)] / rhoL - 0.5 * uL * uL;
        pL = block_.eos_.ComputePressure(rhoL, internal_eL);
        hL = eL + pL / rhoL;

        rhoR = reconstructed_conservatives_r[EIndex(FI::EquationEum::Mass)];
        uR = reconstructed_conservatives_r[EIndex(FI::EquationEum::MomentumX)] / rhoR;
        eR = reconstructed_conservatives_r[EIndex(FI::EquationEum::Energy)] / rhoR;
        internal_eR = reconstructed_conservatives_r[EIndex(FI::EquationEum::Energy)] / rhoR - 0.5 * uR * uR;
        pR = block_.eos_.ComputePressure(rhoR, internal_eR);
        hR = eR + pR / rhoR;

        alpha = 1 / (std::sqrt(std::abs(rhoL)) + std::sqrt(std::abs(rhoR)));

        uu = (std::sqrt(std::abs(rhoL))*uL + std::sqrt(std::abs(rhoR))*uR)*alpha;
        hh = (std::sqrt(std::abs(rhoL))*hL + std::sqrt(std::abs(rhoR))*hR)*alpha;
        aa = std::sqrt(std::abs((block_.gamma_ - 1) * (hh - 0.5 * uu*uu)));
        ps = std::abs(aa + uu);

        for (auto& e : FI::Equations){
            flux[e][i] = 0.5 * (ConvertConservativesToFlux(reconstructed_conservatives_l, e) + ConvertConservativesToFlux(reconstructed_conservatives_r, e));
            flux[e][i] -= 0.5 * ps * (reconstructed_conservatives_r[e] - reconstructed_conservatives_l[e]);
        }
    }
}

double RusanovRiemannSolver::ConvertConservativesToFlux(const double (&conservatives)[FI::EN()], const std::underlying_type<FI::EquationEum>::type equation ) {
    if (equation == FI::Equations[EIndex(FI::EquationEum::Mass)]) 
        return conservatives[EIndex(FI::EquationEum::Mass)];
    else if (equation == FI::Equations[EIndex(FI::EquationEum::MomentumX)]) {
        const double pressure = (block_.gamma_ - 1) * (conservatives[EIndex(FI::EquationEum::Energy)] - 0.5 * conservatives[EIndex(FI::EquationEum::MomentumX)] * conservatives[EIndex(FI::EquationEum::MomentumX)] / conservatives[EIndex(FI::EquationEum::Mass)]);
        return (conservatives[EIndex(FI::EquationEum::MomentumX)] * conservatives[EIndex(FI::EquationEum::MomentumX)]) / conservatives[EIndex(FI::EquationEum::Mass)] + pressure;
    }
    else if (equation == FI::Equations[EIndex(FI::EquationEum::Energy)]) {
        const double pressure = (block_.gamma_ - 1) * (conservatives[EIndex(FI::EquationEum::Energy)] - 0.5 * conservatives[EIndex(FI::EquationEum::MomentumX)] * conservatives[EIndex(FI::EquationEum::MomentumX)] / conservatives[EIndex(FI::EquationEum::Mass)]);
        return (conservatives[EIndex(FI::EquationEum::MomentumX)] * conservatives[EIndex(FI::EquationEum::Energy)] + pressure * conservatives[EIndex(FI::EquationEum::MomentumX)]) / conservatives[EIndex(FI::EquationEum::Mass)];
    }
    else
        throw std::invalid_argument("Not a valid equation!");
}