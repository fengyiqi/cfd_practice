#ifndef HLLC_RIEMANN_SOLVER_H
#define HLLC_RIEMANN_SOLVER_H

#include "user_specification.h"
#include "block.h"
#include <cmath>
#include <algorithm>
#include "utils.h"

template<class ReconstructionScheme>
class HllcRiemannSolver
{
    ReconstructionScheme stencil_;
    Block &block_;

public:
    HllcRiemannSolver() = delete;
    HllcRiemannSolver(Block &block) : block_(block) {};
    ~HllcRiemannSolver() = default;
    void ComputeCellFaceFlux()
    {
        double reconstructed_conservatives_l[FI::EN()];
        double reconstructed_conservatives_r[FI::EN()];
        double rhoL, uL, internal_eL, pL, aL;
        double rhoR, uR, internal_eR, pR, aR;
        double SL, SR, SP, PLR;
        double Ds[3] = {0.0, 1.0, 0.0};
        double stencil_array[stencil_.StencilSize()];

        for (unsigned int i = 0; i < GI::ICX() + 1; i++)
        {
            for (auto &e : FI::Equations)
            {
                // copy appropriate number of values into stencil (from left to right)
                for (unsigned int j = 0; j < stencil_.StencilSize(); j++)
                {
                    stencil_array[j] = block_.conservative_buffer_[e][i + stencil_.ReconstructionStart() + j];
                }
                reconstructed_conservatives_l[e] = stencil_.Apply(stencil_array);
                // copy appropriate number of values into stencil (from right to left)
                for (unsigned int j = 0; j < stencil_.StencilSize(); j++)
                {
                    stencil_array[j] = block_.conservative_buffer_[e][i + stencil_.ReconstructionStart() - j + stencil_.StencilSize()];
                }
                reconstructed_conservatives_r[e] = stencil_.Apply(stencil_array);
            }

            rhoL = reconstructed_conservatives_l[EIndex(FI::EquationEnum::Mass)];
            uL = reconstructed_conservatives_l[EIndex(FI::EquationEnum::MomentumX)] / rhoL;
            internal_eL = reconstructed_conservatives_l[EIndex(FI::EquationEnum::Energy)] / rhoL - 0.5 * uL * uL;
            pL = block_.eos_.ComputePressure(rhoL, internal_eL);
            aL = block_.eos_.ComputeSpeedOfSound(rhoL, pL);

            rhoR = reconstructed_conservatives_r[EIndex(FI::EquationEnum::Mass)];
            uR = reconstructed_conservatives_r[EIndex(FI::EquationEnum::MomentumX)] / rhoR;
            internal_eR = reconstructed_conservatives_r[EIndex(FI::EquationEnum::Energy)] / rhoR - 0.5 * uR * uR;
            pR = block_.eos_.ComputePressure(rhoR, internal_eR);
            aR = block_.eos_.ComputeSpeedOfSound(rhoR, pR);

            SL = std::min(uL, uR) - std::max(aL, aR);
            SR = std::min(uL, uR) + std::max(aL, aR);
            SP = (pR - pL + rhoL * uL * (SL - uL) - rhoR * uR * (SR - uR)) / (rhoL * (SL - uL) - rhoR * (SR - uR));
            PLR = 0.5 * (pL + pR + rhoL * (SL - uL) * (SP - uL) + rhoR * (SR - uR) * (SP - uR));
            Ds[2] = SP;
            // there is a more effective implementation for HLLC, see Nico's thesis
            if (SL >= 0)
            {
                for (auto &e : FI::Equations)
                    block_.cell_face_flux_temp_[e][i] = Utils::ConvertConservativesToFlux(reconstructed_conservatives_l, pL, e);
            }
            else if (SR <= 0)
            {
                for (auto &e : FI::Equations)
                    block_.cell_face_flux_temp_[e][i] = Utils::ConvertConservativesToFlux(reconstructed_conservatives_r, pR, e);
            }
            else if (SP >= 0 && SL <= 0)
            {
                for (auto &e : FI::Equations)
                    block_.cell_face_flux_temp_[e][i] = (SP * (SL * reconstructed_conservatives_l[e] - Utils::ConvertConservativesToFlux(reconstructed_conservatives_l, pL, e)) + SL * PLR * Ds[e]) / (SL - SP);
            }
            else if (SP <= 0 && SR >= 0)
            {
                for (auto &e : FI::Equations)
                    block_.cell_face_flux_temp_[e][i] = (SP * (SR * reconstructed_conservatives_r[e] - Utils::ConvertConservativesToFlux(reconstructed_conservatives_r, pR, e)) + SR * PLR * Ds[e]) / (SR - SP);
            }
        }
    }
};

#endif