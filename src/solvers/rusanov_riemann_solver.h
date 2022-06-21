#ifndef RUSANOV_RIEMANN_SOLVER_H
#define RUSANOV_RIEMANN_SOLVER_H

#include "user_specification.h"
#include "block.h"
#include <cmath>
#include <algorithm>
#include "utils.h"

template<class ReconstructionScheme>
class RusanovRiemannSolver
{
    ReconstructionScheme stencil_;
    Block &block_;

public:
    RusanovRiemannSolver() = delete;
    RusanovRiemannSolver(Block &block) : block_(block){};
    ~RusanovRiemannSolver() = default;

    void ComputeCellFaceFlux()
    {
        double reconstructed_conservatives_l[FI::EN()];
        double reconstructed_conservatives_r[FI::EN()];
        double rhoL, uL, eL, pL, hL, internal_eL;
        double rhoR, uR, eR, pR, hR, internal_eR;
        double uu, hh, aa, alpha, ps;
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
            eL = reconstructed_conservatives_l[EIndex(FI::EquationEnum::Energy)] / rhoL - 0.5 * uL * uL;
            pL = block_.eos_.ComputePressure(rhoL, eL);
            hL = reconstructed_conservatives_l[EIndex(FI::EquationEnum::Energy)] / rhoL + pL / rhoL;

            rhoR = reconstructed_conservatives_r[EIndex(FI::EquationEnum::Mass)];
            uR = reconstructed_conservatives_r[EIndex(FI::EquationEnum::MomentumX)] / rhoR;
            eR = reconstructed_conservatives_r[EIndex(FI::EquationEnum::Energy)] / rhoR - 0.5 * uR * uR;
            pR = block_.eos_.ComputePressure(rhoR, eR);
            hR = reconstructed_conservatives_r[EIndex(FI::EquationEnum::Energy)] / rhoR + pR / rhoR;

            alpha = 1 / (std::sqrt(std::abs(rhoL)) + std::sqrt(std::abs(rhoR)));

            uu = (std::sqrt(std::abs(rhoL)) * uL + std::sqrt(std::abs(rhoR)) * uR) * alpha;
            hh = (std::sqrt(std::abs(rhoL)) * hL + std::sqrt(std::abs(rhoR)) * hR) * alpha;
            aa = std::sqrt(std::abs((block_.gamma_ - 1) * (hh - 0.5 * uu * uu)));
            ps = std::abs(aa + uu);

            for (auto &e : FI::Equations)
            {
                block_.cell_face_flux_temp_[e][i] = 0.5 * (ConvertConservativesToFlux(reconstructed_conservatives_l, e) + ConvertConservativesToFlux(reconstructed_conservatives_r, e));
                block_.cell_face_flux_temp_[e][i] -= 0.5 * ps * (reconstructed_conservatives_r[e] - reconstructed_conservatives_l[e]);
            }
        }
    }

    double ConvertConservativesToFlux(const double (&conservatives)[FI::EN()], const std::underlying_type<FI::EquationEnum>::type equation)
    {
        if (equation == FI::Equations[EIndex(FI::EquationEnum::Mass)])
            return conservatives[EIndex(FI::EquationEnum::MomentumX)];
        else if (equation == FI::Equations[EIndex(FI::EquationEnum::MomentumX)])
        {
            const double pressure = (block_.gamma_ - 1) * (conservatives[EIndex(FI::EquationEnum::Energy)] - 0.5 * conservatives[EIndex(FI::EquationEnum::MomentumX)] * conservatives[EIndex(FI::EquationEnum::MomentumX)] / conservatives[EIndex(FI::EquationEnum::Mass)]);
            return (conservatives[EIndex(FI::EquationEnum::MomentumX)] * conservatives[EIndex(FI::EquationEnum::MomentumX)]) / conservatives[EIndex(FI::EquationEnum::Mass)] + pressure;
        }
        else if (equation == FI::Equations[EIndex(FI::EquationEnum::Energy)])
        {
            const double pressure = (block_.gamma_ - 1) * (conservatives[EIndex(FI::EquationEnum::Energy)] - 0.5 * conservatives[EIndex(FI::EquationEnum::MomentumX)] * conservatives[EIndex(FI::EquationEnum::MomentumX)] / conservatives[EIndex(FI::EquationEnum::Mass)]);
            return (conservatives[EIndex(FI::EquationEnum::MomentumX)] * conservatives[EIndex(FI::EquationEnum::Energy)] + pressure * conservatives[EIndex(FI::EquationEnum::MomentumX)]) / conservatives[EIndex(FI::EquationEnum::Mass)];
        }
        else
            throw std::invalid_argument("Not a valid equation!");
    }
};

#endif