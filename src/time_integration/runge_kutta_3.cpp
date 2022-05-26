#include "runge_kutta_3.h"


void RungeKutta3::Advance(const unsigned int& stage, Block& block) {
        for (unsigned int cs = 0; cs < FI::CN(); cs++){
            for (unsigned int i = 0; i < GI::ICX(); i++){
                block.conservative_buffer_next_[cs][i] = coef_[stage][0] * block.conservative_buffer_[cs][i] + 
                                                         coef_[stage][1] * block.conservative_buffer_next_[cs][i] + 
                                                         coef_[stage][2] * block.rhs_hand_side_[cs][i] * block.t_step_;
            }
        }

    }