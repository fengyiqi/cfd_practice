#include "computation_module.h"
#include <iostream>
#include <fstream>

ComputationModule::ComputationModule(Block& block, HllcRiemannSolver& riemann_solver, RungeKutta3& runge_kutta_3) : 
                    block_(block), 
                    runge_kutta_3_(runge_kutta_3), 
                    riemann_solver_(riemann_solver),
                    start_time_(0.0)
{}

void ComputationModule::UpdateRightHandSide() {
    riemann_solver_.ComputeCellFaceFlux(block_.conservative_buffer_, block_.cell_face_flux_temp_);
    for (unsigned int cs = 0; cs < FI::CN(); cs++)
        for (unsigned int i = 0; i < GI::ICX(); i++)
            block_.rhs_hand_side_[cs][i] = - (block_.cell_face_flux_temp_[cs][i+1] - block_.cell_face_flux_temp_[cs][i]) / block_.cell_size_;
}

void ComputationModule::TimeIntegration() {
    for (unsigned int stage = 0; stage < runge_kutta_3_.GetTotalStages(); stage++){
        runge_kutta_3_.Advance(stage, block_);
        block_.boundary_condition_.Apply(block_.conservative_buffer_next_, States::Conservatives);
    }
}

void ComputationModule::Solve() {

    while (start_time_ < block_.t_end_) {
        UpdateRightHandSide();
        TimeIntegration();
        start_time_ += block_.t_step_;
//        std::cout << start_time_ << std::endl;
        for (unsigned int cs = 0; cs < FI::CN(); cs++)
            for (unsigned int i = 0; i < GI::TCX(); i++)
                block_.conservative_buffer_[cs][i] = block_.conservative_buffer_next_[cs][i];
    }
    block_.ConvertConservativeToPrimitiveStates();
    std::cout << "Done" <<std::endl;
}