#ifndef COMPUTATION_MODULE_H
#define COMPUTATION_MODULE_H

#include "block.h"
#include "time_integration/runge_kutta_3.h"
#include <iostream>

template<typename RiemannSolverType>
class ComputationModule{
    Block& block_;
    RiemannSolverType riemann_solver_;
    RungeKutta3 runge_kutta_3_;
    double start_time_;
    
public:
    ComputationModule() = delete;
    explicit ComputationModule(Block& block):
                    block_(block),  
                    riemann_solver_(block_),
                    start_time_(0.0) {}
    ~ComputationModule() = default;
    void UpdateRightHandSide(){
    riemann_solver_.ComputeCellFaceFlux();
    for (auto& e : FI::Equations)
        for (unsigned int i = 0; i < GI::ICX(); i++)
            block_.rhs_hand_side_[e][i] = - (block_.cell_face_flux_temp_[e][i+1] - block_.cell_face_flux_temp_[e][i]) / block_.cell_size_;
}
    void TimeIntegration() {
    for (unsigned int stage = 0; stage < runge_kutta_3_.GetTotalStages(); stage++){
        runge_kutta_3_.Advance(stage, block_);
        block_.boundary_condition_.Apply(block_.conservative_buffer_next_, FI::States::Conservatives);
    }
}
    void Solve() {
            while (start_time_ < block_.t_end_) {
        UpdateRightHandSide();
        TimeIntegration();
        start_time_ += block_.t_step_;
        for (auto& e : FI::Equations)
            for (unsigned int i = 0; i < GI::TCX(); i++)
                block_.conservative_buffer_[e][i] = block_.conservative_buffer_next_[e][i];
    }
    block_.ConvertConservativeToPrimitiveStates();
    std::cout << "Done" <<std::endl;
}
};

#endif