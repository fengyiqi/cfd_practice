#ifndef COMPUTATION_MODULE_H
#define COMPUTATION_MODULE_H

class ComputationModule{
public:
    ComputationModule() = delete;
    void UpdateRightHandSide();
    void TimeIntegration();
    void Solve();
};

#endif