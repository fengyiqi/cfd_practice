#ifndef EQUATION_OF_STATE_H
#define EQUATION_OF_STATE_H

class EquationOfState {
protected:
    double gamma_;
    double background_pressure_;
public:
    EquationOfState() = default;
    EquationOfState(const double gamma) : gamma_(gamma), background_pressure_(0) {}
    EquationOfState(const double gamma, const double background_pressure): gamma_(gamma), background_pressure_(background_pressure) {}
    virtual ~EquationOfState() = default;
    virtual double ComputePressure(const double rho, const double e) = 0;
    virtual double ComputeSpeedOfSound(const double rho, const double pressure) = 0;
};


#endif