//
// Created by yiqif on 2022/5/27.
//

#ifndef CASE_PARAMETERS_H
#define CASE_PARAMETERS_H

#include "boundary_condition/symmetric.h"
#include "boundary_condition/fixed_value.h"
#include <cmath>
#include "eos/stiffened_gas.h"

class CaseSpecification {

    // define the x length
    double x_start_ = 0.0;
    double x_end_ = 0.0;
    // define time constrol
    double t_start_ = 0.0;
    double t_end_ = 0.0;
    double t_step_ = 0.0;
    // define gamma
    double gamma_ = 0.0;
    
public:
    CaseSpecification() = delete;
    CaseSpecification(const double x_start, const double x_end, const double t_start, const double t_end, const double t_step, const double gamma)
                    : x_start_(x_start),    x_end_(x_end),      t_start_(t_start),    t_end_(t_end),      t_step_(t_step),     gamma_(gamma) {};
    virtual ~CaseSpecification() = default;
    // define boundary condition
    virtual BoundaryCondition& GetBoundaryCondition() = 0;
    // define equation of state
    virtual EquationOfState& GetEquationOfState() = 0;
    virtual void DefineInitialPrimitiveStates(const double (&x)[GI::TCX()], double (&primitives)[FI::PN()][GI::TCX()]) const = 0;

    virtual inline double GetXStart() const { return x_start_; }
    virtual inline double GetXEnd() const { return x_end_; }
    virtual inline double GetTStart() const { return t_start_; };
    virtual inline double GetTEnd() const { return t_end_; };
    virtual inline double GetTStep() const { return t_step_; };
    virtual inline double GetGamma() const { return gamma_; };
    virtual void DefineCoordinate(double (&x)[GI::TCX()]) const {
        double cell_size = (x_end_ - x_start_) / GI::ICX();
        double x_start_cell_center = cell_size * 0.5;
        double x_ghost_cell_start = x_start_cell_center - GI::GCX() * cell_size;
        for (unsigned int i = 0; i < GI::TCX(); i++)
            x[i] = x_ghost_cell_start + i * cell_size;
    }
};

class SodShockTube : public CaseSpecification {
    // define the x length
    static constexpr double x_start = 0.0;
    static constexpr double x_end = 1.0;
    // define time constrol
    static constexpr double t_start = 0.0;
    static constexpr double t_end = 0.2;
    static constexpr double t_step = 0.0001;
    // define gamma
    static constexpr double gamma = 1.4;
    SymmetricBoundaryCondition boundary_condition;
    StiffendGas eos = StiffendGas(gamma);
public:
    SodShockTube() : CaseSpecification(x_start, x_end, t_start, t_end, t_step, gamma) {}
    ~SodShockTube() = default;

    BoundaryCondition& GetBoundaryCondition() { return boundary_condition; };
    EquationOfState& GetEquationOfState() { return eos; }

    void DefineInitialPrimitiveStates(const double (&x)[GI::TCX()], double (&primitives)[FI::PN()][GI::TCX()]) const override {
        for (unsigned int i = GI::FICX(); i < GI::FRGX(); i++) {

            // define density
            if (x[i] < 0.5)
                primitives[PIndex(FI::PrimeStateEum::Density)][i] = 1.0;
            else
                primitives[PIndex(FI::PrimeStateEum::Density)][i] = 0.125;

            // define velocity
            primitives[PIndex(FI::PrimeStateEum::VelocityX)][i] = 0.0;

            // define pressure
            if (x[i] < 0.5)
                primitives[PIndex(FI::PrimeStateEum::Pressure)][i] = 1.0;
            else
                primitives[PIndex(FI::PrimeStateEum::Pressure)][i] = 0.1;
        }
    }
};



class ShuOsher : public CaseSpecification {
    // define the x length
    static constexpr double x_start = 0.0;
    static constexpr double x_end = 10.0;
    // define time constrol
    static constexpr double t_start = 0.0;
    static constexpr double t_end = 1.8;
    static constexpr double t_step = 0.0001;
    // define gamma
    static constexpr double gamma = 1.4;
    FixedValueBoundaryCondition boundary_condition = FixedValueBoundaryCondition(
        {3.857143, 2.629369, 10.33333},
        {0.9735296499804454, 0.0, 1.0}
    );
    StiffendGas eos = StiffendGas(gamma);
public:
    ShuOsher() : CaseSpecification(x_start, x_end, t_start, t_end, t_step, gamma) {}
    ~ShuOsher() = default;

    BoundaryCondition& GetBoundaryCondition() { return boundary_condition; };
    EquationOfState& GetEquationOfState() { return eos; }

    void DefineInitialPrimitiveStates(const double (&x)[GI::TCX()], double (&primitives)[FI::PN()][GI::TCX()]) const override {
        for (unsigned int i = GI::FICX(); i < GI::FRGX(); i++) {
            if (x[i] < 1.0){
                primitives[PIndex(FI::PrimeStateEum::Density)][i] = 3.857143;
                primitives[PIndex(FI::PrimeStateEum::VelocityX)][i] = 2.629369;
                primitives[PIndex(FI::PrimeStateEum::Pressure)][i] = 10.33333;
            } 
            else{
                primitives[PIndex(FI::PrimeStateEum::Density)][i] = 1 + 0.2 * std::sin(5 * (x[i] - 5));
                primitives[PIndex(FI::PrimeStateEum::VelocityX)][i] = 0.0;
                primitives[PIndex(FI::PrimeStateEum::Pressure)][i] = 1.0;
            }
        }
    }
};

#endif
