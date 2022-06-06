//
// Created by yiqif on 2022/5/27.
//

#ifndef CASE_PARAMETERS_H
#define CASE_PARAMETERS_H

#include "boundary_condition/symmetric.h"
#include "boundary_condition/fixed_value.h"
#include <cmath>
#include "eos/stiffened_gas.h"

struct CaseSpecification {

    double x_start_ = 0.0;
    double x_end_ = 0.0;
    
    CaseSpecification() = delete;
    CaseSpecification(double x_start, double x_end) : x_start_(x_start), x_end_(x_end) {};
    virtual ~CaseSpecification() = default;

    virtual BoundaryCondition& GetBoundaryCondition() = 0;
    virtual EquationOfState& GetEquationOfState() = 0;
    virtual double GetTStart() const = 0;
    virtual double GetTEnd() const = 0;
    virtual double GetTStep() const = 0;
    virtual double GetGamma() const = 0;

    virtual void DefineInitialPrimitiveStates(const double (&x)[GI::TCX()], double (&primitives)[FI::PN()][GI::TCX()]) = 0;
    virtual void DefineCoordinate(double (&x)[GI::TCX()]) const {
        double cell_size = (x_end_ - x_start_) / GI::ICX();
        double x_start_cell_center = cell_size * 0.5;
        double x_ghost_cell_start = x_start_cell_center - GI::GCX() * cell_size;
        for (unsigned int i = 0; i < GI::TCX(); i++)
            x[i] = x_ghost_cell_start + i * cell_size;
    }
};

struct SodShockTube : public CaseSpecification {
    static constexpr double t_start = 0.0;
    static constexpr double t_end = 0.2;
    static constexpr double t_step = 0.0001;
    static constexpr double x_start = 0.0;
    static constexpr double x_end = 1.0;
    static constexpr double gamma = 1.4;
    SymmetricBoundaryCondition boundary_condition;
    StiffendGas eos;

    SodShockTube() : CaseSpecification(x_start, x_end), eos(gamma) {}
    ~SodShockTube() = default;

    BoundaryCondition& GetBoundaryCondition() { return boundary_condition; };
    EquationOfState& GetEquationOfState() { return eos; }

    double GetTStart() const { return t_start; }
    double GetTEnd() const { return t_end; }
    double GetTStep() const { return t_step; }
    double GetGamma() const { return gamma; }

    void DefineInitialPrimitiveStates(const double (&x)[GI::TCX()], double (&primitives)[FI::PN()][GI::TCX()]) override {
        // define density
        for (unsigned int i = GI::FICX(); i < GI::FHHX(); i++) {
            if (x[i] < 0.5)
                primitives[PrimeStatePool::Density][i] = 1;
            else
                primitives[PrimeStatePool::Density][i] = 0.125;
        }

        // define velocity
        for (unsigned int i = GI::FICX(); i < GI::FHHX(); i++) {
            primitives[PrimeStatePool::VelocityX][i] = 0.0;
        }

        // define pressure
        for (unsigned int i = GI::FICX(); i < GI::FHHX(); i++) {
            if (x[i] < 0.5)
                primitives[PrimeStatePool::Pressure][i] = 1.0;
            else
                primitives[PrimeStatePool::Pressure][i] = 0.1;
        }
    }
};

struct ShuOsher : public CaseSpecification {
    static constexpr double t_start = 0.0;
    static constexpr double t_end = 1.8;
    static constexpr double t_step = 0.0001;
    static constexpr double x_start = 0.0;
    static constexpr double x_end = 10.0;
    static constexpr double gamma = 1.4;
    FixedValueBoundaryCondition boundary_condition = FixedValueBoundaryCondition(
        {3.857143, 2.629369, 10.33333},
        {0.9735296499804454, 0.0, 1.0}
    );
    StiffendGas eos;

    ShuOsher() : CaseSpecification(x_start, x_end), eos(gamma) {}
    ~ShuOsher() = default;

    BoundaryCondition& GetBoundaryCondition() { return boundary_condition; };
    EquationOfState& GetEquationOfState() { return eos; }
    double GetTStart() const { return t_start; }
    double GetTEnd() const { return t_end; }
    double GetTStep() const { return t_step; }
    double GetGamma() const { return gamma; }

    void DefineInitialPrimitiveStates(const double (&x)[GI::TCX()], double (&primitives)[FI::PN()][GI::TCX()]) override {
        // define density
        for (unsigned int i = GI::FICX(); i < GI::FHHX(); i++) {
            if (x[i] < 1.0)
                primitives[PrimeStatePool::Density][i] = 3.857143;
            else
                primitives[PrimeStatePool::Density][i] = 1 + 0.2 * std::sin(5 * (x[i] - 5));
        }

        // define velocity
        for (unsigned int i = GI::FICX(); i < GI::FHHX(); i++) {
            if (x[i] < 1.0)
                primitives[PrimeStatePool::VelocityX][i] = 2.629369;
            else
                primitives[PrimeStatePool::VelocityX][i] = 0.0;
        }

        // define pressure
        for (unsigned int i = GI::FICX(); i < GI::FHHX(); i++) {
            if (x[i] < 1.0)
                primitives[PrimeStatePool::Pressure][i] = 10.33333;
            else
                primitives[PrimeStatePool::Pressure][i] = 1.0;
        }
    }
};


#endif
