//
// Created by yiqif on 2022/5/27.
//

#ifndef CASE_PARAMETERS_H
#define CASE_PARAMETERS_H

#include "boundary_condition/symmetric.h"

struct CaseSpecification {
    double t_start = 0.0;
    double t_end = 0.0;
    double t_step = 0.0;
    double x_start = 0.0;
    double x_end = 0.0;
    double gamma = 0.0;

    CaseSpecification() = default;
    virtual ~CaseSpecification() = default;
    virtual void DefineInitialPrimitiveStates(const double (&x)[GI::TCX()], double (&primitives)[FI::PN()][GI::TCX()]) {;};
    virtual void DefineCoordinate(double (&x)[GI::TCX()]) const {
        double cell_size = (x_end - x_start) / GI::ICX();
        double x_start_cell_center = cell_size * 0.5;
        double x_ghost_cell_start = x_start_cell_center - GI::GCX() * cell_size;
        for (unsigned int i = 0; i < GI::TCX(); i++)
            x[i] = x_ghost_cell_start + i * cell_size;
    }
    void SetXStart(double x_s) { x_start = x_s; }
    void SetXEnd(double x_e) { x_end = x_e; }
    void SetTStart(double t_s) { t_start = t_s; }
    void SetTEnd(double t_e) {t_end = t_e; }
    void SetTStep(double t_st) {t_step = t_st; }
    void SetGamma(double ga) { gamma = ga; }
};

struct SodShockTube : public CaseSpecification {
    double t_start = 0.0;
    double t_end = 0.2;
    double t_step = 0.0001;
    double x_start = 0.0;
    double x_end = 1.0;
    double gamma = 1.4;

    SodShockTube() {
        SetXStart(x_start);
        SetXEnd(x_end);
        SetTStart(t_start);
        SetTEnd(t_end);
        SetTStep(t_step);
        SetGamma(gamma);
    }
    ~SodShockTube() override = default;
    void DefineCoordinate(double (&x)[GI::TCX()]) const override {
        CaseSpecification::DefineCoordinate(x);
    }
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

#endif
