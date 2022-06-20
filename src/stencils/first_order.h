#ifndef FIRST_ORDER_H
#define FIRST_ORDER_H

#include "stencil.h"

class FirstOrder : public Stencil<FirstOrder> {

    friend Stencil;
    static constexpr unsigned int stencil_size_ = 1;
    static constexpr unsigned int reconstruction_start_ = GI::FICX() - 1;
    
public:
    explicit constexpr FirstOrder() = default;
    ~FirstOrder() = default;
private:
    constexpr double ApplyImplementation(const double (&array)[stencil_size_]) const {
        return array[0];
    }
};
#endif