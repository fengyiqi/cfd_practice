#ifndef FIRST_ORDER_H
#define FIRST_ORDER_H

class FirstOrder {
    static constexpr unsigned int stencil_size_ = 1;
    static constexpr unsigned int reconstruction_start_ = GI::FICX() - 1;
    
public:
    FirstOrder() = default;
    ~FirstOrder() = default;
    static constexpr unsigned int GetStart() { return reconstruction_start_; }
    static constexpr unsigned int GetStencilSize() { return stencil_size_; }
    constexpr double Apply(const double (&array)[stencil_size_]) {
        return array[0];
    }
};
#endif