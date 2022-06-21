#ifndef FIRST_ORDER_H
#define FIRST_ORDER_H

class FirstOrder {

    static constexpr unsigned int stencil_size_ = 1;
    static constexpr unsigned int reconstruction_start_ = GI::FICX() - 1;
    
public:
    FirstOrder() = default;
    ~FirstOrder() = default;
    static constexpr inline unsigned int StencilSize() { return stencil_size_; }
    static constexpr inline unsigned int ReconstructionStart() { return reconstruction_start_; }
    double Apply(const double (&array)[stencil_size_]) {
        return array[0];
    }
};
#endif