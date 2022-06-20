#ifndef STENCIL_H
#define STENCIL_H

template<typename DerivedStencil>
class Stencil {
    friend DerivedStencil;
    explicit constexpr Stencil() = default;
    
public:
    ~Stencil() = default;
    static constexpr unsigned int StencilSize() { return DerivedStencil::stencil_size_; }
    static constexpr unsigned int ReconstructionStart() { return DerivedStencil::reconstruction_start_; }

    template<typename S>
    constexpr double Apply(const double (&array)[S::StencilSize()]) const {
        return static_cast<DerivedStencil const&>( *this ).ApplyImplementation(array);
    };
};

#endif