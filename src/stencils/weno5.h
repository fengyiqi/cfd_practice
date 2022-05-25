#ifndef WENO5_H
#define WENO5_H

class WENO5 {
    static constexpr double weight_0_ = 0.1;
    static constexpr double weight_1_ = 0.6;
    static constexpr double weight_2_ = 0.3;

    static constexpr double beta_1st_ = 13.0 / 12.0;
    static constexpr double beta_2nd_ = 0.25;

    static constexpr double stencil_coef_00_ = 1.0 / 3.0;
    static constexpr double stencil_coef_01_ = - 7.0 / 6.0;
    static constexpr double stencil_coef_02_ = 11.0 / 6.0;

    static constexpr double stencil_coef_10_ = - 1.0 / 6.0;
    static constexpr double stencil_coef_11_ = 5.0 / 6.0;
    static constexpr double stencil_coef_12_ = 1.0 / 3.0;

    static constexpr double stencil_coef_20_ = 1.0 / 3.0;
    static constexpr double stencil_coef_21_ = 5.0 / 6.0;
    static constexpr double stencil_coef_22_ = - 1.0 / 6.0;

    static constexpr double epsilon_ = 1e-6;

    static constexpr unsigned int stencil_size_ = 5;
public:
    WENO5() = default;
    ~WENO5() = default;

    constexpr double Apply(const double (&array)[stencil_size_]);
    void Test();
};

#endif