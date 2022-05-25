#include "weno5.h"
#include <iostream>

constexpr double WENO5::Apply(const double (&array)[stencil_size_]) {

    const double v0 = array[0];
    const double v1 = array[1];
    const double v2 = array[2];
    const double v3 = array[3];
    const double v4 = array[4];

    const double beta0 = beta_1st_ * (v0 - 2.0 * v1 + v2) * (v0 - 2.0 * v1 + v2) + beta_2nd_ * (v0 - 4.0 * v1 + 3.0 * v2) * (v0 - 4.0 * v1 + 3.0 * v2);
    const double beta1 = beta_1st_ * (v1 - 2.0 * v2 + v3) * (v1 - 2.0 * v2 + v3) + beta_2nd_ * (v1 - v3) * (v1 - v3);
    const double beta2 = beta_1st_ * (v2 - 2.0 * v3 + v4) * (v2 - 2.0 * v3 + v4) + beta_2nd_ * (3.0 * v2 - 4.0 * v3 + v4) * (3.0 * v2 - 4.0 * v3 + v4);

    const double alpha0 = weight_0_ / ((beta0 + epsilon_) * (beta0 + epsilon_));
    const double alpha1 = weight_1_ / ((beta1 + epsilon_) * (beta1 + epsilon_)); 
    const double alpha2 = weight_2_ / ((beta2 + epsilon_) * (beta2 + epsilon_)); 

    const double one_alpha = 1.0 / (alpha0 + alpha1 + alpha2);

    const double w0 = alpha0 * one_alpha;
    const double w1 = alpha1 * one_alpha;
    const double w2 = alpha2 * one_alpha;

    return w0 * (stencil_coef_00_ * v0 + stencil_coef_01_ * v1 + stencil_coef_02_ * v2) + 
           w1 * (stencil_coef_10_ * v1 + stencil_coef_11_ * v2 + stencil_coef_12_ * v3) + 
           w2 * (stencil_coef_20_ * v2 + stencil_coef_21_ * v3 + stencil_coef_22_ * v4);
}

void WENO5::Test(){
    std::cout << "Test WENO5 ...\n";

    std::cout << "\t(";
    double data1[stencil_size_] = {1, 1, 1, 1, 1};
    double reconst1 = Apply(data1);
    for (unsigned int i = 0; i < stencil_size_; i++)
        std::cout << data1[i] << " ";
    std::cout <<")\t\t" << reconst1 << std::endl;

    std::cout << "\t(";
    double data2[stencil_size_] = {1, 1, 1, 0, 0};
    double reconst2 = Apply(data2);
    for (unsigned int i = 0; i < stencil_size_; i++)
        std::cout << data2[i] << " ";
    std::cout <<")\t\t" << reconst2 << std::endl;

    std::cout << "\t(";
    double data3[stencil_size_] = {10, 10, 10, -10, -10};
    double reconst3 = Apply(data3);
    for (unsigned int i = 0; i < stencil_size_; i++)
        std::cout << data3[i] << " ";
    std::cout <<")\t" << reconst3 << std::endl;

    std::cout << "\t(";
    double data4[stencil_size_] = {1, 2, 3, 4, 5};
    double reconst4 = Apply(data4);
    for (unsigned int i = 0; i < stencil_size_; i++)
        std::cout << data4[i] << " ";
    std::cout <<")\t\t" << reconst4 << std::endl;
}