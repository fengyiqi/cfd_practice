#include "weno5.h"
#include <iostream>

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