#include "initial_condition.h"
#include <iostream>
#include <cmath>

InitialCondition::InitialCondition(const double& x_start, const double& x_end, const double& t_end): 
    start_(x_start), end_(x_end), t_end_(t_end), length_(x_end - x_start), cell_size_(length_ / GI::ICX()) 
    {}

void InitialCondition::DefineCoordinate(double (&x)[GI::TCX()]){
    double x_start_cell_center = cell_size_ * 0.5;
    double x_ghost_cell_start = x_start_cell_center - GI::GCX() * cell_size_;
    for (unsigned int i = 0; i < GI::TCX(); i++)
        x[i] = x_ghost_cell_start + i * cell_size_;
}

void InitialCondition::DefineInitialPrimitiveStates(const double (&x)[GI::TCX()], double (&primitives)[FI::PN()][GI::TCX()]) {
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

    boundary_condition_.ApplySymmetricBondaryCondition(primitives);
}

void InitialCondition::TestInitialCondition(){
    std::cout << "Test InitialCondition ... \n";
    std::cout << "\t(x_start, x_end) = (" << start_ << ", " << end_ << ")\n";
    std::cout << "\tx_length: " << length_ << std::endl;
    std::cout << "\tcell_size: " << cell_size_ << std::endl;

    double x[GI::TCX()];
    double primitives[FI::PN()][GI::TCX()];
    double conservatives[FI::PN()][GI::TCX()];
    DefineCoordinate(x);
    DefineInitialPrimitiveStates(x, primitives);

    std::cout << "\tcoordinate: \n\t(";
    for (unsigned int i = 0; i < GI::GCX() + 1; i++)
        std::cout << x[i] << " ";
    std::cout << "... ";
    for (unsigned int i = GI::LIXC(); i < GI::TCX(); i++)
        std::cout << x[i] << " ";
    std::cout << ")\n";

    std::cout << "Primitive States with Symmetric Boundary ... \n";
    for (unsigned int i = 0; i < FI::PN(); i++){
        std::cout << "\t(";
        for (unsigned int j = 0; j < GI::GCX() + 1; j++)
            std::cout << primitives[i][j] << " ";
        std::cout << "... ";
        for (unsigned int j = GI::LIXC(); j < GI::TCX(); j++)
            std::cout << primitives[i][j] << " ";
        std::cout << ")\n"; 
    }
}