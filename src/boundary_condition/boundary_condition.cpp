#include "boundary_condition.h"
#include <iostream>


//void BoundaryCondition::TestBondaryCondition() {
//
//    double buffer[FI::PN()][GI::TCX()];
//    for (unsigned int i = 0; i < FI::PN(); i++)
//        for (unsigned int j = GI::FICX(); j < GI::FHHX(); j++)
//            buffer[i][j] = j - GI::FICX();
//
//    std::cout << "Test PeriodicBondaryCondition ... \n";
//    ApplyPeriodicBondaryCondition(buffer);
//    for (unsigned int i = 0; i < FI::PN(); i++){
//        std::cout << "\t(";
//        for (unsigned int j = 0; j < GI::GCX() + 1; j++)
//            std::cout << buffer[i][j] << " ";
//        std::cout << "... ";
//        for (unsigned int j = GI::LIXC(); j < GI::TCX(); j++)
//            std::cout << buffer[i][j] << " ";
//        std::cout << ")\n";
//
//    }
//
//    std::cout << "Test SymmetricBondaryCondition ... \n";
//    ApplySymmetricBondaryCondition(buffer);
//    for (unsigned int i = 0; i < FI::PN(); i++){
//        std::cout << "\t(";
//        for (unsigned int j = 0; j < GI::GCX() + 1; j++)
//            std::cout << buffer[i][j] << " ";
//        std::cout << "... ";
//        for (unsigned int j = GI::LIXC(); j < GI::TCX(); j++)
//            std::cout << buffer[i][j] << " ";
//        std::cout << ")\n";
//    }
//}