//
// Created by yiqif on 2022/5/27.
//

#include "symmetric.h"

void SymmetricBoundaryCondition::Apply(double (&buffer)[FI::PN()][GI::TCX()], FI::States state) {
    (void)state;
    for (auto& p : FI::PrimeStates) {
        // update left side
        for (unsigned int j = 0; j < GI::GCX(); j++)
            buffer[p][j] = buffer[p][2 * GI::GCX() - 1 - j];
        // update right side
        for (unsigned int j = GI::FRGX(); j < GI::TCX(); j++)
            buffer[p][j] = buffer[p][GI::LIXC() - j + GI::FRGX()];
    }
}