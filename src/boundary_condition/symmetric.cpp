//
// Created by yiqif on 2022/5/27.
//

#include "symmetric.h"

void SymmetricBoundaryCondition::Apply(double (&buffer)[FI::PN()][GI::TCX()], unsigned int state) {
    (void)state;
    for (unsigned int i = 0; i < FI::PN(); i++) {
        // update left side
        for (unsigned int j = 0; j < GI::GCX(); j++)
            buffer[i][j] = buffer[i][2 * GI::GCX() - 1 - j];
        // update right side
        for (unsigned int j = GI::FHHX(); j < GI::TCX(); j++)
            buffer[i][j] = buffer[i][GI::LIXC() - j + GI::FHHX()];
    }
}