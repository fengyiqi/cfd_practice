//
// Created by yiqif on 2022/5/27.
//

#include "periodic.h"

void PeriodicBoundaryCondition::Apply(double (&buffer)[FI::PN()][GI::TCX()], FI::States state) {
    (void)state;
    for (unsigned int i = 0; i < FI::PN(); i++) {
        // update left side
        for (unsigned int j = 0; j < GI::GCX(); j++)
            buffer[i][j] = buffer[i][j+GI::ICX()];
        // update right side
        for (unsigned int j = GI::FRGX(); j < GI::TCX(); j++)
            buffer[i][j] = buffer[i][j - GI::ICX()];
    }
}