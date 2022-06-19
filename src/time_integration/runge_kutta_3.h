#ifndef RUNGE_KUTTA_3_H
#define RUNGE_KUTTA_3_H 

#include <array>
#include "block.h"

class RungeKutta3 {
    unsigned int total_stages_ = 3;
    static constexpr std::array<std::array<double, 3>, 3> coef_= {{
        {1.0,       0.0,       1.0      },
        {0.75,      0.25,      0.25     },
        {1.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0}
    }};

public:
    RungeKutta3() = default;
    ~RungeKutta3() = default;
    void Advance(const unsigned int& stage, Block& block);
    inline unsigned int GetTotalStages() const { return total_stages_; }
};

#endif