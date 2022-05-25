import numpy as np

coef = np.array([
    (1, 0, 1),
    (0.75, 0.25, 0.25),
    (1/3, 2/3, 2/3)
])

def rk3(stage, u, rhs, dt, temp=None):
    return coef[stage-1][0] * u + coef[stage-1][1] * temp + coef[stage-1][2] * dt * rhs
    