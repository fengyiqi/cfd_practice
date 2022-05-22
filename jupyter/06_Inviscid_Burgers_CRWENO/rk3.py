def rk3(stage, u, rhs, dt, temp=None):
    if stage == 1:
        return u + dt * rhs
    if stage == 2:
        return 0.75 * u + 0.25 * temp + 0.25 * dt * rhs
    if stage == 3:
        return (1/3) * u + (2/3) * temp + (2/3) * dt * rhs
    