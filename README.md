# CFD_Julia Python Implementation

This repo is the Python implementation of [CFD_Julia](https://github.com/surajp92/CFD_Julia).

The corrigendum of its paper, [CFD Julia: A Learning Module Structuring an Introductory Course on Computational Fluid Dynamics](https://www.mdpi.com/2311-5521/4/3/159), is also provided.

Note:
 - most recent updates can be found in branch/basic
 - C++ version may be implemented later

## Corrigendum 
see corresponding jupyter notebook:
- 2.2 Runge-Kutta Numerical Scheme

    The third-order Runge-Kutta scheme is given as 

$$
\begin{align*}
u_i^{(1)} =&\ u_i^{(n)} + \frac{\alpha \Delta t}{\Delta x^2} \left( u_{i+1}^{(n)} - 2u_i^{(n)} + u_{i-1}^{(n)} \right)\\
u_i^{(2)} =&\ \frac{3}{4}u_i^{(n)} + \frac{1}{4}u_i^{(1)} + \frac{\alpha \Delta t}{4\Delta x^2}\left( u_{i+1}^{(1)} - 2u_i^{(1)} + u_{i-1}^{(1)} \right)\\
u_i^{(n+1)} =&\ \frac{1}{3}u_i^{(n)} + \frac{2}{3}u_i^{(2)} + \frac{2\alpha \Delta t}{3\Delta x^2}\left( u_{i+1}^{(2)} - 2u_i^{(2)} + u_{i-1}^{(2)} \right)
\end{align*}
$$

- [2.3 Crank-Nicolson Scheme](./jupyter/03_Heat_Equation_CN/cn.ipynb)
- [2.4 Implicit Compact Pade (ICP) Scheme](./jupyter/04_Heat_Equation_ICP/icp.ipynb)
- [3.1 WENO-5 Scheme](./jupyter/05_Inviscid_Burgers_WENO/weno5_dirichlet.ipynb)
- [4.2 Riemann Solver: Rusanov Scheme](./jupyter/08_Inviscid_Burgers_Rieman/rusanov_riemann_periodic.ipynb)
- [5   One-Dimensional Euler Solver](./jupyter/09_Euler_1D_Roe/roe.ipynb)
- [5.2 HLLC Riemann Solver](./jupyter/10_Euler_1D_HLLC/hllc.ipynb)