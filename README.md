# CFD_Julia Python Implementation

This repo is the Python implementation of [CFD_Julia](https://github.com/surajp92/CFD_Julia).

The corrigendum of its paper, [CFD Julia: A Learning Module Structuring an Introductory Course on Computational Fluid Dynamics](https://www.mdpi.com/2311-5521/4/3/159), is also provided.

Personally, I am more interested in compressible flow solvers. Currently, the tutorials untill Rusanov solver have been implemented. Incompressible flow solvers may be implemented in the future.

Note:
 - most recent updates can be found in branch/basic
 - C++ version may be implemented later
 - GitHub supports LaTeX mathematical expressions, finally!

## Content
- [01. 1D heat equation: Forward time central space (FTCS) scheme](./jupyter/01_Heat_Equation_FTCS/ftcs.ipynb)
- [02. 1D heat equation: Third-order Runge-Kutta (RK3) scheme](./jupyter/02_Heat_Equation_RK3/rk3.ipynb)
- [03. 1D heat equation: Crank-Nicolson (CN) scheme](./jupyter/03_Heat_Equation_CN/cn.ipynb)
- [04. 1D heat equation: Implicit compact Pade (ICP) scheme](./jupyter/04_Heat_Equation_ICP/icp.ipynb)
- [05. 1D inviscid Burgers equation: WENO5 with Dirichlet and periodic boundary condition](./jupyter/05_Inviscid_Burgers_WENO/weno5_dirichlet.ipynb)
- [06. 1D inviscid Burgers equation: CRWENO5 with Dirichlet and periodic boundary conditions](./jupyter/06_Inviscid_Burgers_CRWENO/crweno5_dirichlet.ipynb)
- [07. 1D inviscid Burgers equation: Flux-splitting approach with WENO5](./jupyter/07_Inviscid_Burgers_Flux_Splitting/flux_splitting_periodic.ipynb)
- [08. 1D inviscid Burgers equation: Riemann solver approach with WENO5 using Rusanov solver](./jupyter/08_Inviscid_Burgers_Rieman/rusanov_riemann_periodic.ipynb)
- [09. 1D Euler equations: Roe solver, WENO5, RK3 for time integration](./jupyter/09_Euler_1D_Roe/roe.ipynb)
- [10. 1D Euler equations: HLLC solver, WENO5, RK3 for time integration](./jupyter/10_Euler_1D_HLLC/hllc.ipynb)
- [11. 1D Euler equations: Rusanov solver, WENO5, RK3 for time integration](./jupyter/11_Euler_1D_Rusanov/rusanov.ipynb)

## Corrigendum 

### 2.2 Runge-Kutta Numerical Scheme

- The third-order Runge-Kutta scheme is given as 

    $$
    \begin{align*}
    u_i^{(1)} =&\ u_i^{(n)} + \frac{\alpha \Delta t}{\Delta x^2} \left( u_{i+1}^{(n)} - 2u_i^{(n)} + u_{i-1}^{(n)} \right)\\
    u_i^{(2)} =&\ \frac{3}{4}u_i^{(n)} + \frac{1}{4}u_i^{(1)} + \frac{\alpha \Delta t}{4\Delta x^2}\left( u_{i+1}^{(1)} - 2u_i^{(1)} + u_{i-1}^{(1)} \right)\\
    u_i^{(n+1)} =&\ \frac{1}{3}u_i^{(n)} + \frac{2}{3}u_i^{(2)} + \frac{2\alpha \Delta t}{3\Delta x^2}\left( u_{i+1}^{(2)} - 2u_i^{(2)} + u_{i-1}^{(2)} \right)
    \end{align*}
    $$

### 2.3 Crank-Nicolson Scheme
  
- In Listing 3

    ```
    beta: (alpha*dt)/(2*dx*dx)
    ```

### 2.4 Implicit Compact Pade (ICP) Scheme

- Equation 22:

    $$\frac{1}{\alpha}\left( \frac{1}{12}\frac{\partial u}{\partial t}\Bigr|_{i-1} + \frac{10}{12}\frac{\partial u}{\partial t}\Bigr|_{i} + \frac{1}{12}\frac{\partial u}{\partial t}\Bigr|_{i+1} \right) = \frac{u_{i-1} - 2u_i + u_{i+1}}{\Delta x^2}$$

- Equation 23:

    $$\frac{1}{\alpha}\left( \frac{u_{i-1}^{n+1} - u_{i-1}^{n}}{12\Delta t} + 10\frac{u_{i}^{n+1} - u_{i}^{n}}{12\Delta t} + \frac{u_{i+1}^{n+1} - u_{i+1}^{n}}{12\Delta t} \right) = \frac{u_{i-1}^{n+1} - 2u_{i}^{n+1} + u_{i+1}^{n+1} + u_{i-1}^{n} - 2u_{i}^{n} + u_{i+1}^{n}}{2\Delta x^2}$$

- Equation 28:

    $$r_i = \frac{-2}{\alpha\Delta t}\left( u_{i-1}^{(n)} + 10u_{i}^{(n)} + u_{i+1}^{(n)} \right) - \frac{12}{\Delta x^2}\left( u_{i-1}^{(n)} - 2u_{i}^{(n)} + u_{i+1}^{(n)} \right)$$

### 3.1 WENO-5 Scheme

- Equition 40
  
    $$\beta_2 = \frac{13}{12}(u_i-2u_{i+1}+u_{i+2})^2 + \frac{1}{4}(3u_i - 4u_{i+1}+u_{i+2})^2$$

### 4.2 Riemann Solver: Rusanov Scheme

- Equation 59:

    $$f_{1+1/2} = \frac{1}{2}\left( f_{i+1/2}^L + f_{1+1/2}^R \right) - \frac{c_{1+1/2}}{2}\left( u_{1+1/2}^R - u_{i+1/2}^L \right)$$

- Equation 60:

    $$c_{i+1/2} = max(|u_i|, |u_{i+1}|)$$

### 5 One-Dimensional Euler Solver

- In Equation 62, an important step that derives $F$ with respect to $q$ is missed.

    $$F = 
    \begin{Bmatrix}
    q_2 \\
    \frac{q_2^2}{q_1} + p \\
    \frac{q_2 q_3}{q_1} + p \frac{q_2}{q_1}
    \end{Bmatrix}
    $$

    and

    $$p = (\gamma-1)(q_3 - \frac{q_2^2}{2q_1})$$


### 5.2 HLLC Riemann Solver

- Equation 81:

    $$ 
    \begin{align*}
    F^L,\qquad                                                   & if\ S_L \geq 0 \\
    F^R,\qquad                                                   & if\ S_R \leq 0 \\
    \frac{S_*(S_L u_L - F^L) + S_L P_{LR} D_*}{S_L - S_*},\qquad & if\ S_L \leq 0\ and\ S_* \geq 0 \\
    \frac{S_*(S_R u_R - F^R) + S_R P_{LR} D_*}{S_R - S_*},\qquad & if\ S_R \geq 0\ and\ S_* \leq 0  
    \end{align*}
    $$
