# Global Instability Analysis

## Introduction
**Global Instability Analysis** is a powerful tool for predicting the stability of a given state of flow with respect to external perturbations. It has been successful in predicting several flow phenomena, including boundary layer transition. For the purpose of this analysis, the rate of growth of the disturbance is expressed as a function of the existing perturbed state (through Jacobian of Navier Stokes equations). This allows the formulation of an eigenvalue problem, wherein the eigenvalues of the *Flux Jacobian* are calculated. 

![Equation for instability analysis](/Illustrations/Eq2.JPG)

The sign of the eigenvalues pf this matrix dictates whether the corresponding disturbance vector grows or dies with time.

The solver developed for this purpose [(SIM.py)](https://github.com/Abdul-Hannan-Faruqi/Global-Instability-Analysis/tree/master/Solver) finds the specified number of the rightmost eigenvalues of the input matrix. If all eigenvalues lie in the left half-plane (are negative), the flow is considered to be globally stable. However, if even a single eigenvalue is positive, the flow shows instability in that mode (i.e. when perturbed with the corresponding disturbance vector - *eigenvector*).

## Test Cases
1. Lid-driven cavity (LDC) flow: This is a benchmark problem in CFD and it has been used to validate the solver. The flow is stable below a Reynolds no. (Re) ~ 8000 and turns unstable above that. This has been confirmed using the two cases, Re = 100 and Re = 10,000.

2. Natural convection in open heated cavity: This is the problem under investigation. The fluid inside the bottom heated cavity remains stable and shows steady isotherms below a certain Rayleigh number. Once the Rayleigh number exceeds the critical value, the onset of natural convection occurs and the flow becomes unstable.    
