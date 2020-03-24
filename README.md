# Global Instability Analysis

**Global Instability Analysis** is a powerful tool for predicting the stability of a given state of flow with respect to external perturbations. It has been successful in predicting several flow phenomena, including boundary layer transition. For the purpose of this analysis, the rate of growth of the disturbance is expressed as a function of the existing perturbed state (through Jacobian of Navier Stokes equations). This allows the formulation of an eigenvalue problem, wherein the eigenvalues of the *Flux Jacobian* are calculated. 

![Equation for instability analysis](/Eq2.JPG =370x115)

The sign of the eigenvalues pf this matrix dictates whether the corresponding disturbance vector grows or dies with time.

The solver developed for this purpose (SIM.py) finds the specified number of the rightmost eigenvalues of the input matrix.

Onset of natural convection in open heated cavities: The base solution of this problem is obtained using OpenFOAM and is used to calculate the Flux Jacobian. By varying the Rayleigh number (different base solutions), the critical value is established through change in sign of the eigenvalues (obtained from the eigenvalue solver). This instability denotes the onset of natural convection. The effect of the Prandtl number is also studied to establish a relationship between the critical Rayleigh number and the Prandtl number. The results are verified through DNS and Landau theory.
