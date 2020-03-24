# SIM

## Introduction (Subspace Iteration Method):
The Subspace Iteration method generates a subspace of the original matrix space and projects the matrix onto it. The algorithm then proceeds with the calculation of the eigenvalues and eigenvectors of the projected matrix, while updating the basis vectors of the subspace. This process continues until the basis vectors correspond to the eigenvetors of the original matrix and the eigenvalues of the projected matrix approach its eigenvalues.

## Usage:
**SIM (n,m)**
1. Save the *M* and *A* matrices (as numoy arrays) in the Matrices folder.
2. Run 
     python SIM.py
3. Enter the inputs on the prompt
     dimension of matrix    ->  n
     dimension of subspace  ->  m
4. The program saves the solution in the *Solution* folder
> eig.npy ->  array of *m* eigenvalues
>
> R.npy   ->  *nxm* array of eigenvectors
     
    
