## module Gauss Seidel
''' x, num_iter, omega = GaussSeidel(update_step, x, tol, tol = 1.0e-09)
    Gauss Seidel method for solving [A]{x} = {b}.
    Coefficient matrix [A] should be sparse. User should supply the
    function update_step(x, omega) that returns the new {x} given the 
    current {x} where 'omega' is the relaxation factor.
'''

import numpy as np

def GaussSeidel(update_step, x, tol=1.0e-09):
    omega = 1.0
    k = 10
    p = 1
    for i in range(1,501):
        x_old = x.copy()
        x = update_step(x, omega)
        dx = np.sqrt(np.dot(x-x_old, x-x_old))
        if dx < tol : return x, i, omega
        # Compute relaxation factor after k+p iterations
        if i == k:
            dx1 = dx
        if i == k + p:
            dx2 = dx
            omega = 2.0 / (1.0 + np.sqrt(1.0 \
                - (dx2/dx1)**(1.0/p)))
    print("Gauss Seidel failed to converge.\n")
