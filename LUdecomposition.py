## Module LUDecomposition
""" A = LUdecomposition(A)
    LUdecomposition: [L][U] = [A]

    x = LUSolve(A,b)
    Solution phase : solves [L][U]{x} = {b}

    This module contains both the decomposition
    and the solution phase.
    The decomposition phase returns the matrix [L\\U].
    In the solution phase, the contents of b are replaced
    by y during forward substitution.
    Similarly, the back substitution overwrites y with the
    solution x.
"""

import numpy as np

def LUdecomposition(A):
    n = len(A)
    for k in range(0,n-1):
        for i in range(k+1,n):
            if A[i,k] != 0.0:
                lam = A[i,k]/A[k,k]
                A[i,k+1:n] = A[i,k+1:n] - lam*A[k,k+1:n]
                A[i,k] = lam
    return A

def LUSolve(A,b):
    n = len(A)
    # Forward Substitution
    for k in range(1,n):
        b[k] = b[k] - np.dot(A[k,0:k],b[0:k])
    # Back Subsitution
    b[n-1] = b[n-1]/A[n-1,n-1]
    for k in range(n-2,-1,-1):
        b[k] = (b[k] - np.dot(A[k,k+1:n],b[k+1:n]))/A[k,k]
    return b