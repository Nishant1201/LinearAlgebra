""" Solves [A | b] = x by Gauss Elimination,
    A is the 6x6 Vandermode matrix generated
    from the vector,
        v = [1.0 1.2 1.4 1.6 1.8 2.0]'
    and,
        b = [0 1 0 1 0 1]'
    We also evaluate the accuracy of the solution.
    Note - Vandermode matrices tend to be ill conditioned.
"""

import numpy as np
from GaussElimination import *

def Vandermode(v):
    n = len(v)
    A = np.zeros((n,n))
    for j in range(n):          ## 0 to n-1
        A[:,j] = v**(n-j-1)
    return A

v = np.array([1.0, 1.2, 1.4, 1.6, 1.8, 2.0])
b = np.array([0.0, 1.0, 0.0, 1.0, 0.0, 1.0])
A = Vandermode(v)
A_Orig = A.copy()   # Save original coefficient matrix (Deep Copy)
b_Orig = b.copy()   # and the constant vector
x = GaussElimination(A,b)
det = np.prod(np.diagonal(A))
print("x= \n", x)
print("\ndet = ", det)
print("\nCheck result : [A]{x} - b = \n", np.dot(A_Orig, x) - b_Orig)   # Residue = Ax-b is an estimate of error.
input("\nPress return to exit") 