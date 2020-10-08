## Module Choleski
""" L = Choleski(A)
    Choleski decomposition : [L][L]transpose = [A]

    x = CholeskiSol(L,b)
    Solution phase of Choleski's decomposition method.
"""

"""
    The code below is implemented with the following observation:
    A[i,j] appears only in the formula for L[i,j]. 
    Therefore, once L[i,j] has been computed, A[i,j] is no longer needed.
    This makes it possible to write the elements of L over the lower triangular
    portion of A as they are computed. The elements above the
    leading diagonal of A will remain untouched.
    The function Choleski(A) implements Choleski’s decomposition.
    If a negative diagonal term is encountered during decomposition,
    an error message is printed and the program is terminated.
    After the coefficient matrix A has been decomposed, the solution of Ax = b can
    be obtained by the usual forward and back substitution operations. The function
    choleskiSol carries out the solution phase. 
"""

import numpy as np
import math

def Choleski(A):
    n = len(A)
    for k in range(n):
        try:
            A[k,k] = math.sqrt(A[k,k] - np.dot(A[k,0:k],A[k,0:k]))
        except ValueError:
            print('Matrix is not positive definite.')
        for i in range(k+1,n):
            A[i,k] = (A[i,k] - np.dot(A[i,0:k],A[k,0:k]))/A[k,k]
    for k in range(1,n): A[0:k,k] = 0.0
    return A

def CholeskiSol(L,b):
    n = len(b)
    # Solution of [L]{y} = {b}
    for k in range(n):
        b[k] = (b[k]-np.dot(L[k,0:k],b[0:k]))/L[k,k]
    # Solution of [L_transpose]{x} = {y}
    for k in range(n-1,-1,-1):
        b[k] = (b[k]-np.dot(L[k+1:n,k],b[k+1:n]))/L[k,k]
    return b