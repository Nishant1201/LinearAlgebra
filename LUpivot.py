## module LUpivot
''' A,seq = LUdecomp(A,tol=1.0e-9).
    LU decomposition of matrix [A] using scaled row pivoting.
    The returned matrix [A] = contains [U] in the upper
    triangle and the nondiagonal terms of [L] in the lower triangle.
    Note that [L][U] is a row-wise permutation of the original [A];
    the permutations are recorded in the vector {seq}.
    
    x = LUsolve(A,b,seq).
    Solves [L][U]{x} = {b}, where the matrix [A] and the
    permutation vector {seq} are returned from LUdecomp.
'''
import numpy as np
import swap

def LUdecomp(A,tol=1.0e-9):
    n = len(A)
    seq = np.array(range(n))
    
  # Set up scale factors
    s = np.zeros((n))
    for i in range(n):
        s[i] = max(abs(A[i,:]))        
    
    for k in range(0,n-1):
        
      # Row interchange, if needed
        p = np.argmax(np.abs(A[k:n,k])/s[k:n]) + k
        if abs(A[p,k]) <  tol: print('Matrix is singular')
        if p != k:
            swap.swapRows(s,k,p)
            swap.swapRows(A,k,p)
            swap.swapRows(seq,k,p)
            
      # Elimination            
        for i in range(k+1,n):
            if A[i,k] != 0.0:
                lam = A[i,k]/A[k,k]
                A[i,k+1:n] = A[i,k+1:n] - lam*A[k,k+1:n]
                A[i,k] = lam
    return A,seq

def LUsolve(A,b,seq):
    n = len(A)
    
  # Rearrange constant vector; store it in [x]
    x = b.copy()
    for i in range(n):
        x[i] = b[seq[i]]
        
  # Solution
    for k in range(1,n):
        x[k] = x[k] - np.dot(A[k,0:k],x[0:k])
    x[n-1] = x[n-1]/A[n-1,n-1]    
    for k in range(n-2,-1,-1):
       x[k] = (x[k] - np.dot(A[k,k+1:n],x[k+1:n]))/A[k,k]
    return x


