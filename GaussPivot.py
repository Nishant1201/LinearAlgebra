## module gaussPivot
''' x = gaussPivot(A,b,tol=1.0e-12).
    Solves [A]{x} = {b} by Gauss elimination with
    scaled row pivoting. Apart from row swapping,
    the elimination and solution phases are identical
    to GaussElimination module.
'''    
import numpy as np
import swap 

def GaussPivot(A,b,tol=1.0e-12):
    n = len(b)
    
  # Set up scale factors
    s = np.zeros(n)
    for i in range(n):
        s[i] = max(np.abs(A[i,:]))
            
    for k in range(0,n-1):
        
      # Row interchange, if needed
        p = np.argmax(np.abs(A[k:n,k])/s[k:n]) + k
        if abs(A[p,k]) < tol: print('Matrix is singular')
        if p != k:
            swap.swapRows(b,k,p)
            swap.swapRows(s,k,p)
            swap.swapRows(A,k,p)
            
      # Elimination
        for i in range(k+1,n):
            if A[i,k] != 0.0:
                lam = A[i,k]/A[k,k]
                A[i,k+1:n] = A[i,k+1:n] - lam*A[k,k+1:n]
                b[i] = b[i] - lam*b[k]
    if abs(A[n-1,n-1]) < tol: print('Matrix is singular')
                   
  # Back substitution
    b[n-1] = b[n-1]/A[n-1,n-1]
    for k in range(n-2,-1,-1):
        b[k] = (b[k] - np.dot(A[k,k+1:n],b[k+1:n]))/A[k,k]
    return b



        
