## Module PentadiagLUdecomp
## The function PentadiagLUdecomp decomposes a symmetric, pentadiagonal matrix A of the
## form A = [f\e\d\e\f]. The original vectors d, e, and f are destroyed and replaced by
## the vectors of the decomposed matrix. After decomposition, the solution of Ax = b
## can be obtained by PentadiagLUsolve. During forward substitution, the original b is replaced
## by y. Similarly, y is written over by x in the back substitution phase, so that b contains
## the solution vector upon exit from PentadiagLUsolve.

''' d,e,f = PentadiagLUdecomp(d,e,f).
    LU decomposition of symetric pentadiagonal matrix [a], where
    {f}, {e} and {d} are the diagonals of [a]. On output
    {d},{e} and {f} are the diagonals of the decomposed matrix.
    
    x = PentadiagLUsolve(d,e,f,b).
    Solves [a]{x} = {b}, where {d}, {e} and {f} are the vectors
    returned from PentadiagLUdecomp.
    '''

def PentadiagLUdecomp(d,e,f):
    n = len(d)
    for k in range(n-2):
        lam = e[k]/d[k]                         
        d[k+1] = d[k+1] - lam*e[k]              
        e[k+1] = e[k+1] - lam*f[k]              
        e[k] = lam
        lam = f[k]/d[k]
        d[k+2] = d[k+2] - lam*f[k]              
        f[k] = lam
    lam = e[n-2]/d[n-2]
    d[n-1] = d[n-1] -lam*e[n-2]
    e[n-2] = lam
    return d,e,f

def PentadiagLUsolve(d,e,f,b):
    n = len(d)
    b[1] = b[1] - e[0]*b[0]
    for k in range(2,n):
        b[k] = b[k] - e[k-1]*b[k-1] - f[k-2]*b[k-2]
    b[n-1] = b[n-1]/d[n-1]
    b[n-2] = b[n-2]/d[n-2] - e[n-2]*b[n-1]
    for k in range(n-3,-1,-1):
        b[k] = b[k]/d[k] - e[k]*b[k+1] - f[k]*b[k+2]
    return b