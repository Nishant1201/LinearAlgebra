## Module TridiagLUdecomp
## LU decomposition for Tridiagonal matrices.
## This module contains the functions TridiagLUdecomp and TrodiagLUsolve
## for the decomposition and solution phases of a tridiagonal matrix.
## In TridiagLUsolve, the vector y writes over the constant vector b
## during forward substitution. Similarly, the solution vector
## x overwrites y in the back substitution process.
## Therefore, b contains the solution upon exit from TridiagLUsolve.

''' c,d,e = TridiagLUdecomp(c,d,e).
    LU decomposition of tridiagonal matrix [A], where {c}, {d}
    and {e} are the diagonals of [A]. On output
    {c},{d} and {e} are the diagonals of the decomposed matrix.

    x = TridiagLUsolve(c,d,e,b).
    Solution of [A]{x} {b}, where {c}, {d} and {e} are the
    vectors returned from TridiagLUdecomp.
'''

def TridiagLUdecomp(c,d,e):
    n = len(d)
    for k in range(1,n):
        lam = c[k-1]/d[k-1]
        d[k] = d[k] - lam*e[k-1]
        c[k-1] = lam
    return c,d,e

def TridiagLUsolve(c,d,e,b):
    n = len(d)
    for k in range(1,n):
        b[k] = b[k] - c[k-1]*b[k-1]
    b[n-1] = b[n-1]/d[n-1]
    for k in range(n-2,-1,-1):
        b[k] = (b[k] -e[k]*b[k+1])/d[k]
    return b
