import numpy as np
import LUpivot as lup

A = np.array([[2.0,-2.0,6.0], [-2.0,4.0,3.0],[-1.0,8.0,4.0]])
b = np.array([16.0,0.0,-1.0])

lu, seq = lup.LUdecomp(A)
x = lup.LUsolve(lu,b,seq)
print("\nSolution with LU pivot.\n")
print("x= \n", x)