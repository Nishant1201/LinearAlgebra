import numpy as np
import GaussPivot as GP

A = np.array([[2.0,-2.0,6.0], [-2.0,4.0,3.0],[-1.0,8.0,4.0]])
b = np.array([16.0,0.0,-1.0])

x = GP.GaussPivot(A,b)
print("\nSolution with Gauss Pivot.\n")
print("x= \n",x)