import numpy as np
from LUdecomposition import *

A = np.array([[ 3.0, -1.0, 4.0], [-2.0, 0.0, 5.0], [7.0, 2.0, -2.0]])
B = np.array([[6.0, 3.0, 7.0], [-4.0, 2.0, -5.0]])

A = LUdecomposition(A)          # Decompose A
det = np.prod(np.diagonal(A))
print("\nDeterminant = ", det)
print(" ",len(B))
for i in range(len(B)):         
    x = LUSolve(A,B[i])         # one constant vector at a time
    print("x",i+1,"=",x)
input("Press return to exit.") 
