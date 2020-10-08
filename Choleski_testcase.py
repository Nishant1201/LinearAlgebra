import numpy as np
from Choleski import *

A = np.array([[1.44, -0.36, 5.52, 0.0], [-0.36, 10.33, -7.78, 0.0],\
            [5.52, -7.78, 28.4, 9.0], [0.0, 0.0, 9.0, 61.0]])
b = np.array([0.04, -2.15, 0.0, 0.88])
A_Orig = A.copy()
L = Choleski(A)
x = CholeskiSol(L,b)
print("x =", x)
print("\nCheck : A*x =\n", np.dot(A_Orig,x))
input("\nPress return to exit")