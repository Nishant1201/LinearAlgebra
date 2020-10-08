import numpy as np
from TridiagLUdecomp import *

d = np.ones((5))*2.0
c = np.ones((4))*(-1.0)
b = np.array([5.0, -5.0, 4.0, -5.0, 5.0])
e = c.copy()
c, d, e = TridiagLUdecomp(c,d,e)
x = TridiagLUsolve(c,d,e,b)
print("\nx = \n",x)