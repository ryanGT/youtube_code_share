import matplotlib.pyplot as plt
import numpy as np
from numpy import array

import control

# Example 1
n1 = array([ 1.        ,  0.18948775])
d1 = array([ 1.        ,  0.82020392,  2.55995487,  0.        ])

G1ex = control.TransferFunction(n1,d1)

# Example 2
n2 = array([ 1.        ,  0.49701718,  7.74116892])
d2 = array([  1.        ,  33.78936269,   2.21071938,   0.33958128,   0.        ])

G2ex = control.TransferFunction(n2,d2)

# Example 3
n3 = array([   1.        ,  188.67345402,   47.0211672 ])
d3 = array([  1.        ,  11.4049788 ,  18.77565203])

G3ex = control.TransferFunction(n3,d3)
