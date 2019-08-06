import numpy as np
import scipy.interpolate as interp
import matplotlib.path

x = (1, 3, 5, 7)
y = (2, 4, 6, 8)
vals = (2, 12, 30, 56)
newx = (1, 2, 3, 4, 5, 6)
newy = (1, 2, 3, 4, 5, 6)

nearest = interp.NearestNDInterpolator((x, y), vals)
newvals = nearest(2, 3)

print(newvals)