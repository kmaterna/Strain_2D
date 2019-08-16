"""
Driver program for strain calculation
"""

import method_directory

# Options: hammond, delaunay, gpsgridder, visr, spline,.... many others. 
# note: tape runs seperately

# strain_method="spline"
# strain_method="visr"
# strain_method="gpsgridder"
# strain_method="hammond"
# strain_method="delaunay"
strain_method="ND_interp"


method_directory.method_directory(strain_method);