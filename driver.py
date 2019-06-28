"""
Driver program for strain calculation
"""

import method_directory

# Options: hammond, delaunay, gpsgridder, visr, spline, tape,.... many others. 
# strain_method="spline"
# strain_method="visr"
# strain_method="gpsgridder"
strain_method="hammond"
# strain_method="delaunay"

method_directory.method_directory(strain_method);