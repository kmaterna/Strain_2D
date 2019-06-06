This library contains several methods to compute geodetic strain from a GPS velocity field.  It is intended to be used as a means to compare various strain modeling techniques (see https://files.scec.org/s3fs-public/0930_Sandwell_UCERF_strain.pdf for an example showing the variety of strain modeling techniques). Several of the strain methods are borrowed from other authors, github repositories, and papers.  Credit should be given to the original authors accordingly.  

Inside this toolbox, we will have strain methods from: 

1.  Delaunay Triangulation: the simplest method. The equations can be found in many papers, including Cai and Grafarend (2007), Journal of Geodynamics, 43, 214-238. 

2.  The "Hammond" method is a generalization of the Delaunay Triangulation for a spherical earth. This Python implementation is based on Savage et al., JGR October 2001, p.22,005, courtesy of Bill Hammond's initial matlab implementation. 

3.  The "VISR" method is a fortran code for the interpolation scheme of Zheng-Kang Shen et al., Strain determination using spatially discrete geodetic data, Bull. Seismol. Soc. Am., 105(4), 2117-2127, doi: 10.1785/0120140247, 2015. http://scec.ess.ucla.edu/~zshen/visr/visr.html.  You can download the source code, which must be compiled and linked, for example by : 
"gfortran -c voronoi_area_version.f90
gfortran visr.f voronoi_area_version.o -o visr.exe"

4.  The "GPS gridder" method is based on a thin-sheet elastic interpolation scheme from Sandwell, D. T., and P. Wessel (2016), Interpolation of 2-D vector data using constraints from elasticity, GRL.  The implementation of the code is in GMT. 

5.  The "Spline" method is based on Python's built-in interpolation. I have found that it doesn't always produce reasonable results.  

6.  The "Tape" method is a wavelet-based matlab program. It is not implemented yet. 

7.  The "nearest neighbor" method is not implemented yet. 

8.  The "VDoHS" method is not ipmlemented yet. 

9.  The "Gaussian Weighted Interpolation" method is not implemented yet. 

10.  The "Kreemer" method is not implemented yet. 


The program is controlled from the 2D_Strain/ directory by calling "python driver.py". The input and output parameters are controlled in common_io_functions.py.  Output strain components are written as text files and plotted in GMT.  

Requirements and Versions: 
You must add the 2D_Strain/Strain_Code/ directory to your python path. You must add the GMT directory to your bash path. 
The code is written in Python3, and requires GMT5.  It may require Matlab or a Fortran compiler depending on the choice of strain technique.

Example Output: 

![strain](https://github.com/kmaterna/2D_Strain/Example_data/strain_example.png)

