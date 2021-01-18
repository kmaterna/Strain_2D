
## 2D Strain Rate Calculators
This library contains several methods to compute geodetic strain from a GPS velocity field.  It is intended to be used as a means to compare various strain modeling techniques on the same input data (see https://files.scec.org/s3fs-public/0930_Sandwell_UCERF_strain.pdf for an example showing the variability of strain modeling techniques). Several of the strain methods are borrowed from other authors, github repositories, and papers.  Credit should be given to the original authors accordingly.  

### Requirements:
* Python3, numpy, scipy, matplotlib
* GMT5
* Tectonic_Utils - https://github.com/kmaterna/Tectonic_Utils
* GNSS_Timeseries_Viewers - https://github.com/kmaterna/GNSS_TimeSeries_Viewers (for the moment)
* Third-party Matlab or Fortran codes may be required depending on the specific strain technique you select, as detailed further below

To install this library, clone the repo and add the parent directory that contains Strain_2D/ to your $PYTHONPATH. 
 
 ### Usage: 
The program is controlled using a config file (see example in test/testing_data/) that specifies inputs/outputs, general strain options, and any required parameters for various strain rate techniques. 

Input velocities must be in one of several formats (still under development). I typically use a format similar to the UNAVCO Plate Boundary Observatory velocity fields. 
 
The main executable is strain_driver.py in the Strain_2D/ directory. An example run-string would be: 
```python [path-to-code]/strain_driver.py config.txt```

Output strain components and derived quantities (invariants, eigenvectors) are written as grd files or text files and plotted in GMT.  


### Contributing
If you're using this library and have suggestions, let me know!  I'm happy to work together on the code and its applications. 

### Supported strain methods:  

1.  <ins>delaunay_flat</ins>: Delaunay Triangulation, the simplest method. The equations can be found in many papers, including Cai and Grafarend (2007), Journal of Geodynamics, 43, 214-238. No parameters are required to use this method.  

2.  <ins>delaunay</ins>: a generalization of the Delaunay Triangulation for a spherical earth. This Python implementation is based on Savage et al., JGR October 2001, p.22,005, courtesy of Bill Hammond's matlab implementation. No parameters are required to use this method. 

3.  <ins>visr</ins>: The "VISR" method is a fortran code for the interpolation scheme of Zheng-Kang Shen et al., Strain determination using spatially discrete geodetic data, Bull. Seismol. Soc. Am., 105(4), 2117-2127, doi: 10.1785/0120140247, 2015. http://scec.ess.ucla.edu/~zshen/visr/visr.html.  You can download the source code, which must be compiled and linked on your own system, for example by : 
```gfortran -c voronoi_area_version.f90 ``` / ```gfortran visr.f voronoi_area_version.o -o visr.exe```.
Four additional config parameters are required to use this method. 

4.  <ins>gps_gridder</ins>: based on a thin-sheet elastic interpolation scheme from Sandwell, D. T., and P. Wessel (2016), Interpolation of 2-D vector data using constraints from elasticity, GRL.  The implementation of the code is in GMT. Three additional config parameters are required to use this method. 

5. <ins>huang</ins>: the weighted nearest neighbor algorithm of Mong-Han Huang. Two additional config parameters are required to use this method.

### Pending methods:
1.  <ins>tape</ins>: a wavelet-based matlab program from Tape, Muse, Simons, Dong, Webb, "Multiscale estimation of GPS velocity fields," Geophysical Journal International, 2009 (https://github.com/carltape/compearth). It is not fully integrated yet.
  
### Not included methods:
If you have another strain method that you'd be willing to contribute, I would love to work with you to include it!  More methods results in a more robust estimate of strain rate variability.
I have not included the following techniques for the reasons given:

1.  <ins>Spline</ins>: based on Python's built-in spline interpolation of the velocity field. I have found that it doesn't always produce reasonable results.
2.  <ins>NDInterp</ins>: based on Numpy's linear interpolation of the velocity field. It turns out to be basically the same as Delaunay.  
3.  <ins>Least squares collocation/kriging</ins>: This method for velocity interpolation has been published by various names (El-Fiky et al., 1997; Kato et al., 1998; Wu et al., 2011; Shen et al. 1996 and 2015 are also variants of the method) but all of them are essentially equivalent to the geostatistical kriging method, where spatial covariance is used to determine optimal weights to apply to nearby stations to estimate a value at an unknown location. To compute the spatial covariance, at minimum a correlation length scale and maximum variance must be assumed. A pure geostatistical approach would be to estimate these directly from the data itself, but other approaches have been used, such as cross-validation. 

### Units and Sign Conventions
This library uses the following units and sign conventions: 
* station velocities: units of mm/yr in the internal format
* strain rates (exx, etc): units of nanostrain / year (1e-9 / yr) 
* Rotation: units of (1e-3 radians)/yr, or radians / Ka
* Second invariant: units of (nanostrain/yr)^2. These are very unusual units. 
* Dilatation: units of nanostrain / year. Positive means shortening and negative means extension.
* Max Shear: units of nanostrain / year


### Example Computations: 

![strain](https://github.com/kmaterna/2D_Strain/blob/master/sample_plots/front_page_four_maps.png)

