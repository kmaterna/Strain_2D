
## 2D Strain Rate Calculators
This library contains several methods to compute geodetic strain from a GPS velocity field.  It is intended to be used as a means to compare various strain modeling techniques on the same input data (see https://files.scec.org/s3fs-public/0930_Sandwell_UCERF_strain.pdf for an example showing the variability of strain modeling techniques). Several of the strain methods are borrowed from other authors, github repositories, and papers.  Credit should be given to the original authors accordingly.  

### Requirements:
* Python3, numpy, scipy, matplotlib
* GMT5, PYGMT
* Tectonic_Utils - https://github.com/kmaterna/Tectonic_Utils (can be installed by pip)
* Third-party Matlab or Fortran codes may be required depending on the specific strain technique you select, as detailed further below

To install this library, clone the repo and add the parent directory that contains Strain_2D/ to your $PYTHONPATH. 
 
 ### Usage: 
The program is controlled using a config file that specifies inputs/outputs, general strain options, and any required parameters for various strain rate techniques. You can print a sample config file from the main executable ```strain_rate_compute.py```.  

Input velocities must be in a text file. They must be a space-separated table with format as follows: 
```
# lon(deg) lat(deg) VE(mm) VN(mm) VU(mm) SE(mm) SN(mm) SU(mm) name(optional)
``` 
 
The main executable is strain_rate_compute.py in the Strain_2D/Strain_Tools/bin directory. An example run-string would be: 
```python [path-to-code]/strain_rate_compute.py config.txt```

Output strain components and derived quantities (invariants, eigenvectors) are written as grd files or text files and plotted in GMT.  


### Contributing
If you're using this library and have suggestions, let me know!  I'm happy to work together on the code and its applications. See the section on API below if you'd like to contribute a method. 

### Supported strain methods:  

1.  <ins>delaunay_flat</ins>: Delaunay Triangulation, the simplest method. The equations can be found in many papers, including Cai and Grafarend (2007), Journal of Geodynamics, 43, 214-238. No parameters are required to use this method.  

2.  <ins>delaunay</ins>: a generalization of the Delaunay Triangulation for a spherical earth. This Python implementation is based on Savage et al., JGR October 2001, p.22,005, courtesy of Bill Hammond's matlab implementation. No parameters are required to use this method. 

3.  <ins>visr</ins>: The "VISR" method is a fortran code for the interpolation scheme of Zheng-Kang Shen et al., Strain determination using spatially discrete geodetic data, Bull. Seismol. Soc. Am., 105(4), 2117-2127, doi: 10.1785/0120140247, 2015. http://scec.ess.ucla.edu/~zshen/visr/visr.html.  You can download the source code, which must be compiled and linked on your own system, for example by : 
```gfortran -c voronoi_area_version.f90 ``` / ```gfortran visr.f voronoi_area_version.o -o visr.exe```.
Four additional config parameters are required to use this method. 

4.  <ins>gpsgridder</ins>: based on a thin-sheet elastic interpolation scheme from Sandwell, D. T., and P. Wessel (2016), Interpolation of 2-D vector data using constraints from elasticity, GRL.  The implementation of the code is in GMT. Three additional config parameters are required to use this method. 

5. <ins>huang</ins>: the weighted nearest neighbor algorithm of Mong-Han Huang. Two additional config parameters are required to use this method.

6.  <ins>tape</ins>: a wavelet-based matlab program from Tape, Muse, Simons, Dong, Webb, "Multiscale estimation of GPS velocity fields," Geophysical Journal International, 2009 (https://github.com/carltape/surfacevel2strain). Needs a matlab installation and some manual run steps.
  
### Not included methods:
If you have another strain method that you'd be willing to contribute, I would love to work with you to include it!  More methods results in a more robust estimate of strain rate variability.
I have not included the following techniques for the reasons given:

1.  <ins>Spline</ins>: based on Python's built-in spline interpolation of the velocity field. I have found that it doesn't always produce reasonable results.
2.  <ins>NDInterp</ins>: based on Numpy's linear interpolation of the velocity field. It turns out to be basically the same as Delaunay.  

### Units and Sign Conventions
This library uses the following units and sign conventions: 
* station velocities: units of mm/yr in the internal format
* strain rates (exx, etc): units of nanostrain / year (1e-9 / yr) 
* Rotation: units of (1e-3 radians)/yr, or radians / Ka
* Second invariant: units of (nanostrain/yr)^2. These are very unusual units. 
* Dilatation: units of nanostrain / year. Positive means shortening and negative means extension.
* Max Shear: units of nanostrain / year

### Internal Library API
If you're interested in contributing your own type of strain computation, I am happy to add new methods to the library.  You would need to build a Python 'compute' function that matches the API of the other compute functions in the repository (see models/strain_2d.py for template). 

```
class your_strain_method(strain_2d.Strain_2d):
    def __init__(self, Params):
        # perform various types of parameter-setting and input-verification, 
        # different for each technique
        return

    def compute(self, myVelfield):
        [lons, lats, rot_grd, exx_grd, exy_grd, eyy_grd] = your_strain_method_compute(myVelfield... etc);
        return [lons, lats, rot_grd, exx_grd, exy_grd, eyy_grd];
```

where:
* myVelfield is a 1D list of StationVel objects (very easy to use).
* lons is a 1D array of longitudes, in increasing order 
* lats is a 1D array of latitudes, in increasing order
* rot_grd, exx_grd, exy_grd, eyy_grd are 2D arrays that correspond to lons/lats grids
    * exx/exy/eyy have units of nanostrain
    * rot has units of radians/Ka

### Example:
To reproduce the figure in the README:  
1. In the example/ directory, run: ```../Strain_Tools/bin/strain_rate_compute.py --print_config``` to print an example config file (similar to the one provided here in the repo).
2. Run: ```../Strain_Tools/bin/strain_rate_compute.py example_strain_config.txt``` to compute delaunay strain rate using the parameters in the newly-created config file. 
3. Change the method in the config file to try other methods, and re-run each time for different strain calculations with different parameters.
4. Run: ```../Strain_Tools/bin/strain_rate_comparison.py example_strain_config.txt ``` to compute average strain rate maps from several results.
5. Run: ```Display_output/comparison_rows_example.sh -125/-121/38/42``` to view a GMT plot with several strain rate calculations in Northern California together. 

![strain](https://github.com/kmaterna/2D_Strain/blob/master/example/Display_output/output_rows.png)

