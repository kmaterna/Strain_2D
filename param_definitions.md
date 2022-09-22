
# Explanation of config parameters for Strain_2D calculators

### [general]
* ```method```: string, one of the methods below
* ```output_dir```: string, path to output parent-directory
* ```input_vel_file```: string, path to text file wtih input velocities

### [strain]
* ```range_strain```: float/float/float/float, representing the target region for strain rate to be calculated upon, in W/E/S/N degrees longitude and latitude
* ```range_data```: float/float/float/float, representing the region of velocity data that will be used for the calculation, in W/E/S/N degrees longitude and latitude
* ```inc```: float/float, representing x/y grid spacing of resulting strain grid, in degrees longitude and latitude 

### [delaunay]/ [delaunay_flat]
* no parameters

### [visr]
* [See Native Documentation](http://scec.ess.ucla.edu/~zshen/visr/visr.html)
* ```distance_weighting```: string, either 'gaussian' or 'quadratic'. 
  * If Gaussian: Distance weights (Shen et al., 2015) 

$$ L_i = e^{- \Delta R_i^2 / D^2} $$ 

  * If Quadratic: Distance weights (Shen et al., 2015)
  
$$ L_i = \frac {1} {1 + \Delta R_i^2 / D^2} $$ 

  * ```spatial_weighting```: string, either 'voronoi' or 'azimuth'. Corresponds to choosing a *Z_i* function in Shen et al. (2015)
  * ```min_max_inc_smooth```: float/float/float, representing minimum, maximum, and incremental spatial smoothing constants (I think that's *D*) (km).  the 'Dist' in the output will be calculated from one of these choices (min + n*inc). 
The result is less block-like if you use finer increments. 
  * ```weighting_threshold```: float, weighting threshold, *Wt*, from Shen et al., 2015; it ends up determining the optimal smoothing distance *D*. It's very important.
  * ```uncertainty_threshold```: float, cutoff uncertainty for velocity data input. The VISR examples have 0.5. I've been putting it too low (0.05) and I didn't have bad effects... not sure what units it has. 
  * ```executable```: string, path to location of compiled fortran executable, visr.exe or similar

  * Based on the VISR documentation, I'm not sure where Shen et al. (2015)'s *L_0* (distance weighting threshold for ignoring the weight) gets defined. It might not get defined in the config file. 

### [gpsgridder]
* [See Native Documentation](http://gmt.soest.hawaii.edu/doc/latest/supplements/potential/gpsgridder.html) 
* ```poisson```: float, poisson's ratio used in gpsgridder -S argument
* ```fd```: float, fudge factor used in gpsgridder -Fd argument. The GMT default value is 0.01
* ```eigenvalue```: float, ratio of the smallest eigenvalue used in the fit to the largest eigenvalue, placed in gpsgridder -C argument 


### [loc_avg_grad]
* ```EstimateRadiusKm```: float, radius of searching neighborhood in km
* ```nstations```: integer, minimum number of stations within search radius 

### [wavelets]

* [See Native Documentation](https://github.com/carltape/surfacevel2strain/blob/master/USER_INFO/surfacevel2strain_manual.pdf)
* ```code_dir```: string, path to location where surfacevel2strain on your computer system, such as '/Users/usrname/Documents/Software/surfacevel2strain'
* ```qmin```: integer, minimum scale wavelength for the computation 
* ```qmax```: integer, maximum scale wavelength for the computation
* ```qsec```: integer, scale wavelength for the secular velocity field. Must be between qmin and qmin. 
* ```output_tag```: string, path to the directory with preferred Matlab strain results from compearth. This directory is usually created by the Matlab run. 

### [geostats]
* ```model_type```: string, one of [Gaussian, Exponential, Nugget].
* ```sill_east```: float
* ```range_east```: float
* ```nugget_east```: float
* ```sill_north```: float
* ```range_north```: float
* ```nugget_north```: float
* ```trend```: float
