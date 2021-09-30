
## Explanation of config parameters for strain rate calculators

#### General
* **range_strain**: float/float/float/float, representing the target region for strain rate to be calculated upon, in W/E/S/N degrees longitude and latitude
* **range_data**: float/float/float/float, representing the region of velocity data that will be used for the calculation, in W/E/S/N degrees longitude and latitude
* **inc**: float/float, representing x/y grid spacing of resulting strain grid, in degrees longitude and latitude 

#### delaunay/ delaunay_flat: 
* no parameters

#### visr:

See: http://scec.ess.ucla.edu/~zshen/visr/visr.html

* **distance_weighting**: string
* **spatial_weighting**: string
* **min_max_inc_smooth**: float/float/float
* **executable**: string, relative path to location of compiled fortran executable, visr.exe or similar 

#### gpsgridder

See: http://gmt.soest.hawaii.edu/doc/latest/supplements/potential/gpsgridder.html

* **poisson**: float
* **fd**: float
* **eigenvalue**: float 


#### huang
* **EstimateRadiusKm**: float, radius of searching neighborhood in km
* **nstations**: integer, minimum number of stations within search radius 

#### Tape
* **code_dir**: string 
* **qmin**: integer
* **qmax**: integer
* **qsec**: integer

#### geostats
* **model_type**: string, one of [Gaussian, Exponential, Nugget].
* **sill_east**: float
* **range_east**: float
* **nugget_east**: float
* **sill_north**: float
* **range_north**: float
* **nugget_north**: float
* **trend**: float
