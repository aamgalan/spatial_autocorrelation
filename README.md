# Fast Spatial Autocorrelation

**Contact**

Anar Amgalan (anar.amgalan@stonybrook.edu) and Steven S. Skiena (skiena@cs.stonybrook.edu)

An example of a variable/feature values with spatial coordinates:

![County variables](plot-pub-75_spatial_county_coordinates_3D_Mollweide_3-variables.png)

![The trace of within-cluster squared deviations](plot-pub-75-15_spatial_skiena_trace_3-panel_lin-trans.png)

This repository contains the simple algorithm for fast computation of spatial autocorrelation. 
The statistic and the algorithm are introduced in ICDM 2020 (archive link or DOI identifier). 


## Dependencies


* python >= 3.6.8
* numpy >= 1.14.3
* scipy >= 1.1.0


## Datasets


The datasets used in the paper are 


* PPI
* PPI-large (a larger version of PPI)
* Reddit
* Flickr
* Yelp
* Amazon
* ogbn-products
* ... (more to be added)

