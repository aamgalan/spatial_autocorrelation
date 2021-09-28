# Fast Spatial Autocorrelation

**Contacts**

Anar Amgalan (anar.amgalan@alumni.stonybrook.edu) and Steven S. Skiena (skiena@cs.stonybrook.edu)

Several examples of variable/feature values with spatial coordinates:

![County variables](plot-pub-75_spatial_county_coordinates_3D_Mollweide_3-variables.png)

![The trace of within-cluster squared deviations](plot-pub-75-15_spatial_skiena_trace_3-panel_lin-trans.png)

This repository contains a simple algorithm for fast computation of spatial autocorrelation. 
The statistic and the algorithm are introduced in ICDM 2020 (archive link: https://arxiv.org/abs/2010.08676). 


## Dependencies


* python3 
* numpy 
* scipy 
* disjoint_set


## Datasets


The county-level dataset used in the paper is mainly from and sources cited therein and in our publication:

Can Twitter be used to predict county excessive alcohol consumption rates?
Curtis B, Giorgi S, Buffone AEK, Ungar LH, Ashford RD, et al. (2018) 
Can Twitter be used to predict county excessive alcohol consumption rates?. 
PLOS ONE 13(4): e0194290. https://doi.org/10.1371/journal.pone.0194290

The dataset used in the example of time-series features with preserved coordinates is from:

Lilianne R. Mujica-Parodi, Anar Amgalan, Syed Fahad Sultan, Botond Antal, Xiaofei Sun, Steven Skiena, Andrew Lithen, Noor Adra, Eva-Maria Ratai, Corey Weistuch, Sindhuja Tirumalai Govindarajan, Helmut H. Strey, Ken A. Dill, Steven M. Stufflebeam, Richard L. Veech, Kieran Clarke
Proceedings of the National Academy of Sciences Mar 2020, 117 (11) 6170-6177; DOI: https://doi.org/10.1073/pnas.1913042117 
