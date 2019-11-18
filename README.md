3D Couple Resemblance
======

This repository contains data and R analysis scripts for **__Distribution of facial resemblance in romantic couples suggests both positive and negative assortative processes influence human mate choice__**. Please note that due to data protection reasons, we cannot share original 3D face scans

Files
------

*Please note that some files in this repository are managed by [Git Large File Storage](https://git-lfs.github.com/ "Install Git LFS")*

* 1_preppingScanData: R code for submitting delineated and re-warped face scans to Generalized Procrustes Analysis

* 2_mainAnalysisScript: R code and output file for all analyses and plots reported in manuscript, as well as additional visualizations and analyses. The main output file can also be directly accessed in a browser at: http://facelab.org/3D-couple-resemblance/results.html - please allow time for loading!

Please see 
[session_info_20191118.txt](../master/session_info_20191118.txt) for details on used R packages and versions


Acknowledgments
------
Many of the custom functions used in our analyses were adapted from other packages

* **Select PCs according to broken stick model**: adapted from function `evplot`

   Francois Gillet, http://adn.biol.umontreal.ca/~numericalecology/numecolR/
   
* **Re-create 3D mesh from point cloud**: adapted from `nat::as.mesh3d.ashape3d`

  Jefferis, G. S. X. E. & Manton, J. D. (2014). NeuroAnatomy Toolbox v1.5.2. ZENODO. https://doi.org/10.5281/zenodo.10171
  
* **Plot 3D principal components**: adapted from `geomorph::plotTangentSpace`

   Adams, D. C., Collyer, M. L. & Kaliontzopoulou, A. (2019) Geomorph: Software for geometric morphometric analyses. R package version 3.1.0. https://cran.r-project.org/package=geomorph
