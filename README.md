# Fall-2022-Research

This repository contains the files/code necessary for the analysis done in the Fall 2022 research. The following packages/codes are required:

* numpy
* scipy
* pygravlens

Abstract: 

Current lensing models assume a uniform density along a line-of-sight except at a mass distribution (i.e.,
lens). The implication is that light rays will only deflect near such distributions, traveling in a straight line
everywhere else. However, large-scale cosmological structures create density fluctuations, which can produce
additional lensing effects. Our project seeks to quantify these effects. In particular, we want to compute
the variance in deflection angles for nearby line-of-sights. This information will tell us the degree to which
current lens models are accurate and whether additional parameters are needed for the lens equation. For
this analysis, we use the TNG 100-3 simulation and pygravlens.


Keywords: Gravitational Lensing Models, Large-Scale Structures, Deflection Angles, TNG 100-3 Simulations, pygravlens

Files:

The file data_ext.py contains the code that is used to extract data from the Illustris TNG 100-3 simulation. This code is run in the Illustris TNG JupyterLab workspace. To use this file, use the following code:

            import numpy as np
            import data_ext as dext
            
Then, to extract data using this file, use the following code:

            rotmat = dext.rotator() # It is important that this code is only run once so that snapshots of the same number get the same rotation matrix applied.

            d = dext.dEx(snap, snapclo, snapfar, sbox_size, particle_type, rotmat) # Make sure to fill in the individual inputs as needed. 
            d.z_restrictor_bounds()
            d.data_extractor()
            d.data_reducer1()
            d.hist()

Note that the amount of time it takes the data to load may vary substantially from the examples based on the server speed. I have recorded a maximum time of close to 20 minutes and a minimum time of 30 seconds. Once all of the required files have been saved (mass densities for dark matter, gas, stars, and the centers for each snapshot), download and upload them (I suggest using Jupyter Notebook for this). Then, we can use the file, conv_comp.py. First, make sure the file: pygravlens.py (credit: Professor Keeton) is downloaded. Next, run the following code in a Jupyter notebook:

            import conv_comp as cc
            
            mtoconv = cc.m_to_conv(snap)
            mtoconv.dosnap2()

Note that the input, snap, is an array of snapshot numbers. E.g., if we just want to calculate the lensing maps for snapshot 52, the proper input is snap = [52]. If, however, we wanted to calculate the lensing maps for snapshots x, ..., y, we need to use snap = np.arange(x, y + 1, 1).

The following site links to an (ongoing) summary of the research methodologies used: 
https://www.overleaf.com/read/qcwdxgsgmncs

Credits:

* Professor Keeton: https://github.com/chuckkeeton/pygravlens
* Illustris TNG:
  * First results from the IllustrisTNG simulations: a tale of two elements - chemical evolution of magnesium and europium: Naiman, Jill (et al)
  * First results from the IllustrisTNG simulations: the stellar mass content of groups and clusters of galaxies:  Pillepich, Annalisa (et al)
  * First results from the IllustrisTNG simulations: radio haloes and magnetic fields: Marinacci, Federico (et al)
  * First results from the IllustrisTNG Simulations: the galaxy color bimodality: Nelson, Dylan (et al)
  * First results from the IllustrisTNG simulations: matter and galaxy clustering: Springel, Volker (et al)
