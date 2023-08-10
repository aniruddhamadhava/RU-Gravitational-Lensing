# RU Gravitational Lensing (2022 - 2025)
This repository contains the files/code necessary for the analysis done in the Fall 2022/Spring 2023 research. The following packages/codes are required (packages necessary only for data extraction are marked by "DE"):

* numpy
* scipy
* pygravlens
* shapely
* illustris_python (DE)
* astropy
* scipy

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

To compute deflection statistics for a cluster, we use the clensmod package. This file contains two classes: (i) lens, and (ii) clusterclass. The first class sets up the multiplane lensmodel containing all snapshots. The next class computes deflection statistics for a galaxy cluster. To set up the lensmodel, use the following code: 

            import clensmod as cl
            
            lensmodel = cl.lens(data_dictionary_file, kappa_file_directory, box_size)
            model, mybox = lensmodel.lensmodel(start_snap, end_snap, fov)

data_dictionary_file indicates the location where the data_dictionary (output from conv_comp) is stored. kappa_file_directory is the directory to where the kappa files (output from conv_comp) is stored. box_size is the simulation box size, expressed in Mpc/h. start_snap is the starting snapshot number, and end_snap is the ending snapshot number. fov is the field-of-view in degrees. If we want the code to also output the deflection array and kappa arrays, we pass an additional argument, "output = True". The default output is the multiplane lensmodel, and lensmodel box. Once we run these code, we can compute deflection statistics for a cluster as follows:

            cluster = cl.clusterclass(cluster_image_file, z_clus, clusname)
            info = cluster.DefStats(model,Nsamp= 10000, box = mybox)
            cluster.plot_DefStats()

The cluster_image_file is the .txt file containing image positions, z_clus is the cluster redshift, and clusname is the name of the cluster (e.g., A2744, M0416). info is a dictionary containing the covariance matrices:

            info = {'LOS': [defavg_los, defcov_los], 'FG': [defavg_fg, defcov_fg], 'LOSFG': [defcov_los_fg]}

cluster.plot_DefStats() plots the actual deflection statistics. On average, computing deflection statistics for each cluster takes ~ 30-40 minutes for 10,000 samples. 

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
