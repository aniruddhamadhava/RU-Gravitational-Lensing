# Fall-2022-Research

This repository contains the files/code necessary for the analysis done in the Fall 2022 research.
The file data_ext.py contains the code that is used to extract data from the Illustris TNG 100-3 simulation. To use this file, use the following code:

            import numpy as np
            import data_ext as dext
            
Then, to extract data using this file, use the following code:

            rotmat = dext.rotator() # It is important that this code is only run once so that snapshots of the same number get the same rotation matrix applied.

            d = dext.dEx(snap, snapclo, snapfar, sbox_size, particle_type, rotmat) # Make sure to fill in the individual inputs as needed. 
            d.z_restrictor_bounds()
            d.data_extractor()
            d.data_reducer1()
            d.hist()

Once all of the required files have been saved (mass densities for dark matter, gas, stars, and the centers for each snapshot), download and upload them (I suggest use Jupyter Notebook for this). Then, we can use the file, conv_comp.py. First, make sure the file: pygravlens.py (credit: Prof. Charles Keeton) is downloaded. Next, run the following code in a Jupyter notebook:

            import conv_comp as cc
            
            mtoconv = cc.m_to_conv(snap)
            mtoconv.dosnap2()

Note that the input, snap, is an array of snapshot numbers. E.g., if we just want to calculate the lensing maps for snapshot 52, the proper input is snap = [52]. If, however, we wanted to calculate the lensing maps for snapshots x, ..., y, we need to use snap = np.arange(x, y + 1, 1).
