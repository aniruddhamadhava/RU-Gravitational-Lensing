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

Once all of the required files have been saved, download and upload them (I suggest use Jupyter Notebook for this). Then, we can use the MtoK.ipynb file. First, make sure the file: pygravlens.py (credit: Prof. Charles Keeton) is downloaded. Next, make sure to run all of the functions in the MtoK.ipynb file. Then, to convert the mass density map of a given snapshot to the lensing maps, use the following code:

            dosnap(snap) # Change snap to get the required snapshot number

Edit: A new .ipynb notebook has been uploaded to demonstrate how the new data_ext.py code works.
