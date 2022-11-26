#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import illustris_python as il
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import binned_statistic_2d
from astropy.cosmology import Planck15 as cosmo
import astropy.units as u
from astropy.constants import c, G
from numpy.random import randint

def rotator():
    """
    Line 5 loads the .npy file that contains the rotation matrices. Then, each snapshot (50 to 99) is randomly assigned
    a rotation matrix. This function should only be run once so that all three data slices for each snapshot have the 
    same rotation matrix applied to them.
    """
    allrot = np.load('rotation.npy', allow_pickle = True)
    rotmat = {}
    for i in range(50, 100):
        rotmat[i] = allrot[randint(0, 23)]
    return rotmat

class dEx:
    def __init__(self, snap, snapclo, snapfar, Sbox_size, particle_type, rotmat):
        """
        snap: Snapshot number
        snapclo: Snapshot number with z = 0.0
        snapfar: Snapshot number with z = 1.0
        particle_type: Three options: (1) Dark Matter, (2) Gas, and (3) Stars
        rotmat: Rotation matrix output from rotator() function.
        This function needs to be run each time the particle type and/or snapshot number is changed.
        """
        self.snap = snap
        self.snapclo = snapclo
        self.snapfar = snapfar
        self.Sbox_size = Sbox_size
        self.particle_type = particle_type
        self.rotmat = rotmat
    def z_restrictor_bounds(self):
        """
        (1) Creates an array of redshifts (i.e., reds) (going from snapshot 49 - 99). This array is specific to the Illustris TNG 100-3 simulation.
        (2) Creates a dictionary (i.e., zsnap), linking each snapshot number with its respective redshift.
        (3) Creates a dictionary (i.e., Lsnap), linking each spapshot number to its associated distance in comoving kpc
        (4) Finds the upper and lower bound for the z-values (will be used to restrict the data)
        """
        reds = [1.04, 1.0,0.95,0.92,0.89,0.85,0.82,0.79,0.76,0.73,0.70,0.68,0.64,0.62,0.60,0.58,0.55,0.52,0.50,0.48,0.46,0.44,0.42,0.40,0.38,0.36,0.35,0.33,0.31,0.30,0.27,0.26,0.24,0.23,0.21,0.20,0.18,0.17,0.15,0.14,0.13,0.11,0.10,0.08,0.07,0.06,0.05,0.03,0.02,0.01,0.00]
        zsnap = {}
        for i in np.arange(49, 100,1):
            zsnap[i] = reds[i - 49]
        Lsnap = {}
        for i in zsnap.keys():
            Lsnap[i] = cosmo.comoving_distance(zsnap[i]).to('kpc').value*cosmo.h
        if self.snap < self.snapclo and self.snap > self.snapfar:
            zplo = 0.5 * (Lsnap[self.snap + 1] - Lsnap[self.snap])
            zphi = 0.5 * (Lsnap[self.snap -1] - Lsnap[self.snap])
        elif self.snap == self.snapclo:
            zplo = 0.5* Lsnap[self.snap]
            zphi = 0.5 * (Lsnap[self.snap - 1] -Lsnap[self.snap])
        elif self.snap == self.snapfar:
            zplo = 0.5 * (Lsnap[self.snap + 1] - Lsnap[self.snap])
            zphi = 0.5 * Lsnap[self.snap]
        zplo += self.Sbox_size/2
        zphi += self.Sbox_size/2
        self.zplo = zplo
        self.zphi = zphi
    def data_extractor(self):
        """
        Extracts the positions and masses of particles for each particle type and snapshot number. This section is
        specific to the Illustris TNG 100-3 simulation.
        """
        header = il.groupcat.loadHeader('../sims.TNG/TNG100-3/output/', self.snap)
        r = self.rotmat[self.snap]
        if self.particle_type == 'dm':
            data = il.snapshot.loadSubset('../sims.TNG/TNG100-3/output/', self.snap, 'dm', fields = ['Coordinates'])
            data = data@r.T
            mass = np.full(len(data),0.0323567549719664) 
            self.coords = data
            self.mass = mass
            self.header = header
        elif self.particle_type == 'stars':
            data = il.snapshot.loadSubset('../sims.TNG/TNG100-3/output/', self.snap, 'stars', fields = ['Coordinates', 'Masses'])
            data['Coordinates'] = data['Coordinates']@r.T
            self.coords = data['Coordinates']
            self.mass = data['Masses']
            self.header = header
        else:
            data = il.snapshot.loadSubset('../sims.TNG/TNG100-3/output/', self.snap, 'gas', fields = ['Coordinates', 'Masses'])
            data['Coordinates'] = data['Coordinates']@r.T
            self.coords = data['Coordinates']
            self.mass = data['Masses']
            self.header = header
    def data_reducer1(self):
        """
        Uses the upper and lower bounds for the z-values obtained from the z_restrictor_bounds() function and chooses
        points in each simulation whose z-values lie within this range. This is to avoid double-counting the same region
        when calculating the mass densities.
        """
        if self.particle_type == 'dm':
            if np.mean(self.coords[:, 2]) < 0:
                self.flag1 = (self.coords[:,2]+self.Sbox_size>=self.zplo)
                self.flag2 = (self.coords[:,2]+self.Sbox_size<=self.zphi)
                self.restricted_data = self.coords[self.flag1*self.flag2]
                self.mass_rest = np.full(len(self.restricted_data),0.0323567549719664) 
            else:
                self.flag1 = (self.coords[:,2]>=self.zplo)
                self.flag2 = (self.coords[:,2]<=self.zphi)
                self.restricted_data = self.coords[self.flag1*self.flag2]
                self.mass_rest = np.full(len(self.restricted_data),0.0323567549719664)
            self.x_rest =self.restricted_data[:, 0]
            self.y_rest =self.restricted_data[:, 1] 
        else:
            if np.mean(self.coords[:, 2]) < 0:
                self.flag1 = (self.coords[:, 2] + self.Sbox_size >= self.zplo)
                self.flag2 = (self.coords[:, 2] + self.Sbox_size <= self.zphi)
                self.restricted_data = self.coords[self.flag1 * self.flag2]
                self.mass_rest = self.mass[:len(self.restricted_data):]
                self.x_rest =self.restricted_data[:, 0]
                self.y_rest =self.restricted_data[:, 1]
            else:
                self.flag1 = (self.coords[:, 2]>= self.zplo)
                self.flag2 = (self.coords[:, 2]<= self.zphi)
                self.restricted_data = self.coords[self.flag1 * self.flag2]
                self.mass_rest = self.mass[:len(self.restricted_data):]
                self.x_rest =self.restricted_data[:, 0]
                self.y_rest =self.restricted_data[:, 1] 
    def hist(self):
        """
        Calculates the mass density for each snapshot, saves the density and the x-, y-centers to be used in the 
        Fourier Analysis. 
        """
        centers = {}
        grid, xedge, yedge, _ = binned_statistic_2d(self.x_rest, self.y_rest,self.mass_rest, 'sum', bins=[600, 600], range=[[np.min(self.x_rest),np.max(self.x_rest)],[np.min(self.y_rest),np.max(self.y_rest)]])
        xcen = (xedge[:-1]+xedge[1:])/2
        ycen = (yedge[:-1]+yedge[1:])/2

        pxSize = self.header['BoxSize'] / 600 # code units
        pxSize_kpc = pxSize * self.header['Time'] / self.header['HubbleParam']
        pxArea = pxSize_kpc**2

        grid_log_msun_per_kpc2 = grid * 1e10 / self.header['HubbleParam'] / pxArea
        if self.particle_type == 'dm':
            np.save(f'density_{self.snap}_dm.npy', grid_log_msun_per_kpc2)
        elif self.particle_type == 'stars':
            np.save(f'density_{self.snap}_stars.npy', grid_log_msun_per_kpc2)
        else:
            np.save(f'density_{self.snap}_gas.npy', grid_log_msun_per_kpc2)
        del grid_log_msun_per_kpc2
        if self.particle_type == 'dm':
            centers['x'] = xcen
            centers['y'] = ycen
            np.save(f'center_{self.snap}.npy', centers)
        del xcen, ycen

