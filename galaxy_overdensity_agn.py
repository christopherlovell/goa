"""
Calculate galaxy overdensity
- can calculate for random regions, or for each galaxy. If random, set random flag to 'True'.
"""

from collections import Counter

import sys

import pandas as pd
import numpy as np
import itertools as it
from scipy.spatial.distance import cdist

from hightolowz import distance

# https://github.com/patvarilly/periodic_kdtree
from periodic_kdtree import PeriodicCKDTree

print sys.argv[1:]

selection_str = sys.argv[1]
redshift_str = sys.argv[2]
random = bool(int(sys.argv[3]))
print random

n = 100

# selection_str = 'sfr'
# redshift_str = '3p10'
# random = False

print "Reading galaxy data..."

#directory = '/lustre/scratch/astro/cl478/protoclusters_data/henriques2015a_z%s_%s.csv' % (redshift_str, selection_str)
directory = '~/sussex/protoclusters/data/r200/henriques2015a_z%s_%s_r200.csv' % (redshift_str, selection_str)
out_directory = '~/sussex/protoclusters/data/r200'

# selection_str = selection_str+'10'
# print selection_str

print directory
sys.stdout.flush()


# gals = pd.read_csv(directory, skiprows=104, skipfooter=1, engine='python')
gals = pd.read_csv(directory, skiprows=107, skipfooter=1, engine='python')


print "Filling in NaN values..."
gals.ix[np.isnan(gals['z0_haloId']), 'z0_haloId'] = -1
gals.ix[np.isnan(gals['z0_centralId']), 'z0_centralId'] = -1
gals.ix[np.isnan(gals['z0_central_mcrit200']), 'z0_central_mcrit200'] = 0

L = 480.279

if random:
    print "Initialising random regions..."
    N = 100000
    coods = pd.DataFrame(np.random.rand(N,3) * L, columns=['zn_x','zn_y','zn_z'])
    location_str = 'random'
else:
    print "Copying galaxy coordinates..."
    coods = gals[['zn_x','zn_y','zn_z']].copy()
    location_str = 'gals'

dimensions = np.array([L, L, L])

T = PeriodicCKDTree(dimensions, gals[['zn_x','zn_y','zn_z']])

r = [20, 12.5, 7.5, 5]
r_str = ['20', '12.5', '7.5', '5']

ngal = {'20': [None] * len(coods), '12.5': [None] * len(coods), '7.5': [None] * len(coods), '5': [None] * len(coods)}
nagn = {'20': [None] * len(coods), '12.5': [None] * len(coods), '7.5': [None] * len(coods), '5': [None] * len(coods)}
dgal = {'20': [None] * len(coods), '12.5': [None] * len(coods), '7.5': [None] * len(coods), '5': [None] * len(coods)}
max_fraction = {'20': [None] * len(coods), '12.5': [None] * len(coods), '7.5': [None] * len(coods), '5': [None] * len(coods)}
max_fraction_mass = {'20': [None] * len(coods), '12.5': [None] * len(coods), '7.5': [None] * len(coods), '5': [None] * len(coods)}
n_cluster_desc = {'20': [None] * len(coods), '12.5': [None] * len(coods), '7.5': [None] * len(coods), '5': [None] * len(coods)}
frac_cluster_desc = {'20': [None] * len(coods), '12.5': [None] * len(coods), '7.5': [None] * len(coods), '5': [None] * len(coods)}

# AGN constants
epsilon = 0.1
c = 2.97 * 10**8 # m s^-1

Mdot = (gals['zn_quasarAccretionRate'] + gals['zn_radioAccretionRate']) * 1.989 * 10**30 / (365.25*24*60*60)

Lbol = (epsilon * Mdot * c**2) / 1e-7  # erg s^-1
gals['agn_lum'] = Lbol

agn_lim = 1e46


print "Counting galaxies..."
sys.stdout.flush()

# can't calculate distances all in one go, so need to chunk
#for i,gals in z6_galaxies_mstar.groupby(np.arange(len(z6_galaxies_mstar))//n):
for j,c in coods.groupby(np.arange(len(coods))//n):
    
    if j % 80 == 0:
        print round(float(c.shape[0] * (j+1)) / coods.shape[0] * 100, 2), '%'

    for R, R_str in zip(r, r_str):

        gal_index = T.query_ball_point(c, r=R)

        start_index = (j*n)
        end_index = (j*n)+len(gal_index)

        ngal[R_str][start_index:end_index] = [len(x) for x in gal_index]

        for i in range(len(gal_index)):
            
            counter = Counter(gals.ix[gal_index[i]]['z0_central_mcrit200']).most_common()
          
            if len(counter) == 0:
                max_fraction_mass[R_str][start_index+i] = 0
                n_cluster_desc[R_str][start_index+i] = 0
                frac_cluster_desc[R_str][start_index+i] = 0
                max_fraction[R_str][start_index+i] = 0
                nagn[R_str][start_index+i] = 0
            else:
                total = sum([x[1] for x in counter])
    
                nagn[R_str][start_index+i] = np.sum(gals.ix[gal_index[i]]['agn_lum'] > agn_lim)
                max_fraction_mass[R_str][start_index+i] = counter[0][0]
                cluster_descendants = [(x[0] > 1e4) for x in counter]
                n_cluster_desc[R_str][start_index+i] = np.sum(cluster_descendants)
                frac_cluster_desc[R_str][start_index+i] = float(np.sum([x[1] for x in counter if (x[0] > 1e4)])) / total
                max_fraction[R_str][start_index+i] = float(counter[0][1]) / total


for R, R_str in zip(r, r_str):
    
    # print np.where([x is None for x in ngal[R_str]]) 

    print "Saving data..."
    print "R: ", R_str

    avg = float(gals.shape[0]) / L**3 * 4./3 * np.pi * R**3

    print "Average density: ", avg, "\n ........ "

    # delta_galaxy
    dgal[R_str][:] = (np.array(ngal[R_str]) - avg) / avg

    df = pd.DataFrame(np.array([dgal[R_str], ngal[R_str], max_fraction[R_str], max_fraction_mass[R_str], n_cluster_desc[R_str], nagn[R_str]]).T,
                     columns=('delta_gal_%s' % R_str,
                              'ngal_%s' % R_str,
                              'max_fraction_%s' % R_str,
                              'max_fraction_mass_%s' % R_str,
                              'n_cluster_desc_%s' % R_str,
                              'n_agn_%s' % R_str))

    df.to_csv('%s/dgal_%s_%s_r%s_%s.csv' % (out_directory, selection_str, redshift_str, R_str, location_str), index=False)

    print 'Saved to %s/dgal_%s_%s_r%s_%s.csv' % (out_directory, selection_str, redshift_str, R_str, location_str)

print "Complete!"
