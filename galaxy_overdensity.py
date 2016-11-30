"""
Calculate galaxy overdensity
- can calculate for random regions, or for each galaxy. If random, set random flag to 'True'.
"""

random = True

import sys

import pandas as pd
import numpy as np
import itertools as it
from scipy.spatial.distance import cdist

from hightolowz import distance

# https://github.com/patvarilly/periodic_kdtree
from periodic_kdtree import PeriodicCKDTree

n = 100
selection_str = 'sfr'
redshift_str = '2p07'

print "Reading galaxy data..."
print directory


directory = '/lustre/scratch/astro/cl478/protoclusters_data/henriques2015a_z2p07_mstar.csv'

gals = pd.read_csv(directory, skiprows=104, skipfooter=1, engine='python')

print "Filling in NaN values..."
gals.ix[np.isnan(gals['z0_haloId']), 'z0_haloId'] = -1
gals.ix[np.isnan(gals['z0_centralId']), 'z0_centralId'] = -1
gals.ix[np.isnan(gals['z0_central_mcrit200']), 'z0_central_mcrit200'] = 0

L = 480.279

if random:
    print "Initialising random regions..."
    N = 30000
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
dgal = {'20': [None] * len(coods), '12.5': [None] * len(coods), '7.5': [None] * len(coods), '5': [None] * len(coods)}
max_fraction = {'20': [None] * len(coods), '12.5': [None] * len(coods), '7.5': [None] * len(coods), '5': [None] * len(coods)}
max_fraction_mass = {'20': [None] * len(coods), '12.5': [None] * len(coods), '7.5': [None] * len(coods), '5': [None] * len(coods)}
n_cluster_desc = {'20': [None] * len(coods), '12.5': [None] * len(coods), '7.5': [None] * len(coods), '5': [None] * len(coods)}


print "Counting galaxies..."

# can't calculate distances all in one go, so need to chunk
#for i,gals in z6_galaxies_mstar.groupby(np.arange(len(z6_galaxies_mstar))//n):
<<<<<<< HEAD
for j,c in coods.groupby(np.arange(len(coods))//n):
    
    if j % 5 == 0:
        print round(float(c.shape[0] * (j+1)) / coods.shape[0] * 100, 2), '%'
        sys.stdout.flush()

    # calculate distances
    #dist = np.vstack(c.apply(lambda x: distance(x, gals[['zn_x','zn_y','zn_z']], dimensions), axis=1))

    for R, R_str in zip(r, r_str):

        #gal_index = dist < R
        
        print "Query KDtree..."        
        gal_index = T.query_ball_point(c, r=R)
        
        start_index = (j*len(gal_index))
        end_index = (j*len(gal_index))+len(gal_index)
        
        print "Count gals..."
        ngal[R_str][start_index:end_index] = [len(x) for x in gal_index]
        
        print "Aggregating gals..."
        for i in range(len(gal_index)):

            #n_galaxy = len(gal_index[i])
            #ngal[R_str][j+i] = n_galaxy

            if ngal[R_str][i] == 0:
                m_max=0
                ncd=0
                agg_count=np.array([0])
                max_frac=0
            else:
                agg_mvir = gals.ix[gal_index[i]].groupby('z0_centralId', sort=False).mean()['z0_central_mcrit200']
                agg_count = gals.ix[gal_index[i]].groupby('z0_centralId', sort=False)['z0_centralId'].count().astype(float)

                agg = pd.DataFrame([agg_mvir, agg_count]).T

                m_max = agg.loc[agg['z0_centralId'].idxmax()]['z0_central_mcrit200'] # find mass of most common descendant
                ncd = sum(agg_mvir > 1e4)
                max_frac = agg_count.max() / agg_count.sum()

            max_fraction_mass[R_str][j+i] = m_max
            max_fraction[R_str][j+i] = max_frac
            n_cluster_desc[R_str][j+i] = ncd


for R, R_str in zip(r, r_str):

    print "Saving data..."
    print "R: ", R_str

    avg = float(gals.shape[0]) / L**3 * 4./3 * np.pi * R**3

    print "Average density: ", avg, "\n ........ "

    # delta_galaxy
    dgal[R_str] = (np.array(ngal[R_str]) - avg) / avg

    df = pd.DataFrame(np.array([dgal[R_str], ngal[R_str], max_fraction[R_str], max_fraction_mass[R_str], n_cluster_desc[R_str]]).T,
                     columns=('delta_gal_%s' % R_str,
                              'ngal_%s' % R_str,
                              'max_fraction_%s' % R_str,
                              'max_fraction_mass_%s' % R_str,
                              'n_cluster_desc_%s' % R_str))

    df.to_csv('data/planck1/dgal_%s_%s_r%s_%s.csv' % (selection_str, redshift_str, R_str, location_str), index=False)

    print 'Saved to data/planck1/dgal_%s_%s_r%s_%s.csv' % (selection_str, redshift_str, R_str, location_str)

print "Complete!"
