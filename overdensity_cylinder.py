"""
Calculate galaxy overdensity
- can calculate for random regions, or for each galaxy. If random, set random flag to 'True'.
- finds overdensity in redshift space for a cylinder with radius R, depth deltac
"""

import sys

import pickle as pcl

import numpy as np
import pandas as pd

from collections import Counter

# https://github.com/patvarilly/periodic_kdtree
from periodic_kdtree import PeriodicCKDTree

# TODO: check this is the exact cosmology...
from astropy.cosmology import Planck13

from norm_cood import norm_coods


print sys.argv[1:]

selection_str = sys.argv[1]
redshift_str = sys.argv[2]
random = bool(int(sys.argv[3]))
r = [float(sys.argv[4])]
# half_deltac = [float(sys.argv[5])]

# selection_str = 'sfr'
# redshift_str = '3p10'
# random = False
r_str = str(r[0]).replace('.','p')

print "selection:",selection_str
print "random?:", random
print "r:",r_str

# selection_str = selection_str+'10'
n = 100       # chunk length
N = 200000    # number of random regions
L = 480.279   # box side length

z = float(redshift_str.replace('p','.'))
dimensions = np.array([L, L, L])

directory = '/lustre/scratch/astro/cl478/protoclusters_data/henriques2015a_z%s_%s_r200.csv' % (redshift_str, selection_str)
# directory = '~/sussex/protoclusters/data/r200/henriques2015a_z%s_%s_r200.csv' % (redshift_str, selection_str)
#directory = '~/protoclusters/data/r200/henriques2015a_z%s_%s_r200.csv' % (redshift_str, selection_str)

out_directory = '/lustre/scratch/astro/cl478/protoclusters_data'
# out_directory = '/home/chris/sussex/protoclusters/data/r200'
#out_directory = '~/protoclusters/data/r200'

print "dir:", directory
sys.stdout.flush()

print "Reading galaxy data..."
gals = pd.read_csv(directory, skiprows=122, skipfooter=1, engine='python')

print "filtering by stellar mass..."
print gals.shape
selection_str += str(10)
print selection_str
gals = gals[gals['zn_stellarMass'] > 1].reset_index(drop=True)
print gals.shape

print "Filling in NaN values..."
gals.ix[np.isnan(gals['z0_haloId']), 'z0_haloId'] = -1
gals.ix[np.isnan(gals['z0_centralId']), 'z0_centralId'] = -1
gals.ix[np.isnan(gals['z0_central_mcrit200']), 'z0_central_mcrit200'] = 0

if random:
    print "Initialising random regions..."
    coods = pd.DataFrame(np.random.rand(N,3) * L, columns=['zn_x','zn_y','zn_z'])
    location_str = 'random'
else:
    print "Copying galaxy coordinates..."
    coods = gals[['zn_x','zn_y','zn_z']].copy()
    location_str = 'gals'

# Convert z-axis to redshift space
gal_coods = gals[['zn_x','zn_y','zn_z']].copy()

print Planck13.H(z)
gal_coods['zn_z'] += gals['zn_velZ'] * (1+z) / Planck13.H(z).value

# if np.sum(gal_coods['zn_z'] < 0) > 0:
#     print 'negative zees:',np.sum(gal_coods['zn_z'] < 0)
# if np.sum(gal_coods['zn_z'] > L) > 0:
#     print 'big zees:',np.sum(gal_coods['zn_z'] > L)

# fix z coordinate positions if they fall outside the box
gal_coods.loc[gal_coods['zn_z'] < 0,'zn_z'] = gal_coods.loc[gal_coods['zn_z'] < 0,'zn_z'] + L
gal_coods.loc[gal_coods['zn_z'] > L,'zn_z'] = gal_coods.loc[gal_coods['zn_z'] > L,'zn_z'] - L

print "Building KDtree..."
T = PeriodicCKDTree(dimensions, gal_coods[['zn_x','zn_y','zn_z']])

avg = float(gals.shape[0]) / L**3 # average overdensity cMpc^-3

# r = [4, 7, 11]
half_deltac = [5.1, 3.5, 6.2, 3.3, 5.2, 3.8, 4.4, 1.9, 4.7]

# numpy multidiemnsional array 
# (radii, depth, coordinates, [overdensity, completeness, purity, descendant mass])
out_stats = np.zeros((len(r), len(half_deltac), len(coods), 4), dtype=np.float16)

print "Counting galaxies... (", len(coods),")"
sys.stdout.flush()

 # loop through radii (r)
for Ridx, R in enumerate(r):

    for idxc, dc in enumerate(half_deltac):

        print "R:",R,"| half_deltac:",dc
        sys.stdout.flush()

        # set deltaz equal to radius (can optionally change deltaz)
        # half_deltac = R

        vol_avg = np.pi * R**2 * (2*dc) * avg  # average overdensity in chosen volume

        # can't calculate distances all in one go, so need to chunk
        for j,c in coods.groupby(np.arange(len(coods))//n):

            # print progress
            if j % 100 == 0:
                print round(float(c.shape[0] * (j+1)) / coods.shape[0] * 100, 2), '%'
                sys.stdout.flush()


            # find all galaxies within a sphere of radius the max extent of the cylinder
            gal_index = T.query_ball_point(c, r=(R**2 + dc**2)**0.5)

            # filter by cylinder using norm_coods()
            gal_index = [np.array(gal_index[k])[norm_coods(gal_coods.ix[gal_index[k]].values, c.ix[k + j*n].values, R, dc, L)] for k in range(len(c))]

            # # save start and end indices
            start_index = (j*n)

            out_stats[Ridx, idxc, start_index:(start_index+len(c)), 0] = (np.array([len(x) for x in gal_index]) - vol_avg) / vol_avg

            for i in range(len(gal_index)):

                cluster_ids = Counter(gals.ix[gal_index[i]][gals.ix[gal_index[i]]['z0_central_mcrit200'] > 1e4]['z0_centralId'])

                if len(cluster_ids) > 0:

                    

                    cstats = np.zeros((len(cluster_ids), 2))

                    for k, (q, no) in enumerate(cluster_ids.items()):
                        cluster_gals = gals.ix[gals['z0_centralId'] == q]
                        cstats[k,0] = float(no) / len(cluster_gals)  # completeness
                        cstats[k,1] = float(no) / len(gal_index[i])  # purity


                    # find id of max completeness and purity in cstats array
                    max_completeness = np.where(cstats[:,0] == cstats[:,0].max())[0]
                    max_purity = np.where(cstats[:,1] == cstats[:,1].max())[0]

                    # sometimes multiple clusters can have same completeness or purity in a single candidate
                    # - use the cluster with the highest complementary completeness/purity
                    if len(max_completeness) > 1:

                        # get matches between completeness and purity
                        matches = [x in max_purity for x in max_completeness]

                        if np.sum(matches) > 0:
                            # just use the first one
                            max_completeness = [np.where(matches)[0][0]]
                            max_purity = [np.where(matches)[0][0]]
                        else:
                            max_completeness = [max_completeness[np.argmax(cstats[max_completeness, 1])]]

                    if len(max_purity) > 1:

                        matches = [x in max_completeness for x in max_purity]

                        if np.sum(matches) > 0:
                            max_completeness = [np.where(matches)[0][0]]
                            max_purity = [np.where(matches)[0][0]]

                        else:
                            max_purity = [max_purity[np.argmax(cstats[max_completeness, 0])]]


                    # sometimes the cluster with the highest completeness does not have the highest purity, or vice versa
                    # - use the cluster with the highest combined purity/completeness added in quadrature
                    if max_completeness[0] != max_purity[0]:
                        max_completeness = [np.argmax([pow(np.sum(x**2), 0.5) for x in cstats])]
                        max_purity = max_completeness

                    # save completeness and purity values
                    out_stats[Ridx, idxc, start_index+i, 1] = cstats[max_completeness[0], 0]
                    out_stats[Ridx, idxc, start_index+i, 2] = cstats[max_purity[0], 1]

                    # save descendant mass
                    # filter by cluster id, save z0 halo mass
                    # can use either max_completeness or max_purity, both equal by this point
                    
                    out_stats[Ridx, idxc, start_index+i, 3] = gals.loc[gals['z0_centralId'] == cluster_ids.keys()[max_completeness[0]], 'z0_central_mcrit200'].iloc[0]
                    

                else:
                    out_stats[Ridx, idxc, start_index+i, 1] = 0.
                    out_stats[Ridx, idxc, start_index+i, 2] = 0.
                    out_stats[Ridx, idxc, start_index+i, 3] = np.nan


print "Saving data..."

pcl.dump(out_stats, open('%s/dgal_%s_z%s_r%s_%s.pcl' % (out_directory, selection_str, redshift_str, r_str, location_str), 'wb'))

print 'Saved to %s/dgal_%s_z%s_r%s_%s.pcl' % (out_directory, selection_str, redshift_str, r_str, location_str)
print "Complete!"

