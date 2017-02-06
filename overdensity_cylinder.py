"""
Calculate galaxy overdensity
- can calculate for random regions, or for each galaxy. If random, set random flag to 'True'.
- finds overdensity in redshift space for a cylinder with radius R, depth deltac
"""

import sys

import pandas as pd
import numpy as np

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

# selection_str = 'sfr'
# redshift_str = '3p10'
# random = False

print selection_str
print random

# selection_str = selection_str+'10'
n = 100       # chunk length
N = 100    # number of random regions
L = 480.279   # box side length

z = float(redshift_str.replace('p','.'))
dimensions = np.array([L, L, L])

#directory = '/lustre/scratch/astro/cl478/protoclusters_data/henriques2015a_z%s_%s.csv' % (redshift_str, selection_str)
directory = '~/sussex/protoclusters/data/r200/henriques2015a_z%s_%s_r200.csv' % (redshift_str, selection_str)
out_directory = '~/sussex/protoclusters/data/r200'

print "dir:", directory
sys.stdout.flush()

print "Reading galaxy data..."
gals = pd.read_csv(directory, skiprows=120, skipfooter=1, engine='python')

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

# print gal_coods
print Planck13.H(z)
gal_coods['zn_z'] += gals['zn_velZ'] * (1+z) / Planck13.H(z).value

print "Building KDtree..."
T = PeriodicCKDTree(dimensions, gal_coods[['zn_x','zn_y','zn_z']])

avg = float(gals.shape[0]) / L**3 # average overdensity cMpc^-3

clim = 0.4
plim = 0.4

r = [2.5, 5, 7.5, 10, 15]

label = np.zeros((len(r), len(coods)))
dgal = np.zeros((len(r), len(coods)))
completeness = np.zeros((len(r), len(coods)))
purity = np.zeros((len(r), len(coods)))

print "Counting galaxies... (", len(coods),")"
sys.stdout.flush()

 # loop through radii (r)
for Ridx, R in enumerate(r):

    print "R:",R,"\n ---------------------------------"    

    # set deltaz equal to radius (can optionally change deltaz)
    half_deltac = R

    vol_avg = np.pi * R**2 * half_deltac * avg  # average overdensity in chosen volume

    # can't calculate distances all in one go, so need to chunk
    for j,c in coods.groupby(np.arange(len(coods))//n):

        # print progress
        if j % 10 == 0:
            print round(float(c.shape[0] * (j+1)) / coods.shape[0] * 100, 2), '%'


        # find all galaxies within a sphere of radius the max extent of the cylinder
        # deltac = 2 * max(r) # TEMP while deltac not changing
        gal_index = T.query_ball_point(c, r=(R**2 + half_deltac**2)**0.5)

        # filter by cylinder using norm_coods()
        gal_index = [np.array(gal_index[k])[norm_coods(gal_coods.ix[gal_index[k]].values, c.ix[k + j*n].values, R, half_deltac, L)] for k in range(n)]

        # # save start and end indices
        start_index = (j*n)

        dgal[Ridx, start_index:(start_index+n)] = (np.array([len(x) for x in gal_index]) - vol_avg) / vol_avg

        for i in range(len(gal_index)):

            cluster_ids = Counter(gals.ix[gal_index[i]][gals.ix[gal_index[i]]['z0_central_mcrit200'] > 1e4]['z0_centralId'])

            if len(cluster_ids) > 0:

                cstats = np.zeros((len(cluster_ids), 2))

                for k, (q, no) in enumerate(cluster_ids.items()):
                    cluster_gals = gals.ix[gals['z0_centralId'] == q]
                    cstats[k,0] = float(no) / len(cluster_gals)  # completeness
                    cstats[k,1] = float(no) / len(gal_index)  # purity

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
                completeness[Ridx, start_index+i] = cstats[max_completeness[0], 0]
                purity[Ridx, start_index+i] = cstats[max_purity[0], 1]

                # label candidates
                if cstats[max_completeness[0], 0] >= clim: # if completeness high
                    if cstats[max_purity[0], 1] >= plim: # ..and purity high
                        label[Ridx, start_index+i] =  1 # 'protocluster'
                    else: # ..and purity low
                        label[Ridx, start_index+i] = 2 # 'pc in field'
                else: # if completeness low
                    if cstats[max_purity[0], 1] >= plim:  # ..and purity high
                        label[Ridx, start_index+i] = 3 # 'part protocluster'
                    else: # ..and purity low
                        label[Ridx, start_index+i] = 0 # 'field'

            else:
                completeness[Ridx, start_index+i] = 0.
                purity[Ridx, start_index+i] = 0.
                label[Ridx, start_index+i] = 0 # 'field'



for Ridx, R in enumerate(r):

    print "Saving data (R:",R,")..."

    df = pd.DataFrame(np.array([label[Ridx], dgal[Ridx], completeness[Ridx], purity[Ridx]]).T,
                  columns = ['label_%s'%R, 'dgal_%s'%R,
                             'completeness_%s'%R, 'purity_%s'%R])

    df.to_csv('%s/dgal_%s_%s_r%s_%s.csv' % (out_directory, selection_str, redshift_str, R, location_str), index=False)

    print 'Saved to %s/dgal_%s_%s_r%s_%s.csv' % (out_directory, selection_str, redshift_str, R, location_str)

print "Complete!"
