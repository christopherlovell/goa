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

from periodic_kdtree import PeriodicCKDTree  # https://github.com/patvarilly/periodic_kdtree
from astropy.cosmology import Planck13
from norm_cood import norm_coods

from methods import normalise_coods
from methods import z_distort
from methods import factor_h


def overdensity_cylinder(gals, coods, R, dc, L, pc_stats=False, cluster_mass_lim=1.48e4, n=100, verbose=False):
    """
    Find overdensity statistics over the whole simulation box for cylindrical apertures.

    Args:
        gals - dataframe of galaxy properties
        random - bool, whether to use random positions or center on all galaxies
        R - aperture radius, cMpc
        dc - half aperture depth, cMpc
        cluster_mass_lim - limiting descendant mass above which to classify clusters, z0_central_mcrit200
        n - chunk length
        N - number of random regions
    Returns:
        out_stats - output statistics, numpy array of shape [len(coods), 4], where coods is either the number of galaxies or number of random regions.
    """

    dimensions = np.array([L, L, L])

    if verbose: print "Building KDtree..."
    T = PeriodicCKDTree(dimensions, gals[['zn_x','zn_y','zn_z']])
    
    avg = float(gals.shape[0]) / L**3 # average overdensity cMpc^-3
    
    out_stats = np.zeros((len(coods),4))
    
    vol_avg = np.pi * R**2 * (2*dc) * avg  # average overdensity in chosen volume
    
    for j,c in coods.groupby(np.arange(len(coods))//n): # can't calculate distances all in one go, so need to chunk
    
        if verbose: # print progress
            if j % 100 == 0:
                print round(float(c.shape[0] * (j+1)) / coods.shape[0] * 100, 2), '%'
                sys.stdout.flush()

    
        # find all galaxies within a sphere of radius the max extent of the cylinder
        gal_index = T.query_ball_point(c, r=(R**2 + dc**2)**0.5)
    
        # filter by cylinder using norm_coods()
        gal_index = [np.array(gal_index[k])[norm_coods(gals.iloc[gal_index[k]][['zn_x','zn_y','zn_z']].values, c.ix[k + j*n].values, R=R, half_deltac=dc, L=L)] for k in range(len(c))]
    
        start_index = (j*n)  # save start index

        # calculate dgal
        out_stats[start_index:(start_index+len(c)), 0] = (np.array([len(x) for x in gal_index]) - vol_avg) / vol_avg
    
        if pc_stats:  # calculate completeness and purity statistics
    
            for i in range(len(gal_index)):

                cluster_ids = gals.iloc[gal_index[i]]
                cluster_ids = Counter(cluster_ids[cluster_ids['z0_central_mcrit200'] > cluster_mass_lim]['z0_centralId'])

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
                    out_stats[start_index+i, 1] = cstats[max_completeness[0], 0]  # completeness
                    out_stats[start_index+i, 2] = cstats[max_purity[0], 1]        # purity

                    # save descendant mass
                    # filter by cluster id, save z0 halo mass
                    # can use either max_completeness or max_purity, both equal by this point

                    out_stats[start_index+i, 3] = gals.loc[gals['z0_centralId'] == cluster_ids.keys()[max_completeness[0]], 'z0_central_mcrit200'].iloc[0]


                else:  # if no galaxies in aperture
                    out_stats[start_index+i, 1] = 0.
                    out_stats[start_index+i, 2] = 0.
                    out_stats[start_index+i, 3] = np.nan
    

    return out_stats

