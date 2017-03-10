"""
Calculate protocluster completeness and purity at a given radii for different selections

"""

import matplotlib.pyplot as plt

from palettable.tableau import GreenOrange_12

cmap = GreenOrange_12.hex_colors

# TODO: check this is the exact cosmology...
from astropy.cosmology import Planck13

import pandas as pd
import numpy as np

import sys

# https://github.com/patvarilly/periodic_kdtree
from periodic_kdtree import PeriodicCKDTree

def cluster_stats(gals, L=500, redshift_distort=False, z=None):

    # clusters = gals[gals['z0_central_mcrit200'] > 1e4].groupby('z0_centralId')['z0_central_mcrit200','z0_centralId'].max()
    clusters = gals[gals['z0_central_mcrit200'] > 1e4].groupby('z0_centralId')['z0_central_mcrit200','z0_centralId'].max().reset_index(level=2, drop=True)

    # if len(cluster_ids) == 0:
    #     print "Selecting clusters..."
    #     cluster_ids = pd.unique(gals[gals['z0_central_mcrit200'] > 1e4]['z0_centralId'])

    print "N Clusters: ", len(clusters)

    cluster_stats = [None] * len(clusters)

    dimensions = np.array([L, L, L])

    gal_coods = gals[['zn_x','zn_y','zn_z']].copy()
    
    # Convert z-axis to redshift space
    if redshift_distort:
        gal_coods['zn_z'] += gals['zn_velZ'] * (1+z) / Planck13.H(z).value

        # fix z coordinate positions if they fall outside the box
        gal_coods.loc[gal_coods['zn_z'] < 0,'zn_z'] = gal_coods.loc[gal_coods['zn_z'] < 0,'zn_z'] + L
        gal_coods.loc[gal_coods['zn_z'] > L,'zn_z'] = gal_coods.loc[gal_coods['zn_z'] > L,'zn_z'] - L

    print "Building periodic KDtree..."
    # T = PeriodicCKDTree(dimensions, gals[['zn_x','zn_y','zn_z']])
    T = PeriodicCKDTree(dimensions, gal_coods)

    print "Calculating cluster properties..."
    for i, cid in enumerate(clusters['z0_centralId']):
        # sys.stdout.flush()

        no_pcs = np.sum(gals['z0_centralId'] == cid)  # total number of galaxies in protocluster

        # subset protocluster galaxy coordinates
        # coods = gals[gals['z0_centralId'] == cid][['zn_x','zn_y','zn_z']].copy().reset_index(drop=True)
        coods = gal_coods[gals['z0_centralId'] == cid].reset_index(drop=True)

        # normalise coordinate values
        if np.abs(coods['zn_x'].max() - coods['zn_x'].min()) > L/2:
            coods['zn_x'] = coods['zn_x'] - L
            coods.loc[coods['zn_x'] < -L/2, 'zn_x'] += L 

        if np.abs(coods['zn_y'].max() - coods['zn_y'].min()) > L/2:
            coods['zn_y'] = coods['zn_y'] - L
            coods.loc[coods['zn_y'] < -L/2, 'zn_y'] += L 

        if np.abs(coods['zn_z'].max() - coods['zn_z'].min()) > L/2:
            coods['zn_z'] = coods['zn_z'] - L
            coods.loc[coods['zn_z'] < -L/2, 'zn_z'] += L

        center = np.median(coods, axis=0)  # find protocluster center

        completeness = []
        purity = []
        dgal = []

        # calculate statistics over R
        for R in [float(x)/2 for x in range(61)[1:]]:

            gal_index = T.query_ball_point(center, r=R)

            all_gals_in_R = len(gal_index)
            pcs_in_R = float(sum(gals.ix[gal_index]['z0_centralId'] == cid))
            completeness.append(pcs_in_R / no_pcs)

            avg_dgal = float(gals.shape[0]) / L**3 * 4./3 * np.pi * R**3
            dgal.append((all_gals_in_R - avg_dgal) / avg_dgal)

            if all_gals_in_R == 0:
                purity.append(1)
            else:
                purity.append(pcs_in_R / all_gals_in_R)


        cluster_stats[i] = np.array([completeness, purity, dgal])


    # print "Calculating percentiles..."
    # completeness_percentiles = np.array([np.percentile(y, [84,16]) for y in np.vstack([x[0] for x in cluster_stats]).T])
    # purity_percentiles = np.array([np.percentile(y, [84,16]) for y in np.vstack([x[1] for x in cluster_stats]).T])

    return {'cstats': np.array(cluster_stats), 'clusters': clusters}


if __name__ == "__main__":

    print "Reading data..."
    sys.stdout.flush()
    gals_z6p42_sfr = pd.read_csv('data/henriques2015a_z3p95_mstar.csv', skiprows=111, skipfooter=1, engine='python')

    print "Calculating stats..."
    sys.stdout.flush()
    cstats = cluster_stats(gals_z6p42_sfr, L = 480.279)

    # print cstats

    plt.plot(range(31), np.ma.masked_where(np.vstack([x[1] for x in cstats['cstats'] if x[1]])==0,
                    np.vstack([x[1] for x in cstats['cstats'] if x[1]])).mean(axis=0), c=cmap[4], label='6.42')

    plt.plot(range(31), np.mean(np.vstack([x[0] for x in cstats['cstats']]), axis=0), c='blue')

    plt.fill_between(range(31), cstats['completeness_percentiles'][:,0],
                     cstats['completeness_percentiles'][:,1], alpha=0.5, label='completeness')

    plt.fill_between(range(31), cstats['purity_percentiles'][:,0], cstats['purity_percentiles'][:,1],
                     alpha=0.5, color='green', label='purity')

    plt.show()
