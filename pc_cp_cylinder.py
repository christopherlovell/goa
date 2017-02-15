"""
Calculate protocluster completeness and purity at a given radii for different selections

"""

import matplotlib.pyplot as plt

from palettable.tableau import GreenOrange_12

cmap = GreenOrange_12.hex_colors

import pandas as pd
import numpy as np

import sys

# https://github.com/patvarilly/periodic_kdtree
from periodic_kdtree import PeriodicCKDTree

# custom cython coordinate filter
from norm_cood import norm_coods

def cluster_stats(gals, L=500):
    """
    Calculate overdensity, completeness and purity of the protocluster
    population in volumes with varying raddi and delta-z.

    NOTE: need to convert z-coordinate to redshift space *before* calling

    Args:
        - gals: dataframe of galaxy properties
        - L: comoving side length of simulation
    """

    # r = [2.5, 5., 7.5, 10., 15.]
    r = [7.5]
    # deltaz = [2.5, 5., 7.5, 10., 15.]
    deltaz = [7.5]

    avg = float(gals.shape[0]) / L**3

    # filter by *group* masses
    clusters = gals[gals['z0_central_mcrit200'] > 1e3].groupby('z0_centralId')['z0_central_mcrit200','z0_centralId'].max()

    print "N Clusters: ", len(clusters)

    cluster_stats = [None] * len(clusters)

    dimensions = np.array([L, L, L])

    print "Building periodic KDtree..."
    T = PeriodicCKDTree(dimensions, gals[['zn_x','zn_y','zn_z']])

    cstats = np.zeros((len(clusters), len(deltaz), len(r), 3))

    print "Calculating cluster properties..."
    for i, cid in enumerate(clusters['z0_centralId']):

        no_pcs = np.sum(gals['z0_centralId'] == cid)  # total number of galaxies in protocluster

        # subset protocluster galaxy coordinates
        coods = gals[gals['z0_centralId'] == cid][['zn_x','zn_y','zn_z']].copy().reset_index(drop=True)

        # normalise coordinate values
        if np.abs(coods['zn_x'].max() - coods['zn_x'].min()) > L/2:
            coods['zn_x'] = coods['zn_x'] - L
            coods.loc[coods['zn_x'] < -L/2, 'zn_x'] = gals[(gals['z0_centralId'] == cid) & (coods['zn_x'] < -L/2)]['zn_x']

        if np.abs(coods['zn_y'].max() - coods['zn_y'].min()) > L/2:
            coods['zn_y'] = coods['zn_y'] - L
            coods.loc[coods['zn_y'] < -L/2, 'zn_y'] = gals[(gals['z0_centralId'] == cid) & (coods['zn_y'] < -L/2)]['zn_y']

        if np.abs(coods['zn_z'].max() - coods['zn_z'].min()) > L/2:
            coods['zn_z'] = coods['zn_z'] - L
            coods.loc[coods['zn_z'] < -L/2, 'zn_z'] = gals[(gals['z0_centralId'] == cid) & (coods['zn_z'] < -L/2)]['zn_z']


        center = np.mean(coods)  # find protocluster center

        # get all galaxies in sphere of radii the max extent of the largest cylinder
        gal_index = T.query_ball_point(center, r=(max(r)**2 + max(deltaz)**2)**0.5)

        for j, dz in enumerate(deltaz):

            for k, R in enumerate(r):

                # filter by cylinder using norm_coods()
                gal_index_temp = np.array(gal_index)[norm_coods(gals.ix[gal_index][['zn_x','zn_y','zn_z']].values, center.values, R, dz, L)]

                all_gals_in_R = len(gal_index_temp)
                pcs_in_R = float(sum(gals.ix[gal_index_temp]['z0_centralId'] == cid))
                cstats[i,j,k,1] = pcs_in_R / no_pcs  # completeness

                avg_dgal = avg * np.pi * R**2 * dz
                cstats[i,j,k,0] = (all_gals_in_R - avg_dgal) / avg_dgal  # dgal

                # purity
                if all_gals_in_R == 0:
                    cstats[i,j,k,2] = 1
                else:
                    cstats[i,j,k,2] = pcs_in_R / all_gals_in_R


    return {'cstats': cstats, 'clusters': clusters}


if __name__ == "__main__":

    print "Reading data..."
    sys.stdout.flush()
    gals_z6p42_sfr = pd.read_csv('data/henriques2015a_z3p95_mstar.csv', skiprows=111, skipfooter=1, engine='python')

    print "Calculating stats..."
    sys.stdout.flush()
    cstats = cluster_stats(gals_z6p42_sfr, L = 480.279)

    plt.plot(range(31), np.ma.masked_where(np.vstack([x[1] for x in cstats['cstats'] if x[1]])==0,
                    np.vstack([x[1] for x in cstats['cstats'] if x[1]])).mean(axis=0), c=cmap[4], label='6.42')

    plt.plot(range(31), np.mean(np.vstack([x[0] for x in cstats['cstats']]), axis=0), c='blue')

    plt.fill_between(range(31), cstats['completeness_percentiles'][:,0],
                     cstats['completeness_percentiles'][:,1], alpha=0.5, label='completeness')

    plt.fill_between(range(31), cstats['purity_percentiles'][:,0], cstats['purity_percentiles'][:,1],
                     alpha=0.5, color='green', label='purity')

    plt.show()
