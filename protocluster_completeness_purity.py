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

def cluster_stats(gals, L):

    cluster_ids = pd.unique(gals[gals['z0_central_mcrit200'] > 1e4]['z0_centralId'])

    print "Clusters: ", len(cluster_ids)

    cluster_stats = [None] * len(cluster_ids)

    dimensions = np.array([L, L, L])

    print "Building periodic KDtree..."
    T = PeriodicCKDTree(dimensions, gals[['zn_x','zn_y','zn_z']])

    for i, cluster in enumerate(cluster_ids):

        # sys.stdout.flush()

        coods = gals[gals['z0_centralId']==cluster][['zn_x','zn_y','zn_z']].copy()
        coods = coods.reset_index(drop=True)

        if np.abs(coods['zn_x'].max() - coods['zn_x'].min()) > L/2:
            coods['zn_x'] = coods['zn_x'] - L
            coods.loc[coods['zn_x'] < -L/2, 'zn_x'] = gals[(gals['z0_centralId'] == cluster) & (coods['zn_x'] < -L/2)]['zn_x']

        if np.abs(coods['zn_y'].max() - coods['zn_y'].min()) > L/2:
            coods['zn_y'] = coods['zn_y'] - L
            coods.loc[coods['zn_y'] < -L/2, 'zn_y'] = gals[(gals['z0_centralId'] == cluster) & (coods['zn_y'] < -L/2)]['zn_y']

        if np.abs(coods['zn_z'].max() - coods['zn_z'].min()) > L/2:
            coods['zn_z'] = coods['zn_z'] - L
            coods.loc[coods['zn_z'] < -L/2, 'zn_z'] = gals[(gals['z0_centralId'] == cluster) & (coods['zn_z'] < -L/2)]['zn_z']

        center = np.mean(coods)
        
        # gal_dist = hightolowz.distance(center, gals[['zn_x','zn_y','zn_z']], dimensions)[0]

        no_pcs = np.sum(gals['z0_centralId'] == cluster)

        completeness = []
        purity = []

        for R in [float(x)/2 for x in range(61)]:

            gal_index = T.query_ball_point(center, r=R)

            all_gals_in_R = len(gal_index)
            pcs_in_R = float(sum(gals.ix[gal_index]['z0_centralId'] == cluster))
            completeness.append(pcs_in_R / no_pcs)

            # all_gals_in_R = np.sum(gal_dist < R)
            # pcs_in_R = float(np.sum(gal_dist[np.array(gals['z0_centralId'] == cluster)] < R))
            # completeness.append(pcs_in_R / no_pcs)

            if all_gals_in_R == 0:
                purity.append(1)
            else:
                purity.append(pcs_in_R / all_gals_in_R)


        cluster_stats[i] = [completeness, purity]

    completeness_percentiles = np.array([np.percentile(y, [90,10]) for y in np.vstack([x[0] for x in cluster_stats]).T])
    purity_percentiles = np.array([np.percentile(y, [90,10]) for y in np.vstack([x[1] for x in cluster_stats]).T])


    return {'cstats': cluster_stats,
            'completeness_percentiles': completeness_percentiles,
            'purity_percentiles': purity_percentiles}


if __name__ == "__main__":

    print "Reading data..."
    sys.stdout.flush()
    gals_z6p42_sfr = pd.read_csv('data/henriques2015a_z3p95_mstar.csv', skiprows=104, skipfooter=1, engine='python')

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
