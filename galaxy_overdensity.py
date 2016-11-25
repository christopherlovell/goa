"""
Calculate galaxy overdensity
- can calculate for random regions, or for each galaxy. If random, set random flag to 'True'.
"""

random = False

import pandas as pd
import numpy as np
import itertools as it
from scipy.spatial.distance import cdist

from hightolowz import distance


print "Reading galaxy data..."

gals = pd.read_csv('data/planck1/henriques2015a_z3p95_mstar.csv', skiprows=101, skipfooter=1, engine='python')

selection_str = 'mstar9'
redshift_str = '3p95'

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

r = [20, 15, 10, 5]
r_str = ['20', '15', '10', '5']

ngal = {'20': [], '15': [], '10': [], '5': []}
dgal = {'20': [], '15': [], '10': [], '5': []}
max_fraction = {'20': [], '15': [], '10': [], '5': []}
max_fraction_mass = {'20': [], '15': [], '10': [], '5': []}
n_cluster_desc = {'20': [], '15': [], '10': [], '5': []}


print "Counting galaxies..."

n = 100

# can't calculate distances all in one go, so need to chunk
#for i,gals in z6_galaxies_mstar.groupby(np.arange(len(z6_galaxies_mstar))//n):
for i,c in coods.groupby(np.arange(len(coods))//n):

    if i % 5 == 0:
        print round(float(c.shape[0] * (i+1)) / coods.shape[0] * 100, 2), '%'

    # calculate distances
    dist = np.vstack(c.apply(lambda x: distance(x, gals[['zn_x','zn_y','zn_z']], dimensions), axis=1))

    for R, R_str in zip(r, r_str):

        gal_index = dist < R

        for i in range(len(gal_index)):

            n_galaxy = np.sum(gal_index[i])
            ngal[R_str].append(n_galaxy)

            if n_galaxy == 0:
                m_max=0
                ncd=0
                agg_count=np.array([0])
                max_frac=0
            else:
                agg_mvir = gals.ix[gal_index[i]].groupby('z0_centralId').mean()['z0_central_mcrit200']
                agg_count = gals.ix[gal_index[i]].groupby('z0_centralId')['z0_centralId'].count().astype(float)

                agg = pd.DataFrame([agg_mvir, agg_count]).T

                m_max = agg.loc[agg['z0_centralId'].idxmax()]['z0_central_mcrit200'] # find mass of most common descendant
                ncd = sum(agg_mvir > 1e4)
                max_frac = agg_count.max() / agg_count.sum()

            max_fraction_mass[R_str].append(m_max)
            max_fraction[R_str].append(max_frac)
            n_cluster_desc[R_str].append(ncd)


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
