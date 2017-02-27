##
## Find protoclusters
##

import pandas as pd
import numpy as np

import sys

import pickle as pcl

import hightolowz


# redshift = '6p42'
redshift = sys.argv[1]
# selection = 'mstar9'
selection = sys.argv[2]

print "z:", redshift, ", selection:", selection
sys.stdout.flush()

directory = '/lustre/scratch/astro/cl478/protoclusters_data'

print "Reading data..."
sys.stdout.flush()

gals = pd.read_csv('%s/henriques2015a_z%s_mstar.csv' % (directory, redshift), skiprows=104, skipfooter=1, engine='python')


dgal_r20 = pd.read_csv('%s/dgal_%s_%s_r20_gals.csv' % (directory, selection, redshift))
dgal_r12p5 = pd.read_csv('%s/dgal_%s_%s_r12.5_gals.csv' % (directory, selection, redshift))
dgal_r7p5 = pd.read_csv('%s/dgal_%s_%s_r7.5_gals.csv' % (directory, selection, redshift))
dgal_r5 = pd.read_csv('%s/dgal_%s_%s_r5_gals.csv' % (directory, selection, redshift))

gals = pd.concat([gals, dgal_r20, dgal_r12p5, dgal_r7p5, dgal_r5], axis=1)

# initialise
print "Intialising arrays..."

L = 480.279

dimensions = np.array([L,L,L])


for R in [20,12.5,7.5,5]:

    print str(R)
    sys.stdout.flush()

    ignore_gals = np.array([True] * len(gals))
    ignore_dgal = np.array([True] * len(gals))

    protoclusters = []
    pc_members = []

    
    dgal = pd.DataFrame({'index': np.argsort(gals['delta_gal_%s' % str(R)])[::-1].astype(int),
                         'dgal': np.sort(gals['delta_gal_%s' % str(R)])[::-1]})
    
    max_dgal = 1    
    i = 0
    while ((sum(~ignore_dgal) != len(ignore_dgal)) & (max_dgal > 0.) & (i < 10000)):

        i += 1

        sys.stdout.flush()
        if i % 20 == 0:
            print "R: ", R, " ", i
    
        # filter dgal by all those points not within 2*R of other protoclusters
        # return the highest dgal available
        temp, max_dgal, max_index = dgal[dgal['index'].isin(np.where(ignore_dgal)[0])].reset_index().ix[0]
    
        # find distance to all galaxies from this overdensity
        dist = np.vstack(hightolowz.distance(gals[['zn_x','zn_y','zn_z']].ix[max_index],
                                             gals[['zn_x','zn_y','zn_z']],
                                             dimensions))[0]
    
    
        if sum((dist < R) & (~ignore_gals)) > 1:
            # Selection regions contains previously selected galaxies. Skipping.
            ignore_array[max_index] = False
    
        else:
            # Isolated region. Saving. 
            # Update ignore array
            ignore_gals[np.array(dist < R)] = False
            ignore_dgal[np.array(dist < 2*R)] = False
    
            # save max index
            protoclusters.extend([max_index.astype(int)])
    
            # save all protocluster members within R
            pc_members.append(np.where(dist < R))
    
#        stdscr.addstr(0, 0, "R: %s" % str(R))
#        stdscr.addstr(1, 0, "Ignore_dgal: {}".format(round(float(sum(~ignore_dgal))/len(ignore_dgal), 4) * 100))
#        stdscr.addstr(2, 0, "delta_gal: {}".format(round(max_dgal, 3)))
#        stdscr.addstr(3, 0, "Galaxies matched: {}".format(round(float(sum(~ignore_gals)) / len(ignore_gals),4) * 100, "%"))
#        stdscr.addstr(4, 0, "Protoclusters identified: {}".format(len(protoclusters)))
#        stdscr.refresh()
        # print round(float(sum(~ignore_gals)) / len(ignore_gals),4) * 100
        # sys.stdout.flush()
    
    
    pcl.dump([protoclusters, pc_members], open('data/planck1/protoclusters_%s_%s_R%s.p' % (redshift, selection, str(R)), 'wb'))


