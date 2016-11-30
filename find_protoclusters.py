##
## Find protoclusters
##

import pandas as pd
import numpy as np

import sys

import pickle as pcl

import hightolowz

redshift = '6p42'
selection = 'mstar9'


print "Reading data..."
gals = pd.read_csv('/lustre/scratch/astro/cl478/protoclusters_data/henriques2015a_z%s_mstar.csv' % redshift, skiprows=104, skipfooter=1, engine='python')


dgal_r20 = pd.read_csv('data/planck1/dgal_mstar9_%s_r20_gals.csv' % redshift)
dgal_r15 = pd.read_csv('data/planck1/dgal_mstar9_%s_r15_gals.csv' % redshift)
dgal_r10 = pd.read_csv('data/planck1/dgal_mstar9_%s_r10_gals.csv' % redshift)
dgal_r5 = pd.read_csv('data/planck1/dgal_mstar9_%s_r5_gals.csv' % redshift)

gals = pd.concat([gals, dgal_r20, dgal_r15, dgal_r10, dgal_r5], axis=1)

# initialise
print "Intialising arrays..."

L = 480.279

dimensions = np.array([L,L,L])


for R in [20,15,10,5]:

    print str(R)

    ignore_gals = np.array([True] * len(gals))
    ignore_dgal = np.array([True] * len(gals))

    protoclusters = []
    pc_members = []

    
    dgal = pd.DataFrame({'index': np.argsort(gals['delta_gal_%s' % str(R)])[::-1].astype(int),
                         'dgal': np.sort(gals['delta_gal_%s' % str(R)])[::-1]})
    
    # print to screen (same lines)
    import curses
    stdscr = curses.initscr()
    curses.noecho()
    #curses.cbreak()
    
    
    while (sum(~ignore_dgal) != len(ignore_dgal)):
    
        # filter dgal by all those points not within 2*R of other protoclusters
        # return the highest dgal available
        temp, max_dgal, max_index = dgal[dgal['index'].isin(np.where(ignore_dgal)[0])].reset_index().ix[0]
    
        # find distance to all galaxies from this overdensity
        dist = np.vstack(hightolowz.distance(gals[['zn_x','zn_y','zn_z']].ix[max_index],
                                             gals[['zn_x','zn_y','zn_z']],
                                             dimensions))[0]
    
    
        if sum((dist < R) & (~ignore_gals)) > 1:
            #print "\rSelection regions contains previously selected galaxies. Skipping."
            ignore_array[max_index] = False
    
        else:
            #print "\rIsolated region. Saving."
            # update ignore array
            ignore_gals[np.array(dist < R)] = False
            ignore_dgal[np.array(dist < 2*R)] = False
    
            # save max index
            protoclusters.extend([max_index.astype(int)])
    
            # save all protocluster members within R
            pc_members.append(np.where(dist < R))
    
        stdscr.addstr(0, 0, "R: %s" % str(R))
        stdscr.addstr(1, 0, "Ignore_dgal: {}".format(round(float(sum(~ignore_dgal))/len(ignore_dgal), 4) * 100))
        stdscr.addstr(2, 0, "delta_gal: {}".format(round(max_dgal, 3)))
        stdscr.addstr(3, 0, "Galaxies matched: {}".format(round(float(sum(~ignore_gals)) / len(ignore_gals),4) * 100, "%"))
        stdscr.addstr(4, 0, "Protoclusters identified: {}".format(len(protoclusters)))
        stdscr.refresh()
        sys.stdout.flush()
    
    
    curses.echo()
    #curses.nocbreak()
    curses.endwin()
    
    pcl.dump([protoclusters, pc_members], open('data/planck1/protoclusters_%s_%s_R%s.p' % (redshift, selection, str(R)), 'wb'))
