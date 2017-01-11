
import numpy as np
import pandas as pd


redshift = '2p07'
selection = 'sfr'


print "Reading data..."
gals = pd.read_csv('data/henriques2015a_z%s_mstar.csv' % redshift, skiprows=104, skipfooter=1, engine='python')

print "Opening file..."
f = open('coods_%s_%s.ascii' % (redshift, selection), 'wb')

print "Writing header..."
f.writelines(['ANDFIELD COORDS\n', '[3,%d]\n' % len(gals), 'BBOX [0,0,0] [2168.64,2168.64,2168.64]\n'])

print "Writing data..."
for i in range(len(gals)):
    f.write('%f %f %f\n' % (gals.ix[i]['zn_x'], gals.ix[i]['zn_y'], gals.ix[i]['zn_z']))


