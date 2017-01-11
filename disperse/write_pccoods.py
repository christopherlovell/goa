
import numpy as np
import pandas as pd
from itertools import product

redshift = ['2p07','3p10','3p95','5p03','6p42']
selection = ['sfr','mstar']

# header is a different size for each file...
len_header = [111] + [104] * 9

for (z,select),lhead in zip(product(redshift, selection), len_header):

    print "Redshift:",z
    print "Selection:",select

    print "Reading data..."
    gals = pd.read_csv('../data/henriques2015a_z%s_%s.csv' % (z,select), 
                        skiprows=lhead, skipfooter=1, engine='python')
    
    print "Looping through protoclusters..."
    for i in range(len(np.unique(gals[gals['z0_central_mcrit200'] > 1e4]['z0_centralId']))):
    
#         print "pc:",i
    
        index = np.int(np.unique(gals[gals['z0_central_mcrit200'] > 1e4]['z0_centralId'])[i])
        coods = gals.loc[(gals['z0_central_mcrit200'] > 1e4) & (gals['z0_centralId'] == index)][['zn_x','zn_y','zn_z']]
    
        header = 'ANDFIELD  COORDS\n[3, %s]\n' % coods.shape[0]
    
        with open('../data/disperse/pccoods_%s_%s_%04d.ascii' % (z, select, i), 'w') as outfile:
    
            outfile.write(header)
    
            coods.to_csv(outfile, sep=' ', index=False, header=None)
    





