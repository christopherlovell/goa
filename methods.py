import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from collections import Counter
from astropy.cosmology import Planck13


def get_protoclusters(gals, L, cluster_lim=1.48e4):
    
    clusters = gals[gals['z0_central_mcrit200'] > cluster_lim].groupby('z0_centralId')['z0_central_mcrit200','z0_centralId'].max()

    coods = np.zeros(len(clusters))
    
#    for i, cid in enumerate(clusters['z0_centralId']):
#        gal_coods = norm_coods(np.array(gals[gals['z0_centralId'] == cid][['zn_x','zn_y','zn_z']]), L)
#         coods[i] = np.mean(gal_coods, axis=0)  # find protocluster center

    coods = np.array([np.mean(normalise_coods(np.array(gals[gals['z0_centralId'] == cid][['zn_x','zn_y','zn_z']]), L), axis=0) for cid in clusters['z0_centralId']])
    
    coods[coods < 0] += L
        
    return pd.DataFrame(coods), clusters


def normalise_coods(coods, L):
    
    original_coods = coods.copy()
    
    for dim in range(3):
        
        if np.abs(coods[:,dim].max() - coods[:,dim].min()) > L/2:
            coods[:,dim] = coods[:,dim] - L
            coods[coods[:,dim] < -L/2,dim] += L #original_coods[coods[:,dim] < -L/2, dim]

    # center = np.median(coods, axis=0)
    # coods = coods - center
    
    return coods


def z_distort(gals, z, L):

    # Convert z-axis to redshift space
    gals['zn_z'] += gals['zn_velZ'] * (1+z) / Planck13.H(z).value
    
    # fix z coordinate positions if they fall outside the box
    gals.loc[gals['zn_z'] < 0,'zn_z'] = gals.loc[gals['zn_z'] < 0,'zn_z'] + L
    gals.loc[gals['zn_z'] > L,'zn_z'] = gals.loc[gals['zn_z'] > L,'zn_z'] - L

    return gals


def factor_h(gals, h):

    gals[['z0_central_mcrit200','zn_x','zn_y','zn_z']] /= h

    return gals


def bhattacharyya(histA, histB):
    """
    Calculate the Bhattacharyya distance between two PDFs

    Args:
        histA (list): normalised PDF values
        histB (list): same length as histA
    """

    if(len(histA) != len(histB)):
        print "\nDistributions are not the same length."
        return -1

    BC = 0;
    for i in range(len(histA)):
        BC += np.sqrt( histA[i] * histB[i] );

    distance = -np.log(BC)
    #angle = math.acos(BC)
    #score = np.sqrt(1 - (1 / np.sqrt(histA.mean() * histB.mean() * len(histA)**2)) * BC)

    return distance, BC


def distance(x0, x1, dimensions):
    """
    Calculate distance on periodic grid
    Source: http://stackoverflow.com/questions/11108869/optimizing-python-distance-calculation-while-accounting-for-periodic-boundary-co

    Args:
        x0 (numpy array): array of coordinates
        x1 (numpy array): array of coordinates
        dimensions (numpy array): length along each dimension
    """

    delta = np.abs(x0 - x1)
    delta = np.where(delta > 0.5 * dimensions, delta - dimensions, delta)

    return [np.linalg.norm(delta, axis=-1)]
    #return [np.sqrt((delta ** 2).sum(axis=-1))]


def z0_halo_properties(pcs, pcmems, gals):
    """
    Calculate descendant details for protocluster candidates

    Args:
        pcs (numpy array): protocluster ids
        pcmems (list of arrays): ids of each protocluster member
        gals (pandas dataframe): high-z galaxy properties
    """

    z0_halos = [None] * len(pcmems)
    halo_ratio = [None] * len(pcmems)

    for i in range(len(pcmems)):
        z0_halos[i] = Counter(gals.ix[pcmems[i]]['z0_centralId'])

        halo_ratio[i] = round(float(z0_halos[i].most_common()[0][1]) \
                  / sum([x[1] for x in z0_halos[i].most_common()]) * 100, 2)

    return z0_halos, halo_ratio


def label(stats, clim, plim, mlim=5e4):
    
    completeness = stats[:,1]
    purity = stats[:,2]
    mass = stats[:,3]
    
    labels = ['proto_lomass','proto_himass','part_lomass','part_himass','pfield_lomass','pfield_himass','field']
    
    # initialise empty label array
    labs = np.array([None] * stats.shape[0])
    
    # assign labels, split by configuration and descendant mass
    labs[(completeness >= clim) & (purity >= plim) & (mass < mlim)] = 'proto_lomass'
    labs[(completeness >= clim) & (purity >= plim) & (mass >= mlim)] = 'proto_himass'
    labs[(completeness < clim) & (purity >= plim) & (mass >= mlim)] = 'part_lomass'
    labs[(completeness < clim) & (purity >= plim) & (mass >= mlim)] = 'part_himass'
    labs[(completeness >= clim) & (purity < plim) & (mass < mlim)] = 'pfield_lomass'
    labs[(completeness >= clim) & (purity < plim) & (mass >= mlim)] = 'pfield_himass'
    labs[(completeness < clim) & (purity < plim)] = 'field'
    
    return labs, labels


def binit(stats, labs, labels, N = 12):
    """
    initialise bins and limits, calculate binned statistics
    
    Args:
        stats: completeness and purity statistics 
        labs: list of labelled regions
        labels: list of label strings
    """

    dgal = stats[:,0] + 1

    binLimits = np.linspace(0, int(np.max(dgal)+1), N)

    lower_bin = binLimits[1] + (binLimits[0]-binLimits[1])/2. 
    upper_bin = binLimits[-1] + (binLimits[0]-binLimits[1])/2.

    bins = np.linspace(lower_bin, upper_bin, N-1)
    
    agg = {x: np.histogram(dgal[labs==x], binLimits)[0] for x in labels}  # save counts for each label
    # agg_total = {x: np.sum(agg[x]) for x in agg}  # find total in each bin
    agg_total = np.sum([v for k,v in agg.iteritems()],axis=0)
    
    fracs = {k: v.astype(float) / agg_total for k,v in agg.iteritems()}

    # fracs = {k: find_fracs(v.astype(float), agg_total) for k,v in agg.iteritems()}  # find fraction of each label
    
    return bins, binLimits, agg, agg_total, fracs


def plotit(ax, stats, axb=None, clim = 0.5, plim = 0.5, N = 12, mlim=5e4, noplot=False):
    """
    
    Args:
        ax - axis object
        selection - selection object
        rid - id of radius selection in selection object
        zis - if od redshift selection in selection object
        axb - bottom axis object. If None, only plots 
    """
    
    colors = ['dimgrey','lightseagreen','lightcoral', 'y']
    
    dgal = stats[:,0] + 1
    # completeness = stats[:,1]
    # purity = stats[:,2]
    mass = stats[:,3]

    labels = ['proto_lomass','proto_himass','part_lomass','part_himass','pfield_lomass','pfield_himass','field']
    labs = label(stats, clim, plim, mlim)
    
    # initialise bins and limits
    binLimits = np.linspace(0, int(np.max(dgal)+1), N)
    
    # print binLimits

    lower_bin = binLimits[1] + (binLimits[0]-binLimits[1])/2. 
    upper_bin = binLimits[-1] + (binLimits[0]-binLimits[1])/2.

    bins = np.linspace(lower_bin, upper_bin, N-1)
    
    # save counts for each label
    agg = {x: np.histogram(dgal[labs==x], binLimits)[0] for x in labels}
    
    agg_total = np.sum(np.vstack([agg[x] for x in agg]), axis=0).astype(float)
    
    # truncate range to where there are at least a couple of samples
    n_limit = 1
    if (np.sum(agg_total < n_limit) > 0):
        
        mask = range(0,np.min(np.where(np.sum(agg,axis=0) < n_limit)))
        
        bins = bins[mask]
        binLimits = binLimits[range(0, np.max(mask)+2)]
        agg_total = agg_total[mask]
    
        for i in range(len(labels)):
            agg[i] = agg[i][mask]
        

    
    # probability density function
    if axb != None:
        
        phiMax = 0.
        
        mask = (label==11) | (label==12) | (label==21) | (label==22)
        phiA = np.histogram(dgal[mask], binLimits, normed=True)[0]

        mask = (label==0)
        phiB = np.histogram(dgal[mask], binLimits, normed=True)[0]
        
        phiMax = np.max(phiA)
        
        DB, BC = bhattacharyya(phiA*np.diff(binLimits), phiB*np.diff(binLimits))
        
        
        mask = (label==12) | (label==22)
        phiC = np.histogram(dgal[mask], binLimits, normed=True)[0]
        
        DB_himass, BC_himass = bhattacharyya(phiC*np.diff(binLimits), phiB*np.diff(binLimits))

        if noplot:
            print "DB(All), DB(High mass)"
            return round(DB, 2), round(DB_himass, 2)
            exit
        
        #axb.text(0.6,0.8, '$D_{B} = %s$'%round(DB, 3), transform=axb.transAxes)

        # axb.step(bins, phiA, color=colors[1], linestyle='dashed')
        axb.step(bins, phiB, color=colors[0], linestyle='solid', linewidth=3)
        
        mask = (label==11)
        if np.sum(mask) > 5:
            phi = np.histogram(dgal[mask], binLimits, normed=True)[0]
            axb.step(bins, phi, color=colors[1], linestyle='solid')
            phiMax = np.max([phiMax, np.max(phi)])
        
        mask = (label==12)
        if np.sum(mask) > 5:
            phi = np.histogram(dgal[mask], binLimits, normed=True)[0]
            axb.step(bins, phi, color=colors[1], linestyle='dashed')
            phiMax = np.max([phiMax, np.max(phi)])
        
        mask = (label==21)
        if np.sum(mask) > 5:
            phi = np.histogram(dgal[mask], binLimits, normed=True)[0]
            axb.step(bins, phi, color=colors[3], linestyle='solid')
            phiMax = np.max([phiMax, np.max(phi)])
        
        mask = (label==22)
        if np.sum(mask) > 5:
            phi = np.histogram(dgal[mask], binLimits, normed=True)[0]
            axb.step(bins, phi, color=colors[3], linestyle='dashed')
            phiMax = np.max([phiMax, np.max(phi)])
            
            
        axb.step(bins, phiA, color=colors[1], linestyle='solid', linewidth=3)
        axb.step(bins, phiA, color=colors[3], linestyle='dotted', linewidth=3)
        
        axb.set_ylim(0, phiMax + 0.1)
    
    
    width = binLimits[2] - binLimits[1]
    
    plt.rcParams['hatch.color'] = 'black'
    plt.rcParams['hatch.linewidth'] = 0.5
    
    ax.bar(bins, agg[1] / agg_total, width=width, 
           label='protocluster ($M<M_{lim}$)', alpha=0.6, color=colors[1])
    
    bar = ax.bar(bins, agg[2] / agg_total, width=width, bottom=agg[1] / agg_total,  
           label='protocluster ($M>M_{lim}$)', alpha=0.8, color=colors[1], hatch='///')

    ax.bar(bins, agg[3] / agg_total, width=width, bottom=np.sum(agg[1:3],axis=0) / agg_total, 
           label='part of a \n protocluster', alpha=0.6, color=colors[3]) 
    
    ax.bar(bins, agg[4] / agg_total, width=width, bottom=np.sum(agg[1:4],axis=0) / agg_total,
           label='part of a \n protocluster', alpha=0.8, color=colors[3], hatch='///') 

    ax.bar(bins, agg[5] / agg_total, width=width, bottom= np.sum(agg[1:5],axis=0) / agg_total,
           label='protocluster \n + field', alpha=0.6, color=colors[2])
    
    ax.bar(bins, agg[6] / agg_total, width=width, bottom= np.sum(agg[1:6],axis=0) / agg_total,
           label='protocluster \n + field', alpha=0.8, color=colors[2], hatch='///')

    ax.bar(bins, agg[0] / agg_total, width=width, 
           bottom= np.sum(agg[1:],axis=0) / agg_total, color=colors[0], 
           label='field', alpha=0.2)
    
    ax.set_xlim(binLimits[0], binLimits[-1])
    if axb:
        axb.set_xlim(binLimits[0], binLimits[-1])
