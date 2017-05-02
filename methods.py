import numpy as np

def norm_coods(coods, L):
    
    original_coods = coods.copy()
    
    for dim in range(3):
        
        if np.abs(coods[:,dim].max() - coods[:,dim].min()) > L/2:
            coods[:,dim] = coods[:,dim] - L
            coods[coods[:,dim] < -L/2,dim] += L #original_coods[coods[:,dim] < -L/2, dim]

    # center = np.median(coods, axis=0)
    # coods = coods - center
    
    return coods

from collections import Counter


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


import matplotlib.pyplot as plt

def plotit(ax, selection, rid, zid, axb=None, clim = 0.5, plim = 0.5, N = 12, noplot=False):
    """
    
    Args:
        ax - axis object
        selection - selection object
        rid - id of radius selection in selection object
        zis - if od redshift selection in selection object
        axb - bottom axis object. If None, only plots 
    
    """
    colors = ['dimgrey','lightseagreen','lightcoral', 'y']
    
    dgal = selection[rid,zid,:,0] + 1
    completeness = selection[rid,zid,:,1]
    purity = selection[rid,zid,:,2]
    mass = selection[rid,zid,:,3]

    # initialise empty label array
    label = np.zeros(selection.shape[2])    
    
    labels = [0,11,12,21,22,31,32]
    
    mlim = 5e4

    # assign labels. split by configuration and descendant mass
    label[(completeness >= clim) & (purity >= plim) & (mass < mlim)] = labels[1]   # protocluster, M < M_lim
    label[(completeness >= clim) & (purity >= plim) & (mass >= mlim)] = labels[2]  # protocluster, M >= M_lim
    label[(completeness < clim) & (purity >= plim) & (mass >= mlim)] = labels[3]   # part protocluster, M < M_lim
    label[(completeness < clim) & (purity >= plim) & (mass >= mlim)] = labels[4]   # part protocluster, M >= M_lim
    label[(completeness >= clim) & (purity < plim) & (mass < mlim)] = labels[5]    # protocluster+field, M < M_lim
    label[(completeness >= clim) & (purity < plim) & (mass >= mlim)] = labels[6]   # protocluster+field, M >= M_lim
    label[(completeness < clim) & (purity < plim)] = labels[0]     # field
    
    # initialise bins and limits
    binLimits = np.linspace(0, int(np.max(dgal)+1), N)

    lower_bin = binLimits[1] + (binLimits[0]-binLimits[1])/2. 
    upper_bin = binLimits[-1] + (binLimits[0]-binLimits[1])/2.

    bins = np.linspace(lower_bin, upper_bin, N-1)
    
    # save counts for each label
    agg = [[] for x in labels]
    
    # agg[0] = np.histogram(dgal, binLimits)[0]    # total count
    
    for i, x in enumerate(labels):
        agg[i] = np.histogram(dgal[label==x], binLimits)[0]
    
    agg_total = np.sum(agg,axis=0).astype(float)
    
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
