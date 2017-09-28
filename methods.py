import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from collections import Counter

from periodic_kdtree import PeriodicCKDTree  # https://github.com/patvarilly/periodic_kdtree
from astropy.cosmology import Planck13

from norm_cood import norm_coods  # cythonised normalisation function



def overdensity_cylinder(gals, coods, R, dc, L, pc_stats=False, cluster_mass_lim=1e4, n=100, verbose=False):
    """
    Find overdensity statistics over the whole simulation box for cylindrical apertures.

    Args:
        gals - dataframe of galaxy properties
        random - bool, whether to use random positions or center on all galaxies
        R - aperture radius, cMpc
        dc - half aperture depth, cMpc
        cluster_mass_lim - limiting descendant mass above which to classify clusters, z0_central_mcrit200
        n - chunk length
        N - number of random regions
    Returns:
        out_stats - output statistics, numpy array of shape [len(coods), 4], where coods is either the number of galaxies or number of random regions.
    """

    dimensions = np.array([L, L, L])

    if verbose: print "Building KDtree..."
    T = PeriodicCKDTree(dimensions, gals[['zn_x','zn_y','zn_z']])
    
    avg = float(gals.shape[0]) / L**3 # average overdensity cMpc^-3
    
    out_stats = np.zeros((len(coods),4))
    
    vol_avg = np.pi * R**2 * (2*dc) * avg  # average overdensity in chosen volume
    
    for j,c in coods.groupby(np.arange(len(coods))//n): # can't calculate distances all in one go, so need to chunk
    
        if verbose: # print progress
            if j % 100 == 0:
                print round(float(c.shape[0] * (j+1)) / coods.shape[0] * 100, 2), '%'
                sys.stdout.flush()

    
        # find all galaxies within a sphere of radius the max extent of the cylinder
        gal_index = T.query_ball_point(c, r=(R**2 + dc**2)**0.5)
    
        # filter by cylinder using norm_coods()
        gal_index = [np.array(gal_index[k])[norm_coods(gals.iloc[gal_index[k]][['zn_x','zn_y','zn_z']].values, c.ix[k + j*n].values, R=R, half_deltac=dc, L=L)] for k in range(len(c))]
    
        start_index = (j*n)  # save start index

        # calculate dgal
        out_stats[start_index:(start_index+len(c)), 0] = (np.array([len(x) for x in gal_index]) - vol_avg) / vol_avg
    
        if pc_stats:  # calculate completeness and purity statistics
    
            for i in range(len(gal_index)):

                cluster_ids = gals.iloc[gal_index[i]]
                cluster_ids = Counter(cluster_ids[cluster_ids['z0_central_mcrit200'] > cluster_mass_lim]['z0_centralId'])

                if len(cluster_ids) > 0:

                    cstats = np.zeros((len(cluster_ids), 2))

                    for k, (q, no) in enumerate(cluster_ids.items()):
                        cluster_gals = gals.ix[gals['z0_centralId'] == q]
                        cstats[k,0] = float(no) / len(cluster_gals)  # completeness
                        cstats[k,1] = float(no) / len(gal_index[i])  # purity


                    # find id of max completeness and purity in cstats array
                    max_completeness = np.where(cstats[:,0] == cstats[:,0].max())[0]
                    max_purity = np.where(cstats[:,1] == cstats[:,1].max())[0]

                    # sometimes multiple clusters can have same completeness or purity in a single candidate
                    # - use the cluster with the highest complementary completeness/purity
                    if len(max_completeness) > 1:

                        # get matches between completeness and purity
                        matches = [x in max_purity for x in max_completeness]

                        if np.sum(matches) > 0:
                            # just use the first one
                            max_completeness = [np.where(matches)[0][0]]
                            max_purity = [np.where(matches)[0][0]]
                        else:
                            max_completeness = [max_completeness[np.argmax(cstats[max_completeness, 1])]]

                    if len(max_purity) > 1:

                        matches = [x in max_completeness for x in max_purity]

                        if np.sum(matches) > 0:
                            max_completeness = [np.where(matches)[0][0]]
                            max_purity = [np.where(matches)[0][0]]

                        else:
                            max_purity = [max_purity[np.argmax(cstats[max_completeness, 0])]]


                    # sometimes the cluster with the highest completeness does not have the highest purity, or vice versa
                    # - use the cluster with the highest combined purity/completeness added in quadrature
                    if max_completeness[0] != max_purity[0]:
                        max_completeness = [np.argmax([pow(np.sum(x**2), 0.5) for x in cstats])]
                        max_purity = max_completeness

                    # save completeness and purity values
                    out_stats[start_index+i, 1] = cstats[max_completeness[0], 0]  # completeness
                    out_stats[start_index+i, 2] = cstats[max_purity[0], 1]        # purity

                    # save descendant mass
                    # filter by cluster id, save z0 halo mass
                    # can use either max_completeness or max_purity, both equal by this point

                    out_stats[start_index+i, 3] = gals.loc[gals['z0_centralId'] == cluster_ids.keys()[max_completeness[0]], 'z0_central_mcrit200'].iloc[0]


                else:  # if no galaxies in aperture
                    out_stats[start_index+i, 1] = 0.
                    out_stats[start_index+i, 2] = 0.
                    out_stats[start_index+i, 3] = np.nan
    

    return out_stats



def fit_func(x, a, b, c, C):
    """
    Fitting function for descendant cluster mass
    
    Args:
        x: (n, 2) array, where first column is galaxy overdensity, second is redshift
        a,b,c,C: parameters
        
    Returns:
        (n, 1) array of estimated descendant mass
        
    """
    return a * (x[1]+1)**b * (x[0]+1)**c + C


def theoretical_ratios(overdensity = False, radii = np.arange(1,50,0.1), a=12.33, b=8.60, c=6.23, N = 200000, Rmax=50):
    """
    Calculate theoretical completeness and purity curves for an arbitrary ellipse with a given overdensity.
    """

    cood_max=np.sqrt(Rmax**2 / 2)

    # uniform distribution
    coordinates = np.reshape(np.random.rand(N*3), (N,3)) * cood_max
    
    in_ellipse = np.array([(coods[0]/a)**2 + (coods[1]/b)**2 + (coods[2]/c)**2 for coods in coordinates])
    
    if overdensity:
        coords_in_ellipse = coordinates[in_ellipse <=1]
        
        for r in range(overdensity):
            coordinates = np.vstack([coordinates, coords_in_ellipse])
            
        in_ellipse = np.array([(coods[0]/a)**2 + (coods[1]/b)**2 + (coods[2]/c)**2 for coods in coordinates])
        

    R = np.array([np.sum(coods**2)**0.5 for coods in coordinates])

    P = [float(np.sum((R <= r) & (in_ellipse <= 1))) / np.sum(R <= r) for r in radii]
    C = [float(np.sum((R <= r) & (in_ellipse <= 1))) / np.sum(in_ellipse <= 1) for r in radii]

    r_sphere = (a*b*c)**(1./3)

    P_sphere = [float(np.sum((R <= r) & (R <= r_sphere))) / np.sum(R <= r) for r in radii]
    C_sphere = [float(np.sum((R <= r) & (R <= r_sphere))) / np.sum(R <= r_sphere) for r in radii]
    
    return P, C, P_sphere, C_sphere

    

def r2(ydata, ypred):
    residuals = ydata - ypred
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((ydata-np.mean(ydata))**2)

    return 1 - (ss_res / ss_tot)


def get_protoclusters(gals, L, cluster_lim=1e4):
    
    clusters = gals[gals['z0_central_mcrit200'] > cluster_lim].groupby('z0_centralId')['z0_central_mcrit200','z0_centralId'].max()

    coods = np.zeros(len(clusters))
    
#    for i, cid in enumerate(clusters['z0_centralId']):
#        gal_coods = norm_coods(np.array(gals[gals['z0_centralId'] == cid][['zn_x','zn_y','zn_z']]), L)
#         coods[i] = np.mean(gal_coods, axis=0)  # find protocluster center

    coods = np.array([np.mean(normalise_coods(np.array(gals[gals['z0_centralId'] == cid][['zn_x','zn_y','zn_z']]), L), axis=0) for cid in clusters['z0_centralId']])
    
    coods[coods < 0] += L
        
    return {'coods': pd.DataFrame(coods), 'clusters': clusters}


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

    gals[['z0_central_mcrit200','zn_x','zn_y','zn_z','z0_x','z0_y','z0_z',
          'z0_mcrit200','z0_central_rcrit200','z0_central_x', 'z0_central_y', 
          'z0_central_z', 'zn_stellarMass','zn_blackHoleMass',
          'zn_centralMvir','zn_mvir','zn_rvir',]] /= h

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


def binit(stats, labs, labels, N = 12, minimum = 1):
    """
    initialise bins and limits, calculate binned statistics
    
    Args:
        stats: completeness and purity statistics 
        labs: list of labelled regions
        labels: list of label strings
    """

    dgal = stats[:,0] + 1
    dgal_max = float(np.max(dgal))
    
    binLimits = np.linspace(0, dgal_max+(float(dgal_max)/N/2), N)

    lower_bin = binLimits[1] + (binLimits[0]-binLimits[1])/2. 
    upper_bin = binLimits[-1] + (binLimits[0]-binLimits[1])/2.

    bins = np.linspace(lower_bin, upper_bin, N-1)
    
    agg_total = np.zeros(len(bins))
    
    agg = {x: np.histogram(dgal[labs==x], binLimits)[0] for x in labels}  # save counts for each label
    agg_total = np.sum([v for k,v in agg.iteritems()],axis=0)   # find total in each bin
    
    # while np.sum(agg_total == 0) > 0:

    bins, binLimits, agg_total = merge_bins(bins, binLimits, agg_total, minimum=minimum)  # merge empty bins
        
    agg = {x: np.histogram(dgal[labs==x], binLimits)[0] for x in labels}  # save counts for each label
    agg_total = np.sum([v for k,v in agg.iteritems()],axis=0)   # find total in each bin
    
    fracs = {k: v.astype(float) / agg_total for k,v in agg.iteritems()}
    
    return bins, binLimits, agg, agg_total, fracs


def merge_bins(bins, binlimits, count, minimum=1):
    """
    merges empty bins with neighbours
    
    Args:
    - bins: length N
    - binlimits: lenght N+1
    - count: length N
    """
    
    # ids of offending bins
    empty = np.where(count < minimum)[0]
    
    while len(empty) > 0:
        
        idx = empty[-1]
        
        if idx == 0:
            print "reached bottom edge!"
            return 0

        # delete bin edge below
        binlimits = np.delete(binlimits, idx)

        # delete bin
        bins = np.delete(bins, idx)

        # create new bin
        bins[idx - 1] = binlimits[idx - 1] + (binlimits[idx] - binlimits[idx-1])/2

        # sum count values
        count[idx - 1] += count[idx]
        count = np.delete(count, idx)
        
        empty = np.where(count < minimum)[0]

        
    return bins, binlimits, count
    

def plotit(ax, stats, axb=None, clim = 0.5, plim = 0.5, N = 12, mlim=5e4, noplot=False, minimum=1):
    """
    
    Args:
        ax - axis object
        selection - selection object
        rid - id of radius selection in selection object
        zis - if od redshift selection in selection object
        axb - bottom axis object. If None, only plots 
    """
    
    colors = ['dimgrey','lightseagreen','lightcoral', 'y']
    
      
    mass = stats[:,3]
    
    dgal = stats[:,0] + 1

    # labels = ['proto_lomass','proto_himass','part_lomass','part_himass','pfield_lomass','pfield_himass','field']
    labs, labels = label(stats, clim, plim, mlim)
    
    bins, binLimits, agg, agg_total, fracs = binit(stats, labs, labels, N = N, minimum=minimum)

    if axb != None:  # probability density function
        
        phiMax = 0.
        
        mask = np.array([labs == lab for lab in ['proto_lomass','proto_himass','part_lomass','part_himass']]).any(axis=0)
        phiA = np.histogram(dgal[mask], binLimits, normed=True)[0]

        mask = (labs == 'field')
        phiB = np.histogram(dgal[mask], binLimits, normed=True)[0]
        
        phiMax = np.max(phiA)
        
        DB, BC = bhattacharyya(phiA*np.diff(binLimits), phiB*np.diff(binLimits))        
        
        mask = np.array([labs == lab for lab in ['proto_himass','part_himass']]).any(axis=0)
        phiC = np.histogram(dgal[mask], binLimits, normed=True)[0]
        
        DB_himass, BC_himass = bhattacharyya(phiC*np.diff(binLimits), phiB*np.diff(binLimits))

        if noplot:
            # print "DB(All), DB(High mass)"
            return round(DB, 3), round(DB_himass, 3)
            exit
        
        
        axb.step(bins, phiB, color=colors[0], linestyle='solid', linewidth=3)
        
        for lab, ls, c in zip(labels[:5], ['solid','dashed','solid','dashed'],[1,1,3,3]):
            mask = (labs==lab)
            if np.sum(mask) > 1:
                phi = np.histogram(dgal[mask], binLimits, normed=True)[0]
                axb.step(bins, phi, color=colors[c], linestyle=ls)
                phiMax = np.max([phiMax, np.max(phi)])
            
            
        axb.step(bins, phiA, color=colors[1], linestyle='solid', linewidth=3)
        axb.step(bins, phiA, color=colors[3], linestyle='dotted', linewidth=3)
        
        axb.set_ylim(0, phiMax + 0.1)
    
    
    width = binLimits[2] - binLimits[1]
    
    plt.rcParams['hatch.color'] = 'black'
    plt.rcParams['hatch.linewidth'] = 0.5
    
    # ax.bar(bins, agg['proto_lomass'] / agg_total, width=width, label='protocluster ($M<M_{lim}$)', alpha=0.8, color=colors[1])
    ax.bar(bins, fracs['proto_lomass'], width=[binLimits[i] - binLimits[i+1] for i in range(len(binLimits)-1)],
           label='protocluster ($M<M_{lim}$)', alpha=0.8, color=colors[1])
    
    # bar = ax.bar(bins, agg['proto_himass'] / agg_total, width=width, bottom=agg[labels[0]] / agg_total, label='protocluster ($M>M_{lim}$)', alpha=0.8, color=colors[1], hatch='///')
    ax.bar(bins, fracs['proto_himass'], width=[binLimits[i] - binLimits[i+1] for i in range(len(binLimits)-1)],
           bottom=fracs[labels[0]], label='protocluster ($M>M_{lim}$)', alpha=0.8, color=colors[1], hatch='///')
                 
    ax.bar(bins, fracs['part_lomass'], width=[binLimits[i] - binLimits[i+1] for i in range(len(binLimits)-1)],
           bottom=np.sum([fracs[key] for key in labels[0:2]],axis=0), label='part of a \n protocluster', alpha=0.8, color=colors[3]) 
    
    ax.bar(bins, fracs['part_himass'], width=[binLimits[i] - binLimits[i+1] for i in range(len(binLimits)-1)],
           bottom=np.sum([fracs[key] for key in labels[0:3]],axis=0), label='part of a \n protocluster', alpha=0.8, color=colors[3], hatch='///') 

    ax.bar(bins, fracs['pfield_lomass'], width=[binLimits[i] - binLimits[i+1] for i in range(len(binLimits)-1)],
           bottom= np.sum([fracs[key] for key in labels[0:4]],axis=0), label='protocluster \n + field', alpha=0.8, color=colors[2])
    
    ax.bar(bins, fracs['pfield_himass'], width=[binLimits[i] - binLimits[i+1] for i in range(len(binLimits)-1)],
           bottom= np.sum([fracs[key] for key in labels[0:5]],axis=0), label='protocluster \n + field', alpha=0.8, color=colors[2], hatch='///')

    ax.bar(bins, fracs['field'], width=[binLimits[i] - binLimits[i+1] for i in range(len(binLimits)-1)],
           bottom= np.sum([fracs[key] for key in labels[0:6]],axis=0), color=colors[0], 
           label='field', alpha=0.2)
    
    ax.set_ylim(0,1)
    ax.set_xlim(binLimits[0], binLimits[-1])
    
    if axb:
        axb.set_xlim(binLimits[0], binLimits[-1])
