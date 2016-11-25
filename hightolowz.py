
import numpy as np
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
