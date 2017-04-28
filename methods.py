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


