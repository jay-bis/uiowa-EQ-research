"""
INPUTS: dist - magnitude of closest points?? (double)
        p - pNorm (integer)
OUTPUTS: weights 

@author: Jack Biscupski
"""

def computeWeights(dist, p):
    reg = 1e-8
    return(p/(dist^(2-p) + reg))
