# -*- coding: utf-8 -*-
"""
.. module:: diffoperators
    :synopsis: Functions concerning differential operators.
    
.. moduleauthor:: Carsten Knoll
"""

import sympy as sp
import numpy as np

def hoderiv(f, x, N=2):
    """
    computes a H igher O rder derivative of the vectorfield f

    Result is a tensor of type (N,0)

    or a  n x L x ... x L (N times) hyper Matrix

    (represented a (N+1)-dimensional numpy array

    """

    import itertools as it

    assert f.shape[1] == 1

    n = f.shape[0]
    L = len(list(x))

    res = np.zeros([n]+[L]*N)

    res= res * sp.Symbol('dummy')

    idx_list = [0]*(N)
    i = 0
    k = 0


    list_of_idcs = list(it.product(*[range(L)]*N))

    # example: [(0, 0), (0, 1), (1, 0), (1, 1)]

    for fi in f:
        #print fi, i
        for idcs in list_of_idcs:
            pos = tuple([i]+list(idcs))

            tmp = fi
            #print pos
            for j in idcs:
                #print x[j]
                tmp = tmp.diff(x[j])

            #print "---\n"
#            res.itemset(pos, k)
            res.itemset(pos, tmp)
            k+=1
            #print "  .. ", i, k

        i+=1

    return res
    
def div(vf, x):
    """divergence of a vector field"""
    vf = list(vf)
    x = list(x)
    assert len(vf) == len(x)

    return sum([c.diff(xi) for c,xi in zip(vf, x)])
    
def gradient(scalar_field, xx):
    """
    # returns a row vector (coverctorfiel)!
    """
    return sp.Matrix([scalar_field]).jacobian(xx)