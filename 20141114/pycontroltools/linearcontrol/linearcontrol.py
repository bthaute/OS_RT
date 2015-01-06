# -*- coding: utf-8 -*-
"""
.. module:: linearcontrol
    :synopsis: Functions concerning linear control algorithms.
    
.. moduleauthor:: Carsten Knoll
"""
import sympy as sp
import numpy as np

import itertools as it

import pycontroltools.auxfuncs.math.matrix
import pycontroltools.auxfuncs.math.miscmath
import pycontroltools.auxfuncs.math.numtools

def cont_mat(A,B):
    """
    Kallmanns controlability matrix
    """
    A = sp.Matrix(A)
    B = sp.Matrix(B)

    assert A.shape[0] == A.shape[1]
    assert A.shape[0] == B.shape[0]
    assert 1 == B.shape[1]

    n = A.shape[0]

    Q = sp.Matrix(B)
    for i in range(n-1):
        Q = Q.row_join(A**i * B)
        Q = Q.applyfunc(sp.expand)

    return Q
    
def is_left_coprime(Ap, Bp=None, eps = 1e-10):
    """
    Test ob Ap,Bp Linksteilerfrei sind
    keine Parameter zulässig

    """
    matrix = pycontroltools.auxfuncs.math.matrix
    miscmath = pycontroltools.auxfuncs.math.miscmath

# folgendes könnte die Berechnung vereinfachen
#    r1, r2 = Ap.shape
#
#    assert r1 <= r2
#
#    minors = all_k_minors(Ap, r1)
#
#    minors = list(set(minors)) # make entries unique


    #print "Achtung, BUG: Wenn ein minor konstant (auch 0) ist kommt ein Fehler"
    r1, r2 = Ap.shape
    if Bp == None:
        # interpret the last columns of Ap as Bp
        Bp = Ap[:, r1:]
        Ap = Ap[:, :r1]
        r1, r2 = Ap.shape


    assert r1 == r2
    r = r1
    r1, m =  Bp.shape
    assert r1 == r
    assert m <= r

    M = (Ap*1).row_join(Bp)

    symbs = list(matrix.matrix_atoms(M, sp.Symbol))
    assert len(symbs) == 1
    symb = symbs[0]

    combinations = it.combinations(range(r+m), r)

    minors = [matrix.col_minor(M, *cols) for cols in combinations]

    nonzero_const_minors = [m for m in minors if (m !=0) and (symb not in m)]

    if len(nonzero_const_minors) > 0:
        return True

    #zero_minors = [m for m in minors if m == 0]

    # polymionors (rows belong together)
    all_roots = [miscmath.roots(m) for m in minors if symb in m]

    # obviously all minors where zeros
    if len(all_roots) == 0:
        return False

    # in general the arrays in that list have differnt lenght
    # -> append a suitable number of roots at inf

    max_len = max([len(a) for a in all_roots])
    root_list = [np.r_[a, [np.inf]*(max_len-len(a))] for a in all_roots]

    all_roots = np.array(root_list)

    # now testing if some finite root is common to all minors
    for i in range(all_roots.shape[0]):
        test_roots = all_roots[i, :]

        other_roots = np.delete(all_roots, i, axis = 0)

        for tr in test_roots:
            if tr == np.inf:
                continue

            min_dist = np.min(np.abs(other_roots-tr), axis = 1)
            if np.all(min_dist < eps):
                # the test root was found in all other minors

                print "critical root:", tr
                return False

    return True
    
def linear_input_trafo(B, row_idcs):
    """
    serves to decouple inputs from each other
    """
    matrix = pycontroltools.auxfuncs.math.matrix
    numtools = pycontroltools.auxfuncs.math.numtools
    
    B = sp.Matrix(B)

    n, m = B.shape

    assert len(row_idcs) == m
    assert len(set(row_idcs)) == m

    Bnew = sp.zeros(m,m)
    for i, idx in enumerate(row_idcs):
        Bnew[i, :] = B[idx, :]

    P = matrix.symbMatrix(m,m)

    leqs = Bnew*P-sp.eye(m)
    res = sp.solve(leqs)

    return numtools.to_np(P.subs(res))