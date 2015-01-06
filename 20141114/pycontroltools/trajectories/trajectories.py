# -*- coding: utf-8 -*-
"""
.. module:: trajectories
    :synopsis: Functions concerning the construction of system trajectories.
    
.. moduleauthor:: Carsten Knoll
"""

import sympy as sp

# avoid name clashes with sage
piece_wise = sp.functions.elementary.piecewise.Piecewise


def make_pw(var, transpoints, fncs):
    """
    """
    transpoints = list(transpoints)
    upper_borders =  list(zip(*transpoints)[0])

    var = sp.sympify(var)

    inf = sp.oo

#    if len(upper_borders) == len(fncs)-1:
#        upper_borders += [inf]
#

    assert len(upper_borders) == len(fncs) -1
    upper_borders += [inf]

    #lower_borders = [-inf] + transpoints
    #fncs+=[fncs[-1]] # use the last fnc beyond the last transpoint

    # generate a list of tuples
    pieces = [(fnc, var < ub) for ub, fnc in zip(upper_borders, fncs)]
    #IPS()
    return piece_wise(*pieces)


def integrate_pw(fnc, var, transpoints):
    """
    due to a bug in sympy we must correct the offset in the integral
    to make the result continious
    """

    F=sp.integrate(fnc, var)

    fncs, conds = zip(*F.args)

    transpoints = list(zip(*transpoints)[0])

    oldfnc = fncs[0]
    new_fncs = [oldfnc]
    for f, tp  in zip(fncs[1:], transpoints):
        fnew = f + oldfnc.subs(var, tp) - f.subs(var, tp)
        new_fncs.append(fnew)
        oldfnc = fnew

    pieces = zip(new_fncs, conds)

    return piece_wise(*pieces)