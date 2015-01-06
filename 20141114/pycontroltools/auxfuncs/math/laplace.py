# -*- coding: utf-8 -*-
"""
.. module:: laplace
    :synopsis: Functions concerning LaPlace.
    
.. moduleauthor:: Carsten Knoll
"""

import sympy as sp

def do_laplace_deriv(laplace_expr, s, t):
    """
    Example:
    laplace_expr = s*(t**3+7*t**2-2*t+4)
    returns: 3*t**2  +14*t - 2
    """

    if isinstance(laplace_expr, sp.Matrix):
        return laplace_expr.applyfunc(lambda x: do_laplace_deriv(x, s,t))

    exp = laplace_expr.expand()

    #assert isinstance(exp, sp.Add)

    P = sp.Poly(exp, s, domain = "EX")
    items = P.as_dict().items()

    res = 0
    for key, coeff in items:
        exponent = key[0] # exponent wrt s

        res += coeff.diff(t, exponent)

    return res