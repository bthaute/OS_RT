# -*- coding: utf-8 -*-
"""
.. module:: taylor
    :synopsis: Functions concerning the construction of Taylor polynomials.
    
.. moduleauthor:: Carsten Knoll
"""

import sympy as sp
import itertools as it

t = sp.var('t')

def multi_taylor(expr, args, x0 = None, order=1):
    """
    compute a multivariate taylor polynomial of a scalar function

    default: linearization about 0 (all args)
    """

    if x0 == None:
        x0 = [0 for a in args]
    x0 = list(x0) # to handle matrices
    assert len(args) == len(x0)

    x0list = zip(args, x0)

    res = expr.subs(x0list)

    arg_idx_list = range( len(args) )

    for o in xrange(1,order+1):

        diff_list = it.product( *([arg_idx_list]*o) )

        for idx_tup in diff_list:

            arg_tup = [args[k] for k in idx_tup]

            prod = sp.Mul( *[args[k]-x0[k] for k in idx_tup] )

            tmp = expr.diff(*arg_tup)/sp.factorial(o)

            res+= tmp.subs(x0list) * prod
    return res


def multi_taylor_matrix(M, args, x0=None, order=1):
    """
    applies multi_taylor to each element
    """

    def func(m):
        return multi_taylor(m, args, x0, order)

    return M.applyfunc(func)
    
def series(expr, var, order):
    """
    taylor expansion at zero (without O(.) )
    """
    if isinstance(expr, sp.Matrix):
        return type(expr)(map(lambda x: series(x, var, order), expr))

    # expr is scalar
    expr = expr.series(var, 0, order).removeO()
    p = sp.Poly(expr, var, domain='EX')
    s = 0

    #!! limit the order (due to a sympy bug this is not already done)
    for i in range(order+1):
        s+= p.nth(i) * t**i

    return s