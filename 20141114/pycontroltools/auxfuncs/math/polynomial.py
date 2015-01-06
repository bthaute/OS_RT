# -*- coding: utf-8 -*-
"""
.. module:: polynomial
    :synopsis: Functions concerning the construction of polynomials.
    
.. moduleauthor:: Carsten Knoll
"""
import sympy as sp
import numpy as np

import pycontroltools.auxfuncs.programming.miscprog
import miscmath

def condition_poly(var, *conditions):
    """
    # this function is intended to be a generalization of trans_poly

    returns a polynomial y(t) that fullfills given conditions

    every condition is a tuple of the following form:

    (t1, y1,  *derivs) # derivs contains cn derivatives

    every derivative (to the highest specified [in each condition]) must be given
    """
    assert len(conditions) > 0

    #assert t1 != t2


    # store the derivs
#    D1 = left[2:]
#    D2 = right[2:]




    # preparations
    cond_lengths = [len(c)-1 for c in conditions] # -1: first entry is t
    condNbr = sum(cond_lengths)
    cn = max(cond_lengths)

    coeffs = map(lambda i: sp.Symbol('a%d' %i), range(condNbr))
    #poly =  (map(lambda i, a: a*var**i, range(condNbr), coeffs))
    #1/0
    poly =  sum(map(lambda i, a: a*var**i, range(condNbr), coeffs))

    Dpoly_list = [poly]+map(lambda i: sp.diff(poly, var, i), range(1,cn+1))

    new_conds = []
    for c in conditions:
        t = c[0]
        for i,d in enumerate(c[1:]):
            new_conds.append((t,d,i))
            # d : derivative at point t (including 0th)
            # i : derivative counter

    # evaluate the conditions

    conds = []

    for t,d,i in new_conds:
        conds += [miscmath.equation(Dpoly_list[i].subs(var, t) , d)]



    sol = miscmath.lin_solve_eqns(conds, coeffs)

    sol_poly = poly.subs(sol)

    return sol_poly




def trans_poly(var, cn, left, right):
    """
    returns a polynomial y(t) that is cn times continous differentiable

    left and right are sequences of conditions for the boundaries

    left = (t1, y1,  *derivs) # derivs contains cn derivatives

    """
    assert len(left) == cn+2
    assert len(right) == cn+2

    t1, y1 = left[0:2]
    t2, y2 = right[0:2]

    assert t1 != t2

    for tmp in (y1, y2):
        assert not isinstance(tmp, (np.ndarray, np.matrix, sp.Symbol) )


    # store the derivs
    D1 = left[2:]
    D2 = right[2:]




    # preparations
    condNbr = 2 + 2*cn

    coeffs = map(lambda i: sp.Symbol('a%d' %i), range(condNbr))
    #poly =  (map(lambda i, a: a*var**i, range(condNbr), coeffs))
    #1/0
    poly =  sum(map(lambda i, a: a*var**i, range(condNbr), coeffs))

    Dpoly = map(lambda i: sp.diff(poly, var, i), range(1,cn+1))


    # create the conditions

    conds = []
    conds += [miscmath.equation(poly.subs(var, t1) , y1)]
    conds += [miscmath.equation(poly.subs(var, t2) , y2)]

    for i in range(cn):
        #

        conds += [miscmath.equation(Dpoly[i].subs(var, t1) , D1[i])]
        conds += [miscmath.equation(Dpoly[i].subs(var, t2) , D2[i])]


    sol = miscmath.lin_solve_eqns(conds, coeffs)

    sol_poly = poly.subs(sol)

    return sol_poly

def poly_degree(expr, var=None):
    """
    returns degree of monovariable polynomial
    """
    miscprog = pycontroltools.auxfuncs.programming.miscprog
    
    var = miscprog.get_expr_var(expr, var)
    if var == None:
        return sp.sympify(0)

    P = sp.Poly(expr, var, domain = "EX")
    return P.degree()
    
def poly_coeffs(expr, var=None):
    """
    returns all (monovariate)-poly-coeffs (including 0s) as a list
    first element is highest coeff.
    """
    miscprog = pycontroltools.auxfuncs.programming.miscprog
    
    var = miscprog.get_expr_var(expr, var)
    if var == None:
        return [expr]

    P = sp.Poly(expr, var, domain="EX")

    pdict = P.as_dict()

    d = P.degree()

    return [pdict.get((i,), 0) for i in reversed(xrange(d+1))]
    
def coeffs(expr, var = None):
    # TODO: besser über as_dict
    # TODO: überflüssig wegen poly_coeffs?
    """if var == None, assumes that there is only one variable in expr"""
    expr = sp.sympify(expr)
    if var == None:
        vars = filter(lambda a:a.is_Symbol, list(expr.atoms()))
        if len(vars) == 0:
            return [expr] # its a constant
        assert len(vars) == 1
        var=vars[0]
        dom = 'RR'
    else:
        dom = 'EX'
    return sp.Poly(expr, var, domain =dom).all_coeffs()
    
def zeros_to_coeffs(*z_list, **kwargs):
    """
    calculates the coeffs corresponding to a poly with provided zeros
    """

    s = sp.Symbol("s")
    p = sp.Mul(*[s-s0 for s0 in z_list])

    real_coeffs = kwargs.get("real_coeffs", True)
    c = np.array(coeffs(p, s), dtype=np.float)

    if real_coeffs:
        c = np.real(c)
    return c

def poly_scalar_field(xx, symbgen, order, poly=False):
    """
    returns a multivariate poly with specified oders
    and symbolic coeffs
    returns also a list of the coefficients
    """
    miscprog = pycontroltools.auxfuncs.programming.miscprog

    if isinstance(order, int):
        orders = [order]
    elif isinstance(order, (list, tuple, sp.Matrix)):
        orders = list(order)

    res = 0
    coeff_list = []
    for i in orders:
        if i == 0:
            c = symbgen.next()
            res += c
            coeff_list.append(c)
            continue

        terms = miscprog.get_diffterms(xx, i)

        for tup in terms:
            c = symbgen.next()
            res += c*sp.Mul(*tup)
            coeff_list.append(c)

    if poly:
        res = sp.Poly(res, *xx, domain='EX')
    return res, coeff_list

def get_order_coeff_from_expr(expr, symb, order):
    """
    example:
        3*s**2 -4*s + 5, s, 3 -> 0
        3*s**2 -4*s + 5, s, 2 -> 3
        3*s**2 -4*s + 5, s, 1 -> -4
        3*s**2 -4*s + 5, s, 9 -> 0
    """
    p = sp.Poly(expr, symb, domain = "EX")
    default = 0
    return p.as_dict().get( (order,), default )

def element_deg_factory(symb):
    """
    returns a function for getting the polynomial degree of an expr. w.r.t.
    a certain symbol
    """
    def element_deg(expr):
        return poly_degree(expr, symb)

    return element_deg

