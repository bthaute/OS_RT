# -*- coding: utf-8 -*-
"""
.. module:: example_lietools
    :synopsis: Involutivity example for testing lietools.
    
.. moduleauthor:: Carsten Knoll
"""
import lietools as lt

import sympy as sp
import numpy as np

import itertools as it

def involutive(f_list, *args):
    r"""
        Checks if a distribution :math:`\Delta = \mathrm{span}(f_1,\dots,f_r)`
        is involutive.
    """
    assert len(args) > 0

    if hasattr(args[0], '__len__'):
        args = args[0]  
        
    r = len(f_list)
    aa = sp.symbols('a1:%i' % (r+1), cls=sp.Dummy)
    
    # for every combination (f1,f2) in f_list
    for f1, f2 in list(it.combinations(f_list, 2)):
        # calculate lie bracket
        bracket = lt.lie_bracket(f1, f2, args)
        
        # building equation systems
        eq = bracket
        for a_symb, elem in zip(aa, f_list):
            eq -= a_symb*elem
        
        # try to solve the equation system
        # if there is no solution,
        # the vector fields are not linearly dependent
        # and the distributin is not involutive
        if len(sp.solve(eq,aa)) == 0:
            return False
        else:
            return True
            

# INVOLUTIVITY EXAMPLE 1
# (see Applied Nonlinear Control by Slotine&Li (Chapter 6.2, Example 6.9))

x1,x2,x3 = sp.symbols('x1 x2 x3')

# involutive distribution defined by the vector fields
f1 = sp.Matrix([4*x3, -1, 0])
f2 = sp.Matrix([-x1, x3**2-3*x2, 2*x3])

# checking involutivity (simple)
# lie bracket [f1,f2]
ff12 = lt.lie_bracket(f1,f2,[x1,x2,x3])

M12 = sp.Matrix(np.hstack([f1,f2,ff12]))
assert sp.det(M12) == 0

# checking involutivity (with function)
assert involutive([f1,f2], x1,x2,x3) == True

## INVOLUTIVITY EXAMPLE 2

# non-involutive distribution defined by the vector fields
f3 = sp.Matrix([4*x3, -1, 0])
f4 = sp.Matrix([-x2, x3**2-3*x2, 2*x1])

assert involutive([f3,f4], x1,x2,x3) == False

# INVOLUTIVITY EXAMPLE 3
# see http://www.roebenack.de/content/example-involutiveintegrable-distribution

# involutive distribution defined by the vector fields
g1 = sp.Matrix([2*x3, -1, 0])
g2 = sp.Matrix([-x1,-2*x2, x3])

assert involutive([g1,g2], x1,x2,x3) == True


# LJAPUNOV'S DIRECT METHOD EXAMPLE
# (see Applied Nonlinear Control by Slotine&Li (Chapter 3.4, Example 3.7))

a,b = sp.symbols('a b')

# vector field
fx = sp.Matrix([x2,-a*sp.sin(x1) - b*x2])
# positive definite scalar function
Vx = 0.5*x2**2 + a*(1 - sp.cos(x1))

# negative semi definite Lie derivative
Vx_dot = lt.lie_deriv(Vx, fx, [x1,x2])
# simplified
Vx_dot = sp.simplify(Vx_dot)
