# -*- coding: utf-8 -*-
"""
.. module:: miscmath
    :synopsis: Miscellaneous mathematical helpfunctions.
    
.. moduleauthor:: Carsten Knoll
"""
#__all__ = ['equation', 'make_eqns', 'fractionfromfloat', 'sp_fff']

import sympy as sp
import numpy as np
import fractions as fr

import pycontroltools.auxfuncs.programming.miscprog
import polynomial

from collections import Counter

Fr = fr.Fraction

class equation(object):
    """#chris:  Klasse equation erstellt Gleichungs-Objekte mittles sympify
                mit Attributen für Lefthandside (lhs) und Righthandside (rhs)
                der Gleichung
    """
    def __init__(self, lhs, rhs = 0):
        self.lhs_ = sp.sympify(lhs)
        self.rhs_ = sp.sympify(rhs)

    def lhs(self):
        return self.lhs_

    def rhs(self):
        return self.rhs_

    def __repr__(self):
        return "%s == %s" % (self.lhs_, self.rhs_)

    def subs(self, *args, **kwargs):
        lhs_  = self.lhs_.subs(*args, **kwargs)
        rhs_  = self.rhs_.subs(*args, **kwargs)
        return type(self)(lhs_, rhs_)

def make_eqns(v1, v2 = None):
    """
    #chris: mehrere lhs,rhs übergeben und daraus Gleichungen erstellen
    """
    if v2 == None:
        v2 = [0]*len(list(v1))
    return [equation(v1[i], v2[i]) for i in range(len(list(v1)))]
    

#chris: denominator = Nenner
#       stellt Gleitkommazahl mit maximalem Nenner als
#       fractions.Fraction-Objekt dar (1.33 --> Fraction(133, 100))
def fractionfromfloat(x_, maxden = 1000):
  """
  fraction from float
  args:
   x
   maxdenominator (default = 1000)
  """

  x = float(x_)
  assert x == x_ # fails intentionally for numpy.complex
  return Fr.from_float(x).limit_denominator(maxden)

def sp_fff(x, maxden):
    """ sympy_fraction from float
    #chris: nimmt anscheinend Objekte vom Typ fractions.Fraction
            (Fraction(133, 10)) und stellt sie als Bruch dar (133/10)
    """
    return sp.Rational(fractionfromfloat(x, maxden))
    
def numer_denom(expr):
    """
    """
    num, denom = expr.as_numer_denom() # resolves problems with multifractions
    return num/denom
    
def uv(n, i):
    """
    unit vectors (columns)
    """
    uv = sp.Matrix([0]*n)
    uv[i-1] = sp.sympify(1)
    return uv

def symbs_to_func(expr, symbs, arg):
    """
    in expr replace x by x(arg)
    where x is any element of symbs
    """
    #TODO: assert all([isinstance(s, sp.Symbol) for s in symbs])
    funcs = [sp.Function(s.name)(arg) for s in symbs]

    return expr.subs(zip(symbs, funcs))
    
def jac(expr, *args):
    """
    Calculates the Jacobian matrix (derivative of a vectorial function)
    using the :func:`jacobian` function from the module
    :mod:`sympy.matrices.matrices.MatrixBase`.
    
    **Advantage:** direct derivation of functions
    
    Jacobian matrix:
    
        .. math::
            
            J_f(a) :=  \\left(\\begin{matrix} 
                        \\frac{\partial f_1}{\partial x_1}(a) &
                        \\frac{\partial f_1}{\partial x_2}(a) &
                        \\ldots &
                        \\frac{\partial f_1}{\partial x_n}(a)\\\\
                        \\vdots & \\vdots & \\ddots & \\vdots \\\\
                        \\frac{\partial f_m}{\partial x_1}(a) &
                        \\frac{\partial f_m}{\partial x_2}(a) & \ldots &
                        \\frac{\partial f_m}{\partial x_n}(a)
                        \\end{matrix}\\right)
                    
    **Parameters**
    
    * expr : expression to derive
        function / row matrix/ column matrix
    * args : coordinates
        separate or as list-like object
        
    **Return**
    
    * returns : Jacobi matrix
    * type : sympy.Matrix
    
    **Examples**
    
        >>> import sympy
        >>> x1,x2,x3 = sympy.symbols('x1 x2 x3')
        
        >>> jac(x1**2+2*x2+x3, x1, x2, x3)
        Matrix([[2*x1, 2, 1]])
        
    .. seealso:: :func:`sympy.jacobian`
    """
    
    if not hasattr(expr, '__len__'):
        expr = [expr]
    return sp.Matrix(expr).jacobian(args)
    
#chris: Erklärung?
def get_coeff_row(eq, vars):
    """
    takes one equation object and returns the corresponding row of
    the system matrix
    """
    if not isinstance(eq, equation):
        # assume its the lhs     and rhs = 0
        eq = equation(eq,0)

    if isinstance(vars, sp.Matrix):
        vars = list(vars)

    get_coeff = lambda var: sp.diff(eq.lhs(), var)
    coeffs =  map(get_coeff, vars)
    rest = eq.lhs() - sum([coeffs[i]*vars[i] for i in range( len(vars) )])
    coeff_row = map(get_coeff, vars) + [eq.rhs() - rest]
    return coeff_row

def lin_solve_all(eqns):
    """
    takes a list of equations and tries to solve wrt. to all
    ocurring symbols
    """
    miscprog = pycontroltools.auxfuncs.programming.miscprog
    eqns = sp.Matrix(eqns)

    Vars = list(miscprog.atoms(eqns, sp.Symbol))

    return lin_solve_eqns(eqns, Vars)


def lin_solve_eqns(eqns, vars):
    """
    takes a list of equation objects
    creates a system matrix of and calls sp.solve
    """
    n = len(eqns)

    vars = list(vars) # if its a vector
    m = len(vars)

    rows = [get_coeff_row(eq, vars) for eq in eqns]

    sysmatrix = sp.Matrix(rows)

    sol = sp.solve_linear_system(sysmatrix, *vars)

    return sol

def lin_solve_eqns_jac(eqns, vars):
    """
    takes a list of equation objects
    creates a system matrix of and calls sp.solve

    # new version !!
    # should replace lin_solve_eqns

    # assumes that eqns is a list of expressions where rhs = 0
    """
    miscprog = pycontroltools.auxfuncs.programming.miscprog    
    
    eqm = sp.Matrix(eqns)

    Jac = eqm.jacobian(vars)
    rhs = -eqm.subs(miscprog.zip0(vars))

    sysmatrix = Jac.row_join(rhs)

    sol = sp.solve_linear_system(sysmatrix, *vars)

    return sol
    
def extract_independent_eqns(M):
    """
    handles only homogeneous eqns

    M Matrix

    returns two lists: indices_of_rows, indices_of_cols

    """

    #S = M # save the symbolical matrix for later use
    M = (np.array(M)!= 0)*1 # Matrix of ones and zeros
    n, m = M.shape

    list_of_occtuples = []

    for i in range(n):
        tmplist = []
        for j in range(m):
            if M[i,j] == 1:
                tmplist.append(j)

        list_of_occtuples.append(tuple(tmplist))

    print list_of_occtuples

    #list_of_occtuples_save = list_of_occtuples[:]

    s0 = set(list_of_occtuples[0])

    list_of_rows=[0]

    while True:
        i=-1
#        print "s0="
#        print s0
#        print
        end = False
        for ot in list_of_occtuples:
            i+=1
            if i in list_of_rows:
                continue

#            print i
            if s0.intersection(ot) != set():
                s0 = s0.union(ot)
                #print " >",i
                list_of_rows.append(i)
                break # restart for loop


#        if end == True:
        if i == len(list_of_occtuples)-1:
            #print "durch"
            break

    s0 = list(s0)
    return sorted(list_of_rows), sorted(s0)
    
def rat_if_close(x, tol=1e-10):
    """
    """
    s = sp.sympify(x)

    maxden = int(tol**-1 / 10.0)
    f  = fractionfromfloat(x, maxden)
    r = sp.Rational(f.numerator, f.denominator)
    if abs(r-x) < tol:
        return r
    else:
        return x

def rationalize_expression(expr, tol=1e-10):
    """
    substitutes real numbers occuring in expr which are closer than tol to a
    rational with a sufficiently small denominator with these rationals

    usefull special case 1.2346294e-15 -> 0

    """
    a = list(expr.atoms(sp.Number))
    b = [rat_if_close(aa, tol) for aa in a]

    return expr.subs(zip(a,b))
    
def roots(expr):
    """
    """
    import scipy as sc
    return sc.roots(polynomial.coeffs(expr))
    
def real_roots(expr):
    """
    """
    import scipy as sc
    r = sc.roots(polynomial.coeffs(expr))
    return np.real( r[np.imag(r)==0] )
    
def fac(i):
    # see also sp.factorial
    if i == 0:
        return 1
    return i * fac(i-1)
    
def trigsimp2(expr):
    """
    sin**2 + cos**2 = 1 in big expressions
    """
    miscprog = pycontroltools.auxfuncs.programming.miscprog
    expr = expr.expand()

    trigterms_sin = list(expr.atoms(sp.sin))
    trigterms_cos = list(expr.atoms(sp.cos))


    trigterm_args = []

    # gucken, ob cos auch vorkommt
    for tts in trigterms_sin:
        arg = tts.args[0]
        if sp.cos(arg) in trigterms_cos:
            trigterm_args.append(arg)



    for s in trigterm_args:
        poly = miscprog.trig_term_poly(expr, s)

        dd = poly.as_dict()


        uncommon_coeff = (dd[(2,0)] - dd[(0,2)]).expand()



        if uncommon_coeff == 0:
            expr += dd[(2,0)] - dd[(2,0)]*sp.sin(s)**2 - dd[(2,0)]*sp.cos(s)**2


        print dd[(2,0)]

    return expr.expand()
    
def multi_series(expr, xx, order, poly=False):
    """
    Reihenentwicklung (um 0) eines Ausdrucks in mehreren Variablen
    """
    miscprog = pycontroltools.auxfuncs.programming.miscprog
    
    xx0 = zip(xx, [0]*len(xx)) # Entwicklungsstelle
    res = 0
    for i in range(order+1):
        if i == 0:
            res += expr.subs(xx0)
            continue
        terms = miscprog.get_diffterms(xx, i)
        for tup in terms:
            cnt = Counter(tup) # returns a dict
            fac_list = [sp.factorial(n) for n in cnt.values()]
            #fac = 1/sp.Mul(*fac_list)
            res += expr.diff(*tup).subs(xx0)*sp.Mul(*tup) / sp.Mul(*fac_list)

    if poly:
        res = sp.Poly(res, *xx, domain="EX")
    return res