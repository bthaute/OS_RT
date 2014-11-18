# -*- coding: utf-8 -*-


"""

useful functions on basis of sympy

"""

#from ipHelp import IPS, ST

import sympy as sp
import numpy as np

from collections import Counter


import random

import itertools as it


# convenience
np.set_printoptions(8, linewidth=300)

piece_wise = sp.functions.elementary.piecewise.Piecewise # avoid name clashes with sage

t = sp.var('t')

zf = sp.numbers.Zero()

class equation(object):

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


class Container(object):

    def __init__(self, **kwargs):
        assert len( set(dir(self)).intersection(kwargs.keys()) ) == 0
        self.__dict__.update(kwargs)


def make_eqns(v1, v2 = None):
    if v2 == None:
        v2 = [0]*len(list(v1))
    return [equation(v1[i], v2[i]) for i in range(len(list(v1)))]

#def pw_curve(var, cn, left, right):
#    """
#    returns a piecewise polynomial y(t) that is cn times continous differentiable
#
#    left and right are sequences of conditions for the boundaries
#
#    left = (t1, y1,  derivs) # derivs contains cn derivatives
#
#    """
#    pass



# Der Inhalt des quickfrac moduls
"""
quick access to the most useful functionality of the fractions module
"""

import fractions as fr

Fr = fr.Fraction

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
    """ sympy_fraction from float"""
    return sp.Rational(fractionfromfloat(x, maxden))



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
        conds += [equation(Dpoly_list[i].subs(var, t) , d)]



    sol = lin_solve_eqns(conds, coeffs)

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
    conds += [equation(poly.subs(var, t1) , y1)]
    conds += [equation(poly.subs(var, t2) , y2)]

    for i in range(cn):
        #

        conds += [equation(Dpoly[i].subs(var, t1) , D1[i])]
        conds += [equation(Dpoly[i].subs(var, t2) , D2[i])]


    sol = lin_solve_eqns(conds, coeffs)

    sol_poly = poly.subs(sol)

    return sol_poly


def make_pw(var, transpoints, fncs):
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






def lie_deriv(sf, vf, x, n = 1):
    """
    lie_deriv of a scalar field along a vector field
    """

    if isinstance(x, sp.Matrix):
        assert x.shape[1] == 1
        x = list(x)

    assert int(n) == n and n >= 0
    if n == 0:
        return sf

    res = jac(sf, x)*vf
    assert res.shape == (1,1)
    res = res[0]

    if n > 1:
        return lie_deriv(res, vf, x, n-1)
    else:
        return res

def lie_bracket(f, g, *args, **kwargs):
    """
    f, g should be vectors (or lists)

    call possibillities:

    lie_bracket(f, g, x1, x2, x3)
    lie_bracket(f, g, [x1, x2, x3])
    lie_bracket(f, g, sp.Matrix([x1, x2, x3]) )



    optional keyword arg n ... order
    """

    assert len(args) > 0


    if isinstance(args[0], sp.Matrix):
        assert args[0].shape[1] == 1
        args = list(args[0])

    if hasattr(args[0], '__len__'):
        args = args[0]
    n = kwargs.get('n', 1) # wenn n nicht gegeben, dann n=1

    if n == 0:
        return g

    assert n > 0 #and isinstance(n, int)
    assert len(args) == len(list(f))

    # Umwandeld in sympy-Matrizen
    f = sp.Matrix(f)
    g = sp.Matrix(g)

    jf = f.jacobian(args)
    jg = g.jacobian(args)


    res = jg * f - jf * g

    if n > 1:
        res = lie_bracket(f, res, *args, n=n-1)

    return res


def lie_deriv_covf(w, f, args, **kwargs):
    """
    Lie derivative of covector fields along vector fields

    w, f should be 1 x n and n x 1 Matrices


    (includes the option to omit the transposition of Dw
    -> transpose_jac = False)
    """

    k,l = w.shape
    m, n = f.shape
    assert  k==1 and n==1
    assert l==m

    if isinstance(args[0], sp.Matrix):
        assert args[0].shape[1] == 1
        args = list(args[0])

    if hasattr(args[0], '__len__'):
        args = args[0]

    assert len(args) == len(list(f))

    n = kwargs.get('n', 1) # wenn n nicht gegeben, dann n=1

    if n == 0:
        return w

    assert n > 0 #and isinstance(n, int)


    # caution: in sympy jacobians of row and col vectors are equal
    # -> transpose is needless (but makes the formula consistent with books)
    jwT = w.T.jacobian(args)

    jf = f.jacobian(args)



    if kwargs.get("transpose_jac", True) == False:
        # stricly this is not a lie derivative
        # but nevertheless sometimes needed
        res = w*jf + f.T * jwT
    else:

        # This is the default case :
        res = w*jf + f.T * jwT.T

    if n > 1:
        res = lie_deriv_covf(res, f, args, n = n-1)

    return res


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


def numer_denom(expr):
    num, denom = expr.as_numer_denom() # resolves problems with multifractions
    return num/denom


def expand(arg):
    """
    sp.expand currently has no matrix support
    """
    if isinstance(arg, sp.Matrix):
        return arg.applyfunc(sp.expand)
    else:
        return sp.expand(arg)



def simplify(arg):
    """
    sp.simplify currently has no matrix support
    """
    if isinstance(arg, sp.Matrix):
        return arg.applyfunc(sp.simplify)
    else:
        return sp.simplify(arg)


def trigsimp(arg, **kwargs):
    """
    sp.trigsimp currently has no matrix support
    """
    if isinstance(arg, (sp.Matrix, sp.ImmutableMatrix)):
        return arg.applyfunc(lambda x: sp.trigsimp(x, **kwargs))
    else:
        return sp.trigsimp(arg, **kwargs)


def ratsimp(arg):
    """
    sp.ratsimp currently has no matrix support
    """
    if isinstance(arg, sp.MatrixBase):
        return arg.applyfunc(sp.ratsimp)
    else:
        return sp.ratsimp(arg)



def uv(n, i):
    """
    unit vectors (columns)
    """
    uv = sp.Matrix([0]*n)
    uv[i-1] = sp.sympify(1)
    return uv

# Todo Unittest (sp.Symbol vs. sp.cls)
def symbMatrix(n, m, s='a', symmetric = 0, cls = sp.Symbol):
    """
    create a matrix with Symbols as each element
    cls might also be sp.Dummy
    """
    A = sp.Matrix(n,m, lambda i,j:cls( s+'%i%i'%(i+1,j+1)) )
    if symmetric == 1:
        subs_list = symmetryDict(A)
        A = A.subs(subs_list)
    return A


def symbs_to_func(expr, symbs, arg):
    """
    in expr replace x by x(arg)
    where x is any element of symbs
    """
    assert all([isinstance(s, sp.Symbol) for s in symbs])
    funcs = [sp.Function(s.name)(arg) for s in symbs]

    # conveniece: type casting
    if isinstance(expr, (tuple, list)):
        expr = sp.Matrix(expr)

    return expr.subs(zip(symbs, funcs))

def getOccupation(M):
    """
    maps (m_ij != 0) to every element
    """
    M = sp.sympify(M)
    n, m = M.shape
    tmp = sp.Matrix(n, m, lambda i,j: 1 if not M[i,j]==0 else 0)
    return tmp

def symmetryDict(M):
    """
    erstellt ein dict, was aus einer beliebigen Matrix M
    mittels M.subs(..) eine symmetrische Matrix macht
    """
    n, m = M.shape
    res = {}
    for i in range(1,n):
        for j in range(i):
            # lower triangle
            res[M[i,j]] = M[j,i]

    return res


# TODO: unit test
def make_global(varList, up_count = 0):
    """
    injects the symbolic variables of a collection to the global namespace
    usefull for interactive sessions

    :up_count: is the number of frames to go back;
    up_count = 0 means up to the upper_most frame
    """

    if not isinstance(varList, (list, tuple)):
        if isinstance(varList, sp.Matrix):
            varList = np.array(varList).flatten()
        else:
            raise TypeError, 'Unexpected type for varList'

    import inspect

    # get the topmost frame
    frame = inspect.currentframe()
    i = up_count
    while True:
        if frame.f_back == None:
            break
        frame = frame.f_back
        i -= 1
        if i == 0:
            break



    # this is strongly inspired by sympy.var
    try:
        for v in varList:
            if getattr(v, 'is_Function', False):
                v = v.func
            if hasattr(v, 'name'):
                frame.f_globals[v.name] = v
            elif hasattr(v, '__name__'):
                frame.f_globals[v.__name__] = v
            else:
                raise ValueError, 'Object %s has no name' % str(v)
    finally:
        # we should explicitly break cyclic dependencies as stated in inspect
        # doc
        del frame


makeGlobal = make_global

def prev(expr, **kwargs):
    """
    sympy preview abbreviation
    """
    KWargs = {'output':'pdf'}
    KWargs.update(kwargs)
    sp.preview(expr, **KWargs)


def mdiff(M, var):
    """
    returns the elementwise derivative of a matrix M w.r.t. var
    """
    return M.applyfunc(lambda elt: sp.diff(elt, var))

# TODO: seems to conflict with zip0
def tup0(xx):
    """
    helper function for substituting.
    takes (x1, x2, x3, ...)
    returns [(x1, 0), (x2, 0), ...]
    """

    return zip(xx, [0]*len(xx))

def jac(expr, *args):
    if not hasattr(expr, '__len__'):
        expr = [expr]
    return sp.Matrix(expr).jacobian(args)


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

def get_rows(A):
    """
    returns a list of n x 1 vectors
    """
    A = sp.Matrix(A)
    n, m = A.shape

    return [A[:,i] for i in range(m)]



# ---------------- July 2010 -----------

def elementwise_mul(M1, M2):
    """
    performs elment wise multiplication of matrices
    """
    assert M1.shape == M2.shape
    return sp.Matrix(np.array(M1)*np.array(M2))


def poly_occ_matrix(expr, arg1, arg2, n = 2):
    p = sp.Poly(expr, arg1, arg2)
    M = sp.matrices.zeros(n+1)
    star = sp.Symbol('#')
    space = sp.Symbol('.')
    for i in range(n+1):
        for j in range(n+1):
            if p.coeff(i,j) != 0:
                M[i,j] = star
            else:
                pass
                M[i,j] = space
    return M



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
    eqns = sp.Matrix(eqns)

    Vars = list(atoms(eqns, sp.Symbol))

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
    eqm = sp.Matrix(eqns)

    Jac = eqm.jacobian(vars)
    rhs = -eqm.subs(zip0(vars))

    sysmatrix = Jac.row_join(rhs)

    sol = sp.solve_linear_system(sysmatrix, *vars)

    return sol

#def lin_solve_eqns(eqns, vars):
#    """
#    takes a list of equation objects
#    creates a system matrix of and calls sp.solve
#    """
#    n = len(eqns)
#
#    vars = list(vars) # if its a vector
#    m = len(vars)
#
#
#
#    #sysmatrix = sp.Matrix(0,0).zeros((n,m+1)) # last col is col of rhs
#
#
#
#    def get_coeff_row(eq):
#        """
#        takes one equation object and returns the corresponding row of
#        the system matrix
#        """
#        get_coeff = lambda var: sp.diff(eq.lhs(), var)
#        coeffs =  map(get_coeff, vars)
#        rest = eq.lhs() - sum([coeffs[i]*vars[i] for i in range( len(vars) )])
#        coeff_row = map(get_coeff, vars) + [eq.rhs() - rest]
#        return coeff_row
#
#    rows = map(get_coeff_row, eqns)
#
#    sysmatrix = sp.Matrix(rows)
#
#
#    sol = sp.solve_linear_system(sysmatrix, *vars)
#
#    return sol


# !! nicht von allgemeinem Interesse

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


def cancel_rows_cols(M, rows, cols):
    """
    cancel rows and cols form a matrix

    rows ... rows to be canceled
    cols ... cols to be canceled
    """

    idxs_col = range(M.shape[1])
    idxs_row = range(M.shape[0])

    rows.sort(reverse=True)
    cols.sort(reverse=True)

    for c in cols:
        idxs_col.pop(c)

    for r in rows:
        idxs_row.pop(r)

#    all_coeffs = sp.Matrix(np.array(all_coeffs)[idxs_col])

    tmp = np.array(M)[idxs_row, :]
    M = sp.Matrix(tmp[:, idxs_col])

    return M


def norm2(v):
    """
    assumes that v is a n x 1 matrix (column vector)
    returns (v.T*v)[0,0] i.e. the squared 2-norm
    """
    assert isinstance(v, (sp.Matrix,))
    return (v.T*v)[0,0]


# TODO: Doctest
def concat_cols(*args):
    """
    takes some col vectors and aggregetes them to a matrix
    """

    col_list = []

    for a in args:
        if a.shape[1] == 1:
            col_list.append( list(a) )
            continue
        for i in xrange(a.shape[1]):
            col_list.append( list(a[:,i]) )
    m = sp.Matrix(col_list).T

    return m

# other name:
col_stack = concat_cols


def col_split(A, *indices):
    """
    returns a list of columns corresponding to the passed indices
    """
    if not indices:
        indices = range(A.shape[1])
    res = [ A[:, i] for i in indices ]
    return res


def crow_split(A, *indices):
    """
    returns a list of rows corresponding to the passed indices
    """
    if not indices:
        indices = range(A.shape[0])
    res = [ A[i, :] for i in indices ]
    return res



# TODO: Doctest
def concat_rows(*args):
    """
    takes some row (hyper-)vectors and aggregetes them to a matrix
    """

    row_list = []

    for a in args:
        if a.shape[0] == 1:
            row_list.append( list(a) )
            continue
        for i in xrange(a.shape[0]):
            row_list.append( list(a[i, :]) )
    m = sp.Matrix(row_list)

    return m

# other name:
row_stack = concat_rows

# geschrieben für Polynommatritzen

def col_minor(A, *cols, **kwargs):
    """
    returns the minor (determinant) of the columns in cols
    """
    n, m = A.shape

    method = kwargs.get('method', "berkowitz")

    assert m >= n
    assert len(cols) == n

    M = sp.zeros(n, n)
    for i, idx in enumerate(cols):
        M[:, i] = A[:, idx]

    return M.det(method = method).expand()


def general_minor(A, rows, cols, **kwargs):
    """
    selects some rows and some cols of A and returns the det of the resulting
    Matrix
    """

    method = kwargs.get('method', "berkowitz")

    Q = row_col_select(A, rows, cols)

    return Q.det(method = method).expand()


def all_k_minors(M, k, **kwargs):
    """
    returns all minors of order k of M

    Note that if k == M.shape[0]

    this computes all "column-minors"
    """
    m, n = M.shape

    assert k<= m
    assert k<= n

    row_idcs = list(it.combinations(range(m), k))
    col_idcs = list(it.combinations(range(n), k))

    rc_idx_tuples = list(it.product(row_idcs, col_idcs))

    method = kwargs.get('method', "berkowitz")

    res = []
    for rr, cc in rc_idx_tuples:
        res.append(general_minor(M, rr, cc, method = method))

    return res

def row_col_select(A, rows, cols):
    """
    selects some rows and some cols of A and returns the resulting Matrix
    """

    Q1 = sp.zeros(A.shape[0], len(cols))
    Q2 = sp.zeros(len(rows), len(cols))

    for i, c in enumerate(cols):
        Q1[:, i] = A[:, c]


    for i, r in enumerate(rows):
        Q2[i, :] = Q1[r, :]

    return Q2




def col_select(A, *cols):
    """
    selects some columns from a matrix
    """
    Q = sp.zeros(A.shape[0], len(cols))

    for i, c in enumerate(cols):
        Q[:, i] = A[:, c]

    return Q



#def is_left_coprime2(Ap, eps=1e-1):
#    """
#    Alternativer Zugang
#    TODO: merge
#    """

#
#    nonzero_numbers = []
#    root_list = []
#    for m in minors:
#        if m.is_number and m != 0:
#            nonzero_numbers.append(m)
#            # we could return true already here!!
#        else:
#            root_list.append(roots(m))
#
#    if len(nonzero_numbers) > 0:
#        return True
#
#    # test is there a common zero in all root-lists?






def is_left_coprime(Ap, Bp=None, eps = 1e-10):
    """
    Test ob Ap,Bp Linksteilerfrei sind
    keine Parameter zulässig

    """

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

    symbs = list(matrix_atoms(M, sp.Symbol))
    assert len(symbs) == 1
    symb = symbs[0]

    combinations = it.combinations(range(r+m), r)

    minors = [col_minor(M, *cols) for cols in combinations]

    nonzero_const_minors = [m for m in minors if (m !=0) and (symb not in m)]

    if len(nonzero_const_minors) > 0:
        return True

    #zero_minors = [m for m in minors if m == 0]

    # polymionors (rows belong together)
    all_roots = [roots(m) for m in minors if symb in m]

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



def get_expr_var(expr, var = None):
    """
    auxillary function
    if var == None returns the unique symbol which is contained in expr:
    if no symbol is found, returns None
    """
    expr = sp.sympify(expr)
    if not var == None:
        assert isinstance(var, sp.Symbol)
        return var
    else: # var == None
        symbs = list(expr.atoms(sp.Symbol))
        if len(symbs) == 0:
            return None
        elif len(symbs) == 1:
            return symbs[0]
        else:
            errmsg = "%s contains more than one variable: %s " % (expr, symbs)
            raise ValueError, errmsg


def poly_degree(expr, var=None):
    """
    returns degree of monovariable polynomial
    """
    var = get_expr_var(expr, var)
    if var == None:
        return sp.sympify(0)

    P = sp.Poly(expr, var, domain = "EX")
    return P.degree()


def poly_coeffs(expr, var=None):
    """
    returns all (monovariate)-poly-coeffs (including 0s) as a list
    first element is highest coeff.
    """
    var = get_expr_var(expr, var)
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



def rat_if_close(x, tol=1e-10):
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

def matrix_with_rationals(A):
    A = sp.Matrix(A)

    def rat(x):
        y = fractionfromfloat(x)
        return sp.Rational(y.numerator, y.denominator)

    A2 = A.applyfunc(rat)

    # error

    diff = A-A2
    a_diff = np.abs(to_np(diff))

    A3 = np.array(A2) # (dtype=object)

    res = np.where(a_diff < 1e-10, A3, to_np(A))

    return sp.Matrix(res)

arr_float = np.frompyfunc(np.float, 1,1)

def to_np(arr, dtype=np.float):
    """ converts a sympy matrix in a nice numpy array
    """
    if isinstance(arr, sp.Matrix):
        symbs = list(matrix_atoms(arr, sp.Symbol))
        assert len(symbs) == 0, "no symbols allowed"

    # because np.int can not understand sp.Integer
    # we temporarily convert to float
    # TODO: make this work with complex numbers..
    arr1 = arr_float( np.array(arr) )
    return np.array( arr1, dtype )

def roots(expr):
    import scipy as sc
    return sc.roots(coeffs(expr))

def real_roots(expr):
    import scipy as sc
    r = sc.roots(coeffs(expr))
    return np.real( r[np.imag(r)==0] )

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


def fac(i):
    # see also sp.factorial
    if i == 0:
        return 1
    return i * fac(i-1)

def div(vf, x):
    """divergence of a vector field"""
    vf = list(vf)
    x = list(x)
    assert len(vf) == len(x)

    return sum([c.diff(xi) for c,xi in zip(vf, x)])

def chop(expr, tol = 1e-10):
    """suppress small numerical values"""

    expr = sp.expand(sp.sympify(expr))
    if expr.is_Symbol: return expr

    assert expr.is_Add

    return sp.Add(*[term for term in expr.as_Add() if sp.abs(term.as_coeff_terms()[0]) >= tol])



# trigsimp


def trigsimp2(expr):
    """
    sin**2 + cos**2 = 1 in big expressions
    """

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
        poly = trig_term_poly(expr, s)

        dd = poly.as_dict()


        uncommon_coeff = (dd[(2,0)] - dd[(0,2)]).expand()



        if uncommon_coeff == 0:
            expr += dd[(2,0)] - dd[(2,0)]*sp.sin(s)**2 - dd[(2,0)]*sp.cos(s)**2


        print dd[(2,0)]

    return expr.expand()




def rev_tuple(tup):
    return [(t[1], t[0]) for t in tup]

def gradient(scalar_field, xx):
    # returns a row vector (coverctorfiel)!
    return sp.Matrix([scalar_field]).jacobian(xx)

def trig_term_poly(expr, s):
    """
    s ... the argument of sin, cos
    """

    X, Y = sp.symbols('tmpX_, tmp_Y')

    poly = sp.Poly( expr.subs([(sp.sin(s), X), (sp.cos(s),Y)]), X,Y, domain='EX')

    return poly


def matrix_atoms(M, *args, **kwargs):
    sets = [m.atoms(*args, **kwargs) for m in list(M)]
    S = set().union(*sets)

    return S

def atoms(expr, *args, **kwargs):
    if isinstance(expr, (sp.Matrix, list)):
        return matrix_atoms(expr, *args, **kwargs)
    else:
        return expr.atoms(*args, **kwargs)


def matrix_count_ops(M, visual=False):
    def co(expr):
        return sp.count_ops(expr, visual)
    return M.applyfunc(co)




#def get_diffterms2(xx, order):
#    assert len(xx) == 2
#    order = int(order)





def get_diffterms(xx, order):
    """
    returns a list such as

    [(x1, x1), (x1, x2), (x1, x3), (x2, x2), (x2, x3), (x3, x3)]

    for xx = (x1, x2, x3) and order = 2

    """
    if order == 0:
        return []

    if len(xx) == 2:
        return [ (xx[0],)*(order-i)+(xx[1],)*(i) for i in range(order+1)]

    if isinstance(order, (list, tuple)):
        return sum([get_diffterms(xx, o) for o in order], [])

    assert isinstance(order, int)

    L1 = list(  it.product( *([xx]*order) )   )
    #L2 = map(list, L1)
    L3 = map(tuple, map(sorted, L1))
    terms = dict(zip(L3, [0]*len(L3))).keys() # remove duplicates
    terms.sort()

    return terms



def multi_series(expr, xx, order, poly=False):
    """
    Reihenentwicklung (um 0) eines Ausdrucks in mehreren Variablen
    """

    xx0 = zip(xx, [0]*len(xx)) # Entwicklungsstelle
    res = 0
    for i in range(order+1):
        if i == 0:
            res += expr.subs(xx0)
            continue
        terms = get_diffterms(xx, i)
        for tup in terms:
            cnt = Counter(tup) # returns a dict
            fac_list = [sp.factorial(n) for n in cnt.values()]
            #fac = 1/sp.Mul(*fac_list)
            res += expr.diff(*tup).subs(xx0)*sp.Mul(*tup) / sp.Mul(*fac_list)

    if poly:
        res = sp.Poly(res, *xx, domain="EX")
    return res


def matrix_series(m, xx, order, poly=False):

    assert isinstance(m, sp.Matrix)
    def appfnc(expr):
        return multi_series(expr, xx, order, poly)

    return m.applyfunc(appfnc)


def linear_input_trafo(B, row_idcs):
    """
    serves to decouple inputs from each other
    """
    B = sp.Matrix(B)

    n, m = B.shape

    assert len(row_idcs) == m
    assert len(set(row_idcs)) == m

    Bnew = sp.zeros(m,m)
    for i, idx in enumerate(row_idcs):
        Bnew[i, :] = B[idx, :]

    P = symbMatrix(m,m)

    leqs = Bnew*P-sp.eye(m)
    res = sp.solve(leqs)

    return to_np(P.subs(res))





def poly_scalar_field(xx, symbgen, order, poly=False):
    """
    returns a multivariate poly with specified oders
    and symbolic coeffs
    returns also a list of the coefficients
    """

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

        terms = get_diffterms(xx, i)

        for tup in terms:
            c = symbgen.next()
            res += c*sp.Mul(*tup)
            coeff_list.append(c)

    if poly:
        res = sp.Poly(res, *xx, domain='EX')
    return res, coeff_list

def np_trunc_small_values(arr, lim = 1e-10):
    assert isinstance(arr, (np.ndarray, np.matrix))

    bool_abs = np.abs(arr) < lim

    res = arr*1
    res[bool_abs] = 0
    return res


def trunc_small_values(expr, lim = 1e-10, n=1):
    expr = ensure_mutable( sp.sympify(expr) )

    a_list = list(atoms(expr, sp.Number))
    subs = []
    for a in a_list:
        if sp.Abs(a) < lim:
            subs.append((sp.Abs(a), 0))
            # substituting Abs(a) circumvents Problems in terms like sin(...)
            if a < 0:
                subs.append((a, 0))

    res = expr.subs(subs)

    if n <= 1:
        return res
    else:
        return trunc_small_values(res, lim, n-1)



def clean_numbers(expr, eps=1e-10):
    """
    trys to clean all numbers from numeric noise
    """

    if isinstance(expr, (list, tuple)):
        return [clean_numbers(elt, eps) for elt in expr]

    expr = trunc_small_values(expr)

    maxden = int(1/eps)
    floats = list(atoms(expr, sp.Float))
    rats = []
    dummy_symbs = []
    symb_gen = sp.numbered_symbols('cde', cls = sp.Dummy)
    for f in floats:
        rat = sp_fff(f, maxden)
        rats.append(rat)
        dummy_symbs.append(symb_gen.next())

    res1 = expr.subs(zip(floats, dummy_symbs))
    res2 = res1.subs(zip(dummy_symbs, rats))

    return res2



def random_equaltest(exp1, exp2,  info = False, integer = False, seed = None, tol = 1e-14, min=-1, max=1):
    """
    serves to check numerically (with random numbers) whether exp1, epx2 are equal
    # TODO: unit test
    """

    if isinstance(exp1, sp.Matrix):
        assert isinstance(exp2, sp.Matrix)
        assert exp1.shape == exp2.shape, "Different shape"
        m,n = exp1.shape

        def func(exp1, exp2):
            return random_equaltest(exp1, exp2, info, integer, seed,
                                                            tol, min, max)
        res = [func(e1, e2) for e1,e2 in zip(list(exp1), list(exp2))]

        if info == True:
            res = [tup[1] for tup in res]
        return sp.Matrix(res).reshape(m,n)

    exp1 = sp.sympify(exp1)
    exp2 = sp.sympify(exp2)


    gen = sp.numbered_symbols('ZZZ', cls=sp.Dummy)

    derivs = exp1.atoms(sp.Derivative).union(exp2.atoms(sp.Derivative))

    # we only look for undefined functions, but not for cos(..) etc
    # the latter would be matched by sp.Function
    func_class = sp.function.AppliedUndef

    funcs = exp1.atoms(func_class).union(exp2.atoms(func_class))

    # replace all functions and derivs with symbols
    SL = [(X, gen.next()) for X in list(derivs)+list(funcs)]


    exp1 = exp1.subs(SL)
    exp2 = exp2.subs(SL)



    a1 = exp1.atoms(sp.Symbol, sp.Dummy)
    a2 = exp2.atoms(sp.Symbol, sp.Dummy)

    r = random
    if seed != None:
        r.seed(seed)

    def get_rand():
        if not integer:
            return (r.random()*(max-min)+min)
        else:
            return r.randint(min, max)

    tuples = [(s, get_rand()) for s in a1.union(a2)]

    if not integer:
        diff = exp1.subs(tuples).evalf() - exp2.subs(tuples).evalf()
    else:
        diff = exp1.subs(tuples) - exp2.subs(tuples)


    if info == False:
        return abs(diff) <= tol
    else:
        return abs(diff) <= tol, diff



def matrix_random_equaltest(M1, M2,  info=False, **kwargs):
    raise DeprecationWarning, "use random_equaltest instead"


# TODO: Funktionen und Ableitungen (aus random_equaltest rauslösen
# und auch hier verwenden )
def subs_random_numbers(expr):
    """
    replaces all symbols in the given expr (scalar or matrx) by random numbers
    and returns the substituted result

    usefull for e.g. for checking the rank at a "generic" point
    """

    a = atoms(expr, sp.Symbol)

    tuples = [(s, random.random()) for s in a]

    return expr.subs(tuples)

def matrix_random_numbers(M):
    raise DeprecationWarning, "use subs_random_numbers"


def zip0(*xx, **kwargs):
    """
    returns a list of tuples like: [(x1, arg), (x2, arg), ...]
    this is very handy for substiution of equilibrium points at zero

    valid kwargs: arg (if a different value than zero is desired)
    """
    arg = kwargs.get('arg', 0)

    # if a sequence was passed instead of single arguments
    if hasattr(xx[0], '__len__'):
        xx = xx[0]

    return zip(xx, [arg]*len(xx))

def aux_make_tup_if_necc(arg):
    """
    checks whether arg is iterable.
    if not return (arg,)
    """
    if not hasattr(arg, '__len__'):
        return (arg,)

    return arg

def expr_to_func(args, expr, modules = 'numpy', **kwargs):
    """
    wrapper for sympy.lambdify to handle constant expressions
    (shall return a numpyfied function as well)

    this function bypasses the following problem:

    f1 = sp.lambdify(t, 5*t, modules = "numpy")
    f2 = sp.lambdify(t, 0*t, modules = "numpy")

    f1(np.arange(5)).shape # -> array
    f2(np.arange(5)).shape # -> int


    Some special kwargs:
    np_wrapper == True:
        the return-value of the resulting function is passed through
        to_np(..) before returning

    eltw_vectorize: allows to handle vectors of piecewise expression

    """

    # TODO: sympy-Matrizen mit Stückweise definierten Polynomen
    # numpy fähig (d.h. vektoriell) auswerten

    # TODO: Unittest

    expr = sp.sympify(expr)
    expr = ensure_mutable(expr)
    expr_tup = aux_make_tup_if_necc(expr)
    arg_tup = aux_make_tup_if_necc(args)

    new_expr = []
    arg_set = set(arg_tup)
    for e in expr_tup:
        assert isinstance(e, sp.Expr)
        # args (Symbols) which are not in that expression
        diff_set = arg_set.difference(e.atoms(sp.Symbol))

        # add and subtract the respective argument such that it occurs
        # without changing the result
        for d in diff_set:
            assert isinstance(d, sp.Symbol)
            e = sp.Add(e, d, -d, evaluate = False)

        new_expr.append(e)

    if not hasattr(expr, '__len__'):
        assert len(new_expr) == 1
        new_expr = new_expr[0]


    # extract kwargs specific for lambdify
    printer = kwargs.get('printer', None)
    use_imps = kwargs.get('use_imps', True)
    func = sp.lambdify(args, new_expr, modules, printer, use_imps)


    if kwargs.get('np_vectorize', False):
        func1 = np.vectorize(func)
    else:
        func1 = func

    if kwargs.get('special_vectorize', False):
        # !! this was introduced for a very special application
        # should be removed from the general code
        def func2(*allargs):
            return to_np(func(*allargs))

        f = np.float
        func3 = np.vectorize(func2, otypes = [f,f,f, f,f,f])
        return func3

    if kwargs.get('eltw_vectorize', False):
        # elementwise vectorization to handle piecewise expressions
        # each function returns a 1d-array
        assert len(new_expr) >=1 # elementwise only makes sense for sequences
        funcs = []
        for e in new_expr:
            func_i = sp.lambdify(args, e, modules, printer, use_imps)
            func_iv = np.vectorize(func_i, otypes=[np.float])
            funcs.append(func_iv)

        def func2(*allargs):
            # each result should be a 1d- array
            results = [to_np(f(*allargs)) for f in funcs]

            # transpose, such that the input axis (e.g. time) is the first one
            return to_np(results).T.squeeze()

        return func2

    if kwargs.get('np_wrapper', False):
        def func2(*allargs):
            return to_np(func1(*allargs))
    elif kwargs.get('list_wrapper', False):
        def func2(*allargs):
            return list(func1(*allargs))
    else:
        func2 = func1
    return func2

def ensure_mutable(arg):
    """
    ensures that we handle a mutable matrix (iff arg is a matrix)
    """
    # TODO: e.g. sp.sympify converts a MutableMatrix to ImmutableMatrix
    # maybe this changes in future sympy releases
    # which might make this function obsolete (?)
    if isinstance(arg, sp.matrices.MatrixBase):
        return as_mutable_matrix(arg)
    else:
        return arg




def as_mutable_matrix(matrix):
    """
    sympy sometimes converts matrices to immutable objects
    this can be reverted by a call to    .as_mutable()
    this function provides access to that call as a function
    (just for cleaner syntax)
    """
    return matrix.as_mutable()

def is_col_reduced(A, symb, return_internals = False):
    """
    tests whether polynomial Matrix A is column-reduced

    optionally returns internal variables:
        the list of col-wise max degrees
        the matrix with the col.-wise-highest coeffs (Gamma)

    Note: concept of column-reduced matrix is important e.g. for
    solving a Polynomial System w.r.t. highest order "derivative"

    Note: every matrix can be made col-reduced by unimodular transformation
    """
    Gamma = as_mutable_matrix(A*0)
    n, m = A.shape

    assert n == m

    A = trunc_small_values(A)

    # degrees:
    A_deg = to_np(matrix_degrees(A, symb), dtype = np.float)
    max_degrees = list(A_deg.max(axis=0)) # columnwise maximum

    # TODO: unit-Test
    # handle zero columns:
    infty = float(sp.oo)
    max_degrees = [int(md) if not md == -infty else md for md in max_degrees]

    # maximum coeffs:
    for j in range(m):
        deg = max_degrees[j]
        for i in range(n):
            Gamma[i,j] = get_order_coeff_from_expr(A[i,j], symb, deg)

    result = Gamma.rank() == m
    if return_internals:
        # some functions might need this information
        internals = Container(Gamma = Gamma, max_degrees = max_degrees)
        return result, internals
    else:
        return result

def is_row_reduced(A, symb, *args, **kwargs):
    """
    transposed Version of is_col_reduced(...)
    """
    res = is_col_reduced(A.T, symb, *args, **kwargs)
    if isinstance(res, tuple):
        C = res[0]
        C.Gamma = C.Gamma.T
    return res


def get_col_reduced_right(A, symb, T = None, return_internals = False):
    """
    Takes a polynomial matrix A(s) and returns a unimod Transformation T(s)
    such that   A(s)*T(s) (i.e. right multiplication) is col_reduced.

    Approach is taken from appendix of the PHD-Thesis of S. O. Lindert (2009)

    :args:
        A:  Matrix
        s:  symbol
        T:  unimod-Matrix from preceeding steps

    -> recursive approach

    :returns:
        Ar: reduced Matrix
        T:  unimodular transformation Matrix
    """

    n, m = A.shape
    assert n == m

    if T == None:
        T = sp.eye(n)
    else:
        assert T.shape == (n, m)
        d = T.berkowitz_det().expand()
        assert d != 0 and not symb in d


    A_work = trunc_small_values(sp.expand(A*T))


    cr_flag, C = is_col_reduced(A_work, symb, return_internals = True)

    # C.Gamma is the matrix with col-wise highest coeff
    if cr_flag:
        # this is the only exit point
        res = A_work.expand(), T
        if return_internals:
            res += (C,)
        return res
    else:
        pass
        # C.Gamma is nonregular

    g = C.Gamma.nullspace()[0]
    non_zero_cols_IDX = to_np(g).flatten() != 0
    # get the max_degrees wrt. to each non-zero component of g
    non_zero_cols_degrees = to_np(C.max_degrees)[non_zero_cols_IDX]

    N = max(non_zero_cols_degrees)
    # construct the diagonal matrix
    diag_list = []
    for i in range(m):
        cd = col_degree(A_work[:, i],symb)
        diag_list.append( symb**int(N-cd) )

    # gamma_p:
    gp = sp.diag(*diag_list)*g


    T1 = unimod_completion(gp, symb)

    TT = trunc_small_values( sp.expand(T*T1) )

    # recall this method with new T

    return get_col_reduced_right(A, symb, TT, return_internals)



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

def matrix_degrees(A, symb):

    element_deg = element_deg_factory(symb)

    return A.applyfunc(element_deg)

def col_degree(col, symb):
    return max(matrix_degrees(col, symb))

def unimod_completion(col, symb):
    """
    takes a column and completes it such that the result is unimodular
    """

    # there must at least one nonzero constant in col:

    n, m = col.shape
    assert m == 1
    element_deg = element_deg_factory(symb)

    idx = None
    for i, c in enumerate(list(col)):
        if c != 0 and element_deg(c) == 0:

        # we want the index of the first non-zero const. of col
            idx = i
            break

    assert not idx == None, "there should have been a nonzero const."


    T = sp.eye(n)

    T[:, idx] = col

    return T







def subs_same_symbs(expr, new_symbs):
    """
    subs_same_symbs(x+y, [x, y])
    returns x+y, where the symbols are taken from the list
    (symbs in exp might be different objects with the same name)

    this functions helps if expr comes from a string

    """

    old_symbs = list(atoms(expr, sp.Symbol))

    string_dict = dict([(s.name, s) for s in new_symbs])


    subs_list = [ (s, string_dict[s.name]) for s in old_symbs]

    return expr.subs(subs_list) # replpace new symbs by old ones



def symm_matrix_to_vect(M):
    """ converts
     a b
     b c
            to      [a, b, c]
    """

    n,m = M.shape
    assert m == n
    assert M == M.T

    res = sp.zeros(int(n+.5*n*(n-1)), 1)
    k = 0
    for i in range(n):
        for j in range(i, n):
            if i == j:
                val = M[i,j]
            else:
                val = 2*M[i,j]
            res[k] = val
            k+=1

    return res


def pinv(M):
    """compute symbolic pseudo inverse"""
    # http://stackoverflow.com/questions/15426624/computing-pseudo-inverse-of-a-matrix-using-sympy
    M = sp.Matrix(M)
    n1, n2 = M.shape
    if n2 > n1:
        # more columns -> assume full row rank
        # right inverse
        pinv = M.T * (M * M.T).inv()
        res = M*pinv

    else:
        #assume full column rank
        # left inverse
        pinv = (M.T * M).inv() * M.T
        res = pinv*M

    #print res
    return pinv

def nullspaceMatrix(M, *args, **kwargs):
    """
    wrapper for the sympy-nullspace method
    returns a Matrix where each column is a basis vector of the nullspace
    additionally it uses the enhanced nullspace function to calculate
    ideally simple (i.e. fraction free) expressions in the entries
    """

    n = enullspace(M, *args, **kwargs)
    return col_stack(*n)


#todo: (idea) extend penalty to rational and complex numbers
def _extended_count_ops(expr):
    """
    extended count ops, which penalizes symbols over numbers
    """
    expr = sp.sympify(expr)
    res = sp.count_ops(expr)
    if res > 0:
        return res
    if expr.is_Symbol:
        return res+.5
    return res


# todo: unit test-example P = sp.Matrix([-x2, -x1, 1]).T (-> )
def enullspace(M, *args, **kwargs):
    """
    enhanced nullspace: calculate basis vectors with ideally simple
    (i.e. fraction free) expressions in the entries

    :kwargs:
        simplify = True (apply simplification before examining denominators)

    """

    vectors = M.nullspace(*args, **kwargs)

    if kwargs.get('simplify', True):
        vectors = [sp.simplify(v) for v in vectors]


    new_vectors = []
    for v in vectors:
        # get the "most complicated" denomitator of the coordinates
        denoms = [ c.as_numer_denom()[1] for c in v]
        #denom_tuples = [(d.count_ops(), d) for d in denoms]
        #denom_tuples.sort()
        denoms.sort(key = _extended_count_ops)

        d = denoms[-1]
        # convert to mutable matrix
        res = sp.Matrix(v*d)
        new_vectors.append( res )

    return new_vectors


def vect_to_symm_matrix(v):
    """ converts
     [a, b, c]

     to    a b
           b c
    """

    v = sp.Matrix(list(v))
    L, m = v.shape
    assert m == 1
    n = -.5 + sp.sqrt(.25+2*L)

    if not int(n) == n:
        raise ValueError, "invalid length"
    n = int(n)

    M = sp.zeros(n,n)
    k = 0
    for i in range(n):
        for j in range(i, n):
            if i == j:
                M[i,j] = v[k]
            else:
                M[i,j] = v[k]/2
                M[j,i] = v[k]/2
            k+=1

    return M


# ! eigentlich eher numtools

def dd(*args):
    """
    dd(a,b,c, ...) = np.dot(a, np.dot(b, np.dot(c, ...)))
    """
    return reduce(np.dot, args)



## !! Laplace specific

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


def simplify_derivs(expr):
    """
    iterates over Derivative-atoms and performs "doit"
    (workarround for a sympy functionality)
    """
    A = list(atoms(expr, sp.Derivative))

    def key_fnc(a):
        """
        get the length of the args to determine the order of the derivative
        """
        return len(a.args)

    # highest first
    A.sort(key = key_fnc, reverse=True)

    SUBS = [(a, a.doit()) for a in A]
#    print SUBS
    return expr.subs(SUBS)


def perform_time_derivative(expr, func_symbols, t_symbol=None, order=1):
    """
    Example: expr = f(a, b). We know that a, b are time-functions: a(t), b(t)
    we want : expr.diff(t) with te appropriate substitutions made
    :param expr:
    :param func_symbols:
    :return:
    """

    if not t_symbol:
        t = sp.Symbol("t")
    else:
        t = t_symbol

    func_symbols = list(func_symbols) # we work with lists here

    funcs = [ symbs_to_func(s, [s], t) for s in func_symbols ]

    derivs1 = [[f.diff(t, ord) for f in funcs] for ord in range(order, 0, -1)]


    def extended_name_symb(base, ord):
        if isinstance(base, sp.Symbol):
            base = base.name
        if base.endswith('_d') or base.endswith('_dd'):
            sep = ''
        else:
            sep = '_'
        new_name = base + sep + 'd'
        if ord == 1:
            return sp.Symbol(new_name)
        else:
            return extended_name_symb(new_name, ord - 1)

    deriv_symbols1 = [ [extended_name_symb(s, ord)
                      for s in func_symbols] for ord in range(order, 0, -1)]

    # flatten the lists:
    derivs = []
    for d_list in derivs1:
        derivs.extend(d_list)

    deriv_symbols = []
    for ds_list in deriv_symbols1:
        deriv_symbols.extend(ds_list)

    subs1 = zip(func_symbols, funcs)

    # important: begin substitution with highest order
    subs2 = zip(derivs + funcs, deriv_symbols + func_symbols)

    expr1 = expr.subs(subs1)
    expr2 = expr1.diff(t, order)
    expr3 = expr2.subs(subs2)

    return expr3


# TODO: (multipl. by integers (especially -1) should be preserved by default)
def symbolify_matrix(M):
    """
    associates a symbol to every matrix entry
    respects equal values, keeps atoms unchanged


    Example: f1, ... f4 are expressions, a1 is a symbol

    Input:

    f1+f2       a1         0    7
     a1     f3*sin(f4)   f1+f2  a1

    Output:

    [(x0, f1+f2), (x1, f1+f2)]

    x0 a1 0  7
    a1 x1 x0 a1
    """
    assert isinstance(M, sp.Matrix)
    L = list(M)

    gen = sp.numbered_symbols(cls = sp.Dummy)
    replaced = []
    new_symbol = []
    result = []


    for e in L:
        if e.is_Atom:
            # leave the entry unchanged
            ns = e
        elif e in replaced:
            # get the old symbol
            ns = new_symbol[replaced.index(e)]
        else:
            replaced.append(e)
            # create a new symbol
            ns = gen.next()
            new_symbol.append(ns)

        result.append(ns)

    res = sp.Matrix(result).reshape(*M.shape)
    replacements = zip(new_symbol, replaced)
    return replacements, res





# eigene Trigsimp-Versuche
# mit der aktuellen sympy-Version (2013-03-29)
# eigentlich überflüssig
# TODO: in Schnipsel-Archiv überführen

def sort_trig_terms(expr):

    expr = expr.expand()
    assert type(expr) == sp.Add

    trig_terms = list(expr.atoms(sp.sin, sp.cos))

    res = {}

    for a in expr.args:
        coeff = a.subs(zip(trig_terms, [1]*len(trig_terms)))


        sign = 1

        # Vorzeichen in auf die Funktion verschieben
        if str(coeff)[0] == "-":
            sign = -1
            coeff*=-1
        else:
            sign = 1

        L = res.get(coeff)
        if L == None:
            L = []
            res[coeff] = L
        L.append( a.subs(coeff, 1) ) # minus is already considered

    return res


def simp_trig_dict(sdict):
    """
    takes a sorted dict, simplifies each value and adds all up
    """

    items = sdict.items()

    res = 0
    for k, v in items:
        res += k*sp.trigsimp(sum(v))

    return res


def my_trig_simp(expr):
    """
    ersetzt größere argumente in trigonom funktionen durch einfache Symbole
    """

    trig_terms = list(expr.atoms(sp.sin, sp.cos))


    args = []

    for tt in trig_terms:
        args.append(tt.args[0])

    args = list(set(args))
    args.sort(key = sp.count_ops)

    symbs = sp.symbols( 'x1:%i' %( len(args)+1) )

    subslist = zip(args, symbs)
    subslist.reverse()

    return expr.subs(subslist), subslist









