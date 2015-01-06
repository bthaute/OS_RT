# -*- coding: utf-8 -*-
"""
.. module:: matrix
    :synopsis: Functions concerning operations on matrices.
    
.. moduleauthor:: Carsten Knoll
"""
#from __future__ import absolute_import

import sympy as sp
import numpy as np

import itertools as it
import random

import pycontroltools.auxfuncs.programming.miscprog

import miscmath
import numtools
import polynomial

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
    if isinstance(arg, sp.Matrix):
        return arg.applyfunc(sp.ratsimp)
    else:
        return sp.ratsimp(arg)
        
def symbMatrix(n, m, s='a', symmetric = 0):
    """
    """
    A = sp.Matrix(n,m, lambda i,j:sp.Symbol( s+'%i%i'%(i+1,j+1)) )
    if symmetric == 1:
        subs_list = symmetryDict(A)
        A = A.subs(subs_list)
    return A
    
#chris: gibt eine Matrix tmp (Abbild der übergebenen Matrix) zurück
#       in der alle belegten Elemente (!=0) mit 1 und alle unbelegten Elemente (=0)
#       mit 0 ersetzt sind
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

def mdiff(M, var):
    """
    returns the elementwise derivative of a matrix M w.r.t. var
    """
    return M.applyfunc(lambda elt: sp.diff(elt, var))
    
def get_rows(A):
    """
    returns a list of n x 1 vectors
    """
    A = sp.Matrix(A)
    n, m = A.shape

    return [A[:,i] for i in range(m)]
    
def elementwise_mul(M1, M2):
    """
    performs elment wise multiplication of matrices
    """
    assert M1.shape == M2.shape
    return sp.Matrix(np.array(M1)*np.array(M2))
    
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

# TODO: Doctest
#chris: selbe funktion wie np.hstack ?
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

# TODO: Doctest
#chris: selbe funktion wie np.vstack ?
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
    
def matrix_with_rationals(A):
    """
    """
    A = sp.Matrix(A)

    def rat(x):
        y = miscmath.fractionfromfloat(x)
        return sp.Rational(y.numerator, y.denominator)

    A2 = A.applyfunc(rat)

    # error

    diff = A-A2
    a_diff = np.abs(numtools.to_np(diff))

    A3 = np.array(A2) # (dtype=object)

    res = np.where(a_diff < 1e-10, A3, numtools.to_np(A))

    return sp.Matrix(res)

arr_float = np.frompyfunc(np.float, 1,1)

def matrix_atoms(M, *args, **kwargs):
    """
    """
    sets = [m.atoms(*args, **kwargs) for m in list(M)]
    S = set().union(*sets)

    return S

def matrix_count_ops(M, visual=False):
    """
    """
    def co(expr):
        return sp.count_ops(expr, visual)
    return M.applyfunc(co)
    
def matrix_series(m, xx, order, poly=False):
    """
    """
    assert isinstance(m, sp.Matrix)
    def appfnc(expr):
        return miscmath.multi_series(expr, xx, order, poly)

    return m.applyfunc(appfnc)

def matrix_random_equaltest(M1, M2,  info=False, **kwargs):
    """
    """
    n1, n2 = M1.shape
    assert M2.shape == (n1, n2)

    m1 = list(M1)
    m2 = list(M2)

    res = [numtools.random_equaltest(e1, e2, info, **kwargs) for e1, e2 in zip(m1, m2)]

    if info:
        res1, res2 = zip(*res)

        res1 = sp.Matrix(res1).reshape(n1, n2)
        res2 = sp.Matrix(res2).reshape(n1, n2)
        res = res1, res2
    else:
        res = sp.Matrix(res).reshape(n1, n2)

    return res

def matrix_subs_random_numbers(M):
    """
    substitute every symbol in M with a random number

    this might be usefull to determine the generic rank of a matrix
    """
    miscprog = pycontroltools.auxfuncs.programming.miscprog

    a = miscprog.atoms(M, sp.Symbol)
    tuples = [(s, random.random()) for s in a]

    return M.subs(tuples)
    
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
    miscprog = pycontroltools.auxfuncs.programming.miscprog

    Gamma = as_mutable_matrix(A*0)
    n, m = A.shape

    assert n == m

    A = numtools.trunc_small_values(A)

    # degrees:
    A_deg = numtools.to_np(matrix_degrees(A, symb), dtype = np.int32)
    max_degrees = list(A_deg.max(axis=0)) # columnwise maximum

    # maximum coeffs:
    for j in range(m):
        deg = max_degrees[j]
        for i in range(n):
            Gamma[i,j] = polynomial.get_order_coeff_from_expr(A[i,j], symb, deg)

    result = Gamma.rank() == m
    if return_internals:
        # some functions might need this information
        internals = miscprog.Container(Gamma = Gamma, max_degrees = max_degrees)
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


    A_work = numtools.trunc_small_values(sp.expand(A*T))


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
    non_zero_cols_IDX = numtools.to_np(g).flatten() != 0
    # get the max_degrees wrt. to each non-zero component of g
    non_zero_cols_degrees = numtools.to_np(C.max_degrees)[non_zero_cols_IDX]

    N = max(non_zero_cols_degrees)
    # construct the diagonal matrix
    diag_list = []
    for i in range(m):
        cd = col_degree(A_work[:, i],symb)
        diag_list.append( symb**int(N-cd) )

    # gamma_p:
    gp = sp.diag(*diag_list)*g


    T1 = unimod_completion(gp, symb)

    TT = numtools.trunc_small_values( sp.expand(T*T1) )

    # recall this method with new T

    return get_col_reduced_right(A, symb, TT, return_internals)
    
def matrix_degrees(A, symb):
    """
    """
    element_deg = polynomial.element_deg_factory(symb)

    return A.applyfunc(element_deg)

def col_degree(col, symb):
    """
    """
    return max(matrix_degrees(col, symb))
    
def unimod_completion(col, symb):
    """
    takes a column and completes it such that the result is unimodular
    """

    # there must at least one nonzero constant in col:

    n, m = col.shape
    assert m == 1
    element_deg = polynomial.element_deg_factory(symb)

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