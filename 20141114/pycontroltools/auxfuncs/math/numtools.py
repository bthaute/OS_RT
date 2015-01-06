# -*- coding: utf-8 -*-
"""
.. module:: numtools
    :synopsis: Numerical tools.
    
.. moduleauthor:: Carsten Knoll
"""
import sympy as sp
import numpy as np
import scipy as sc
import scipy.integrate


import pycontroltools.auxfuncs.programming.miscprog

import matrix
import miscmath
import random

arr_float = np.frompyfunc(np.float, 1,1)
def to_np(arr, dtype=np.float):
    """ converts a sympy matrix in a nice numpy array
    """
    if isinstance(arr, sp.Matrix):
        symbs = list(matrix.matrix_atoms(arr, sp.Symbol))
        assert len(symbs) == 0, "no symbols allowed"

    # because np.int can not understand sp.Integer
    # we temporarily convert to float
    # TODO: make this work with complex numbers..
    arr1 = arr_float( np.array(arr) )
    return np.array( arr1, dtype )
    
def chop(expr, tol = 1e-10):
    """suppress small numerical values"""

    expr = sp.expand(sp.sympify(expr))
    if expr.is_Symbol: return expr

    assert expr.is_Add

    return sp.Add(*[term for term in expr.as_Add() if sp.abs(term.as_coeff_terms()[0]) >= tol])
    
def np_trunc_small_values(arr, lim = 1e-10):
    """
    """
    assert isinstance(arr, (np.ndarray, np.matrix))

    bool_abs = np.abs(arr) < lim

    res = arr*1
    res[bool_abs] = 0
    return res
    
#def np_trunc_small_values(A, eps = 1e-10):
#
#    res = A*1 # leave the original untouched
#    res[np.abs(res) < eps] = 0
#
#    return res
    
def trunc_small_values(expr, lim = 1e-10, n=1):
    
    miscprog = pycontroltools.auxfuncs.programming.miscprog
    
    expr = matrix.ensure_mutable( sp.sympify(expr) )

    a_list = list(miscprog.atoms(expr, sp.Number))
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
    miscprog = pycontroltools.auxfuncs.programming.miscprog

    expr = trunc_small_values(expr)

    maxden = int(1/eps)
    floats = miscprog.atoms(expr, sp.Float)
    rats = []
    for f in floats:
        rat = miscmath.sp_fff(f, maxden)
        rats.append(rat)

    return expr.subs(zip(floats, rats))

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

    a1 = exp1.atoms(sp.Symbol)
    a2 = exp2.atoms(sp.Symbol)

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
    
def dd(*args):
    """
    dd(a,b,c, ...) = np.dot(a, np.dot(b, np.dot(c, ...)))
    """
    return reduce(np.dot, args)
    
# copied from http://mail.scipy.org/pipermail/numpy-discussion/2008-February/031218.html

def matrixrank(A,tol=1e-8):
    s = np.linalg.svd(A,compute_uv=0)

    # if cond == True take 1, else 0:
    return np.sum(  np.where( s>tol, 1, 0 )  )



def zero_crossing_simulation(rhs, zcf, z0, t_values):
    """
    scipy.odeint does not provide a zero crossing function
    naive (and slow) approach

    rhs: rhs function
    zcf: the function whose zerocrossing shall be detected
         takes the state (shape =(n,m) returns shape=n
    z0: initial state
    t_values: time values (up to which the zc event is suspected)
    """

    res = scipy.integrate.odeint(rhs, z0, t_values)

    zerocc = zcf(res) # zero crossing candidate

    test_idx = 2 # ignore numerical noise at the beginning
    try:
        idx0 = np.where(np.sign(zerocc[test_idx:]) != np.sign(zerocc[test_idx]))[0][0]
    except IndexError:
        raise ValueError, "There was no zero crossing"

    idx0+= test_idx

    if zerocc[idx0] == 0:
        idx0+=1


    #IPS()

    t_values = t_values[:idx0]*1 # *1 to prevent referencing

    return t_values, res[:idx0, :]




def cont_continuation(x, stephight, threshold):
    """
    continuous continuation (for 1d-arrays)

    x           .... data


    stephight   ... the expected stephight (e.g 2*pi)

    threshold   .... smallest difference which is considered as a discontinuity
                    which has to be corrected
                    (must be greater than the Lipschitz-Const. of the signal
                     times dt)

    """
    x_d = np.concatenate(  ([0], np.diff(x))  )
    corrector_array = np.sign(x_d) * (np.abs(x_d) > threshold)
    corrector_array = np.cumsum(corrector_array) * -stephight

    return x + corrector_array



# copied from Mailinglist
# http://mail.scipy.org/pipermail/scipy-user/2008-October/018318.html
def extrema(x, max = True, min = True, strict = False, withend = False):
    """
    This function will index the extrema of a given array x.

    Options:
        max        If true, will index maxima
        min        If true, will index minima
        strict        If true, will not index changes to zero gradient
        withend    If true, always include x[0] and x[-1]

    This function will return a tuple of extrema indexies and values
    """

    # This is the gradient
    from numpy import zeros
    dx = zeros(len(x))
    from numpy import diff
    dx[1:] = diff(x)
    dx[0] = dx[1]

    # Clean up the gradient in order to pick out any change of sign
    from numpy import sign
    dx = sign(dx)

    # define the threshold for whether to pick out changes to zero gradient
    threshold = 0
    if strict:
        threshold = 1

    # Second order diff to pick out the spikes
    d2x = diff(dx)

    if max and min:
        d2x = abs(d2x)
    elif max:
        d2x = -d2x

    # Take care of the two ends
    if withend:
        d2x[0] = 2
        d2x[-1] = 2

    # Sift out the list of extremas
    from numpy import nonzero
    ind = nonzero(d2x > threshold)[0]

    return ind, x[ind]


##############################################################################

# control specific

##############################################################################

def pyc2d(a,b,Ts):

    """
    Algorithmus kopiert von Roberto Bucher

    Begründung: man erweitert den Zustand  xneu = (x,u)
    und sagt u_dot = 0 (weil u=konst.)
    Für das neue System bekommt man die zusammengestzte Matrix und pflückt
    sie hinterher wieder auseinander.

    """


    #a,b,c,d = SYS

    n=np.shape(a)[0]
    nb=np.shape(b)[1]
    #nc=np.shape(c)[0]

    ztmp=np.zeros((nb,n+nb))
    tmp=np.hstack((a,b))
    tmp=np.vstack((tmp,ztmp))

    tt = tmp

    tmp=sc.linalg.expm(tmp*Ts)

    A=tmp[0:n,0:n]
    B=tmp[0:n,n:n+nb]

    return A,B



def test_controlability(A, B):

    n,m = B.shape
    assert A.shape == (n,n)

    Ap = np.eye(n)
    elements = []
    for i in range(n):
        elements.append(np.dot(Ap, B))
        Ap = np.dot(Ap, A)

    C = np.column_stack(elements)

    return matrixrank(C) == n



def matrix_power(A, exp, cache):

    assert exp >= 0 and isinstance(exp, int)

    if not 0 in cache:
        cache[0] = np.eye(A.shape[0])
    if not 1 in cache:
        cache[1]=A

    assert A is cache[1]

    if exp in cache:
        return cache[exp]
    else:
        res = np.dot(A, matrix_power(A, exp-1, cache))
        cache[exp] = res
        return res


def mimo_controller_canonical_form_trafo(A, B):

    test_controlability(A, B)
    n,m = B.shape

    # columns of B:
    bi = list(B.T)

    # successive construction of L:
    # add new cols consisting of A**j*b[i] until linear dependence occurs
    # saving controllability indices

    rank = 0
    L = np.zeros((n, 0))
    cache_A = {}
    contr_idcs = []
    i = 0
    while 1:
        j = 0
        while 1:
            new_col = np.dot( matrix_power(A, j, cache_A), bi[i])
            L = np.column_stack((L, new_col))

            if matrixrank(L) == rank:
                # rank was not augmented
                L = L[:, :-1]
                contr_idcs.append(j)
                break

            # rank was augmented
            rank += 1

            j += 1
        i += 1
        if i == m or rank == n:
            break

    sigma_arr = np.cumsum(contr_idcs)

    assert L.shape == (n,n)
    assert sigma_arr[-1] == n

    iL = np.linalg.inv(L)
    #Tracer()()
    q_list = iL[sigma_arr-1, :] # -1 wegen 0-Indizierung

    Q = np.zeros((0, n))

    for i,d in enumerate(contr_idcs):
        for j in range(d):
            new_row = np.dot(q_list[i], matrix_power(A, j, cache_A))
            Q = np.row_stack((Q, new_row.reshape(-1, n)))


    assert Q.shape == (n,n)


    return L, Q

def null(A, eps=1e-10):
    """
    null-space of a Matrix or 2d-array
    """
    n, m = A.shape
    if n > m :
        return null(A.T, eps).T
        #return null(scipy.transpose(A), eps)
    u, s, vh = sc.linalg.svd(A)
    s=sc.append(s,sc.zeros(m))[0:m]
    null_mask = (s <= eps)
    null_space = sc.compress(null_mask, vh, axis=0)
    return null_space.T
