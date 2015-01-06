# -*- coding: utf-8 -*-
"""
.. module:: miscprog
    :synopsis: Miscellaneous helpfunctions concerning programming issues.
    
.. moduleauthor:: Carsten Knoll
"""
#__name__= "auxfuncs.programming.miscprog"
#from __future__ import absolute_import

import sympy as sp
import numpy as np

import itertools as it

import pycontroltools.auxfuncs.math.matrix
import pycontroltools.auxfuncs.math.numtools
    
class Container(object):
    def __init__(self, **kwargs):
        #chris: Assertion Error, falls eines der Keywords bereits vorhanden ist
        assert len( set(dir(self)).intersection(kwargs.keys()) ) == 0
        self.__dict__.update(kwargs)
        
def make_global(varList):
    """
    injects the symbolic variables of a collection to the global namespace
    usefull for interactive sessions
    """
    if not isinstance(varList, (list, tuple)):
        if isinstance(varList, sp.Matrix):
            varList = np.array(varList).flatten()
        else:
            raise TypeError, 'Unexpected type for varList'

    import inspect

    # get the topmost frame
    frame = inspect.currentframe()
    while True:
        if frame.f_back == None:
            break
        frame = frame.f_back


    # this is strongly inspired by sympy.var
    try:
        for v in varList:
            if v.is_Function:
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

#chris: kürzt die sympy.printing.preview-Funktion ab
#       indem es den Output direkt auf 'pdf' setzt
def prev(expr, **kwargs):
    """
    sympy preview abbreviation
    """
    KWargs = {'output':'pdf'}
    KWargs.update(kwargs)
    sp.preview(expr, **KWargs)
    
# TODO: seems to conflict with zip0
def tup0(xx):
    """
    helper function for substituting.
    takes (x1, x2, x3, ...)
    returns [(x1, 0), (x2, 0), ...]
    """

    return zip(xx, [0]*len(xx))
    
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
            
def rev_tuple(tup):
    """
    """
    return [(t[1], t[0]) for t in tup]
    
def trig_term_poly(expr, s):
    """
    s ... the argument of sin, cos
    """

    X, Y = sp.symbols('tmpX_, tmp_Y')

    poly = sp.Poly( expr.subs([(sp.sin(s), X), (sp.cos(s),Y)]), X,Y, domain='EX')

    return poly

def atoms(expr, *args, **kwargs):
    """
    """
    matrix = pycontroltools.auxfuncs.math.matrix
    
    if isinstance(expr, (sp.Matrix, list)):
        return matrix.matrix_atoms(expr, *args, **kwargs)
    else:
        return expr.atoms(*args, **kwargs)
        
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
    
def zip0(xx, arg = 0):
    """ handy for subtituting equilibrium points"""

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

    """
    matrix = pycontroltools.auxfuncs.math.matrix
    numtools = pycontroltools.auxfuncs.math.numtools

    # TODO: sympy-Matrizen mit Stückweise definierten Polynomen
    # numpy fähig (d.h. vektoriell) auswerten


    # TODO: Unittest


    # TODO: only relevant if numpy is in modules

    expr = sp.sympify(expr)
    expr = matrix.ensure_mutable(expr)
    expr_tup = aux_make_tup_if_necc(expr)
    arg_tup = aux_make_tup_if_necc(args)

    new_expr = []
    arg_set = set(arg_tup)
    for e in expr_tup:
        assert isinstance(e, sp.Expr)
        # args (Symbols) which are not in that expression
        diff_set = arg_set.difference(e.atoms(sp.Symbol))

        for d in diff_set:
            assert isinstance(d, sp.Symbol)
            e = sp.Add(e, d, -d, evaluate = False)

        new_expr.append(e)

    if not hasattr(expr, '__len__'):
        assert len(new_expr) == 1
        new_expr = new_expr[0]



    printer = kwargs.get('printer', None)
    use_imps = kwargs.get('use_imps', True)
    func = sp.lambdify(args, new_expr, modules, printer, use_imps)




    if kwargs.get('np_vectorize', False):
        func1 = np.vectorize(func)
    else:
        func1 = func

    if kwargs.get('special_vectorize', False):
        def func2(*allargs):
            return numtools.to_np(func(*allargs))

        f = np.float
        func3 = np.vectorize(func2, otypes = [f,f,f, f,f,f])
        return func3

    if kwargs.get('np_wrapper', False):
        def func2(*allargs):
            return numtools.to_np(func1(*allargs))
    elif kwargs.get('list_wrapper', False):
        def func2(*allargs):
            return list(func1(*allargs))
    else:
        func2 = func1
    return func2

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
    
def simp_trig_dict(sdict):
    """
    takes a sorted dict, simplifies each value and adds all up
    """

    items = sdict.items()

    res = 0
    for k, v in items:
        res += k*sp.trigsimp(sum(v))

    return res