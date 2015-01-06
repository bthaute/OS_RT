# -*- coding: utf-8 -*-
"""
.. module:: lietools
    :synopsis: Functions concerning different types of Lie-derivatives.

.. moduleauthor:: Carsten Knoll, Chris Penndorf
"""

import sympy as sp

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

#chris:     Änderungen am Code:
#           -Check ob type(args) = sympy.Matrix entfernt
def lie_deriv(sf, vf, x, n = 1):
    """
    Calculates the Lie derivative of a scalar field :math:`\lambda(x)` along
    a vector field :math:`f(x)` (e.g. [Isidori]_):

        .. math::
            L_f\lambda(x)   = \\frac{ \partial{\lambda(x)} }{ \partial{x} }f(x)
                            = grad^{T}\lambda(x)\cdot f(x)

    with :math:`L_f^n\lambda(x)=\\frac{ \partial{L_f^{n-1}\lambda(x)} }
    { \partial{x} }f(x)` and :math:`L_f^{0}\lambda(x) := \lambda(x)`

    **Parameters**

    * sf : scalar field to be derived
        function
    * vf : vector field to derive along
        vector
    * x : coordinates for derivation
        list
    * n : number of derivations
        non-negative integer

    **Return**

    * returns : scalar field
    * type : function

    **Examples**

        >>> import sympy

        >>> x1,x2,x3 = sympy.symbols('x1 x2 x3')
        >>> h = x1**2 + 2*x2 + x3
        >>> f = sympy.Matrix([x1*x2, x3, x2])
        >>> x = [x1, x2, x3]

        >>> lie_deriv(h, f, x, n=1)
        2*x1**2*x2 + x2 + 2*x3

    .. seealso:: :func:`lie_bracket`, :func:`lie_deriv_covf`

    """
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
    Calculates the Lie bracket for the vector field :math:`g(x)`
    along the vector field :math:`f(x)` (e.g. [Isidori]_):

        .. math::
            [f,g] = \\frac{ \partial{g(x)} }{ \partial{x} }f(x)
                    - \\frac{ \partial{f(x)} }{ \partial{x} }g(x)
                  = ad_fg(x)

    with :math:`\\quad`:math:`ad_f^ng(x) = [f, ad_f^{n-1}g](x)` :math:`\\quad`
    and :math:`\\quad` :math:`ad_f^0g(x) := g(x)`

    **Parameters**

    * f : vector field (direction for derivation)
        Matrix (shape: (n, 1)) / iterable
    * g : vector field to be derived
        Matrix (shape: (n, 1)) / iterable
    * args : coordinates
        separate scalar symbols or as iterable

    **Keyword Arguments**

    * n : number of derivations
        non-negative integer (default = 1)

    **Exceptions**

    * AssertionError : non-matching shapes of f, g, args

    **Return**

    * returns : vector field
    * type : Matrix

    **Examples**

        >>> import sympy

        >>> x1,x2,x3 = sympy.symbols('x1 x2 x3')
        >>> g = [2*x2, x1**2, 2*x3]
        >>> f = [x1*x2, x3, x2]

        >>> lie_bracket(g, f, x1, x2, x3, n=1)
        Matrix([
        [x1**3 + 2*x2**2 - 2*x3],
        [    -2*x1**2*x2 + 2*x3],
        [          x1**2 - 2*x2]])

    .. seealso:: :func:`lie_deriv`, :func:`lie_deriv_covf`

    """

    assert len(args) > 0

    if hasattr(args[0], '__len__'):
        args = args[0]

    n = kwargs.get('n', 1)  # if not given, then n = 1

    if n == 0:
        return g

    assert n > 0 # and isinstance(n, int)
    assert len(args) == len(list(f))

    # convert in sympy matrices
    f = sp.Matrix(f)
    g = sp.Matrix(g)

    assert f.shape == g.shape
    assert f.shape[1] == 1

    jf = f.jacobian(args)
    jg = g.jacobian(args)

    res = jg * f - jf * g

    if n > 1:
        res = lie_bracket(f, res, *args, n=n-1)

    return res


def lie_deriv_covf(w, f, *args, **kwargs):
    """
    Calculates the Lie derivative of the covector field :math:`\omega(x)`
    along the vector field :math:`f(x)` (e.g. [Isidori]_):

        .. math::
            L_f\omega(x) = f^T(x) \\left( \\frac{ \partial{\omega^T(x)} }
                    { \partial{x} } \\right)^{T} + \omega(x)
                    \\frac{ \partial{f(x)} }
                    { \partial{x} }


    with

       .. math::
            L_f^n\omega(x) = f^T(x) \\left( \\frac{ \partial{(L_f^{n-1}\omega)^T(x)} }
                    { \partial{x} } \\right)^{T} + (L_f^{n-1}\omega)(x)
                    \\frac{ \partial{f(x)} }
                    { \partial{x} }

    and :math:`\\quad` :math:`L_f^0\omega(x) := \omega(x)`

    Includes the option to omit the transposition of :math:`\\quad`
    :math:`\\frac{ \partial{\omega^T(x)} }{ \partial{x} }` :math:`\\quad`
    with :code:`transpose_jac = False`:

        .. math::
            L_f\omega(x) = f^T(x) \\left( \\frac{ \partial{\omega^T(x)} }
                    { \partial{x} } \\right) + \omega(x)\\frac{ \partial{f(x)} }
                    { \partial{x} }


    **Parameters**

    * w : covector field to be derived
        Matrix of shape (1,n)
    * f : vector field (direction of derivation)
        Matrix of shape (n,1)
    * args : coordinates
        separate or as list-like object

    **Keyword Arguments**

    * n : number of derivations
        non-negative integer (default = 1)
    * transpose_jac : transposition of :math:`\\frac{ \partial{\omega^T(x)} }{ \partial{x} }`
        boolean (default = True)(Background: needed for some special applications)

    **Exceptions**

    * AssertionError : non-matching shapes of w, f, args

    **Return**

    * returns : covector field
    * type : sympy.Matrix

    **Examples**

        >>> import sympy

        >>> x1,x2,x3 = sympy.symbols('x1 x2 x3')
        >>> w = sympy.Matrix([[2*x2, x1**2, 2*x3]])
        >>> f = sympy.Matrix([x1, x2, x3])

        >>> lie_deriv_covf(w, f, x1, x2, x3, n=1)
        Matrix([[4*x2, 3*x1**2, 4*x3]])

    .. seealso:: :func:`lie_deriv`, :func:`lie_bracket`
    """
    assert isinstance(w, sp.Matrix) and isinstance(f, sp.Matrix)
    k,l = w.shape
    m, n = f.shape
    assert  k==1 and n==1
    assert l==m

    if hasattr(args[0], '__len__'):
        args = args[0]

    assert len(args) == len(list(f))

    n = kwargs.get('n', 1)  # if not given, then n = 1

    if n == 0:
        return w

    assert n > 0  # and isinstance(n, int)

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

    # TODO: pass all keyword args, not only n
    if n > 1:
        res = lie_deriv_covf(res, f, args, n = n-1)

    return res