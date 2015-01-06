# -*- coding: utf-8 -*-
"""
.. module:: opt_polplatzierung
    :synopsis:  Functions to calculate a robust control matrix for multiple
                input systems.
    
.. moduleauthor:: Carsten Knoll
"""

import sympy as sp
import numpy as np

import pycontroltools.auxfuncs.math.numtools

def ortho_complement(M):
    r"""
    Gets a n,n-matrix M which is assumed to have rank n-1 and 
    returns a "column" v with :math:`v^T M = 0` and :math:`v^T v = 1`.
    
    **Parameters**
    
    * M : matrix of vectors with rank n-1
        sympy.Matrix
        
    **Return**
    
    * returns : orthogonal complement for columns of M
    * type : numpy.array (1d)
    """

    dtype = M.dtype

    M = sp.Matrix(M)
    n, m = M.shape
    assert n == m
    assert M.rank() == n-1

    v = M.T.nullspace()[0]
    v = np.array(np.array(v), dtype = dtype).squeeze()
    v/= np.linalg.norm(v)
    return v


def exchange_all_cols(V, P_list):
    """
    For every column in V:
    Calculates the 1-dimensional basis of the annihilator (:math:`:= a_j`)
    of all the other columns in V and projects :math:`a_j` to its
    correspondent space out of P_list.
    
    Then in V: replaces :math:`v_j` with the new normalized
    projected vector :math:`v_{j,projected}`.
    
    .. seealso:: :func:`opt_place_MI`
    
    **Parameters**
    
    * V : matrix of eigenvectors :math:`V = (v_1,..., v_n)`
        sympy.Matrix
    * P_list : list of spaces :math:`(S_{\mathrm{h}})_i` for :math:`i \in (1,...,n)`
        list
    
    **Return**
    
    * returns : new eigenvector matrix V
    * type : sympy.Matrix
    """

    n = V.shape[0]
    for i in range(n):
        v = V[:,i]
        V[:, i] = v*0

        v2 = ortho_complement(V)
        v2_projected = np.dot(P_list[i], v2)

        norm = np.linalg.norm(v2_projected)
        assert not norm == 0
        v2_projected /= norm

        V[:, i] = v2_projected

    return V

def full_qr(A, only_null_space=False):
    r"""
    Performs the QR numpy decomposition and augments the reduced
    orthonormal matrix :math:`Q_{red}` by its transposed null space
    :math:`\mathrm{null}(Q_{red})^T`
    (such that :math:`Q` is quadratic and regular).
    
        .. math::
                Q = \left(\begin{matrix} 
                        Q_{red} &
                        \mathrm{null}(Q_{red})^T
                    \end{matrix}\right)
    
    **Parameters**
    
    * A : matrix to be QR decomposed
        sympy.Matrix
        
    **Keyword Arguments**
    
    * only_null_space : only the null space of :math:`Q_{red}` will be returned
        boolean (default = False)
        
    **Return**
    
    * returns : Q (quadratic & regular)
    * type : sympy.Matrix
    * returns : r (upper triangular matrix)
    * type : sympy.Matrix
    """
    nt = pycontroltools.auxfuncs.math.numtools
    
    q, r = np.linalg.qr(A)
    n1, n2 = q.shape

    if n2 < n1:
        N = nt.null(q).T
    
    Q = np.hstack((q, N))
    
    n1, n2 = Q.shape
    assert n1 == n2

    if only_null_space:
        return N

    return Q, r

def opt_place_MI(A, B, *eigenvals, **kwargs):
    r"""
    Calculates and returns the optimal control matrix :math:`B_K` for the new
    system matrix :math:`(A + BB_K)` of the closed loop system by the
    algorithm described in [Reinschke14]_.
    
    **Parameters**
    
    * A : state matrix of the open loop
        sympy.Matrix    
    * B : input matrix
        sympy.Matrix    
    * eigenvals : desired eigenvalues for the closed loop system
        separate or as list-like object
        
    **Keyword Arguments**
    
    * rtol : relative tolerance of the change in the last iteration step of the
        resulting determinant of the eigenvector matrix
        
        real number (default = 0.01)
        
    **Return**
    
    * returns :  :math:`B_K`
    * type : sympy.Matrix
    """
    nt = pycontroltools.auxfuncs.math.numtools

    # default values
    rtol = kwargs.get('rtol', 0.01)    
    
    n1, n2 = A.shape
    assert n1 == n2
    n = n2
    
    n2, m = B.shape
    assert n1 == n
    
    assert len(eigenvals) > 0
    
    if hasattr(eigenvals[0], '__len__'):
        eigenvals = eigenvals[0]
    
    assert len(eigenvals) == n
    
    Lambda = np.diag(eigenvals)
    S_list = []
    P_list = [] # projector matrix
    
    # calculating valid subspaces for the correspondent eigenvector
    # (through QR decomposition):
    V0 = np.zeros((n,n))
    Qb, Rb = full_qr(B)
    Qb2 = full_qr(B, only_null_space=True) # "Uh"
    
    for i, s_i in enumerate(np.diag(Lambda)):
        Atmp = np.dot( (A - s_i*np.eye(n)).T, Qb2)
        Sh = full_qr(Atmp, only_null_space=True)
        S_list.append(Sh)

        V0[:,i] = Sh[:, 0] # construct first eigenvector matrix
        
        Q = Sh
        
        P_list.append(np.dot(Q, Q.T))
        
    # perform optimization:
    i = 0
    detV_new = 0
    while True:
        V = exchange_all_cols(V0, P_list)
        
        print i, np.linalg.det(V)
        i += 1
        
        detV_old = detV_new
        detV_new = np.linalg.det(V)
        
        # break if relative change for det(V) from loop before
        # to this loop is less than rtol
        if ((detV_new - detV_old)/detV_new < rtol):
            break
    
    # calculating the control matrix:
    Rinv = np.linalg.inv(Rb)
    Vinv = np.linalg.inv(V)
    Qb1 = Qb[:, :B.shape[1]] # "Uv"
    A_new = nt.dd(V, Lambda, Vinv)
    Bk_new = nt.dd(Rinv, Qb1.T, (A_new - A))
    
    return Bk_new