# -*- coding: utf-8 -*-
"""
Created on Wed Feb 12 10:59:29 2014

.. module:: test_robust_poleplacement
    :synopsis: Unittest for the module robust_poleplacement
    
.. moduleauthor:: Chris Penndorf
"""

import sympy as sp
import numpy as np

import robust_poleplacement as rp
import pycontroltools.auxfuncs.math.numtools

import unittest

class testBadInput(unittest.TestCase):
    #examples    
    #--------
    vec0 = sp.Matrix([0,0,0])
    vec1 = sp.Matrix([1,4,7])
    vec2 = sp.Matrix([12,35,28])
    vec3 = sp.Matrix([23,16,9])
    M_wrongShape1 = np.hstack([vec1,vec2])
    M_wrongShape2 = np.hstack([vec1,vec2,vec3,vec1])
    M_wrongRank1 = np.hstack([vec1,vec2,vec3])
    M_wrongRank2 = np.hstack([vec1,vec0,vec0])
                   
    def testWrongTypeOrtho_complement(self):
        """ ortho_complement should fail with wrong type/shape of input"""
        self.assertRaises(AssertionError, rp.ortho_complement,
                          self.M_wrongShape1)        
        self.assertRaises(AssertionError, rp.ortho_complement,
                          self.M_wrongShape2)
                          
    def testWrongRankOrtho_complement(self):
        """ ortho_complement should fail with rank(M) != n-1"""
        self.assertRaises(AssertionError, rp.ortho_complement,
                          self.M_wrongRank1)        
        self.assertRaises(AssertionError, rp.ortho_complement,
                          self.M_wrongRank2)        
    
class KnownValues(unittest.TestCase):
    nt = pycontroltools.auxfuncs.math.numtools
    
    # examples
    #===========
    vec0 = sp.Matrix([0,0,0])
    
    # vectors for M1
    vec100 = sp.Matrix([1,0,0])
    vec010 = sp.Matrix([0,1,0])
    M1 = nt.to_np(np.hstack([vec100, vec010, vec0]))
    v1 = np.array([0,0,1])
    
    # vectors for M2
    vec1 = sp.Matrix([1,4,7])
    vec2 = sp.Matrix([12,35,28])
    M2 = nt.to_np(np.hstack([vec1, vec2, vec0]))
        
    def testKnownValuesOrtho_complement(self):
        """ ortho_complement should give known result with known input"""
        v1_py = rp.ortho_complement(self.M1)
        self.assertTrue(all(v1_py == self.v1))
        
    def testOrthonormalityOrtho_complement(self):
        """ result of ortho_complement should fullfill orthonormality:
            v.T*M = 0  and  v.T*v = 1
        """
        nt = pycontroltools.auxfuncs.math.numtools
        
        v2 = rp.ortho_complement(self.M2)
        self.assertTrue(np.allclose(np.dot(v2.T,self.M2),
                        nt.to_np(self.vec0), rtol=0, atol=1e-08))
        self.assertTrue(np.allclose(np.dot(v2.T,v2),
                        1, rtol=0, atol=1e-08))
    
if __name__ == "__main__":
	unittest.main()