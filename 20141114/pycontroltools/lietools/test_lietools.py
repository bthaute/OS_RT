# -*- coding: utf-8 -*-
"""
Created on Wed Feb 12 10:59:29 2014

.. module:: test_lietools
    :synopsis: Unittest for the module lietools

.. moduleauthor:: Chris Penndorf
"""

import lietools as lie
import sympy as sp
import unittest

# notice, that the variable "w"
# is the same as the "omega" for the covector field in the equation

x1,x2,x3 = sp.symbols('x1 x2 x3')

class testBadInput(unittest.TestCase):
    #examples
    #--------
    # args
    x = [x1, x2, x3]
    x_wrong = [x1+x2, x2]
    # scalar field
    h = x1**2 + 2*x2 + x3
    # vector fields f & g
    f = sp.Matrix([x1*x2, x3, x2])
    g = sp.Matrix([2*x2, x1**2, 2*x3])
    # covector field w
    w = sp.Matrix([[2*x2, x1**2, 2*x3]])
    w_wrong = sp.Matrix([2*x2, x1**2, 2*x3])
    # matrix
    M = sp.Matrix([x1*x2, x3, x2,
                   2*x2, x1**2, 2*x3,
                   2*x2, x1**2, 2*x3]).reshape(3,3)

    def testWrongTypeJac(self):
        """ jac should fail with wrong type/shape of input"""
        self.assertRaises(TypeError, lie.jac, self.M, self.x)
        self.assertRaises(ValueError, lie.jac, self.h, self.x_wrong)

    def testWrongTypeLie_Deriv(self):
        """ lie_deriv should fail with wrong type/shape of input"""
        self.assertRaises(AssertionError, lie.lie_deriv,
                          self.f, self.f, self.x)
        self.assertRaises(TypeError, lie.lie_deriv,
                          self.M, self.f, self.x)
        self.assertRaises(AssertionError, lie.lie_deriv,
                          self.h, self.M, self.x)

    def testWrongTypeLie_Bracket(self):
        """ lie_bracket should fail with wrong type/shape of input"""
        self.assertRaises(AssertionError, lie.lie_bracket,
                          self.M, self.f, self.x)
        self.assertRaises(AssertionError, lie.lie_bracket,
                          self.f.T, self.g, self.x)
        self.assertRaises(AssertionError, lie.lie_bracket,
                          self.f, self.M, self.x)

    def testWrongTypeLie_Deriv_Covf(self):
        """ lie_deriv_covf should fail with wrong type/shape of input"""
        self.assertRaises(AssertionError, lie.lie_deriv_covf, self.M,
                          self.f, self.x)
        self.assertRaises(AssertionError, lie.lie_deriv_covf, self.w,
                          self.M, self.x)
        self.assertRaises(AssertionError, lie.lie_deriv_covf, self.w_wrong,
                          self.f, self.x)



class KnownValues(unittest.TestCase):
    # examples
    #===========
    # args
    x = [x1, x2, x3]
    # scalar field
    h = x1**2 + 2*x2 + x3
    # scalar fields (harmonic & exponential)
    h_harm = sp.sin(x1) + sp.cos(x2) + sp.acosh(x3)
    h_exp = sp.exp(-x1) + sp.exp(2*x2) + sp.exp(3*x3)
    # vector fields f & g
    f = sp.Matrix([x1*x2, x3, x2])
    g = sp.Matrix([2*x2, x1**2, 2*x3])
    # vector field f (harmonic functions)
    f_harm = sp.Matrix([sp.sin(x1), sp.cos(x2), sp.acosh(x3)])
    f_exp = sp.Matrix([sp.exp(-x1), sp.exp(2*x2), sp.exp(3*x3) ])

    fields = (h, h_harm, h_exp, f, f_harm, f_exp, g)
    # covector field w
    w = sp.Matrix([[2*x2, x1**2, 2*x3]])
    # covector field w (harmonic functions)
    w_harm = sp.Matrix([[sp.sin(x1), sp.cos(x2), sp.acosh(x3)]])
    w_exp = sp.Matrix([[sp.exp(-x1), sp.exp(2*x2), sp.exp(3*x3)]])


    # results calculated in MATLAB R2013a (8.1.0.604)
    # using MATLAB function scripts:
    # lie_deriv.m
    # lie_bracket.m
    # lie_deriv_covf.m
    #-----------------------------------------------

    # for funtion jac()
    h_jac = sp.Matrix([[2*x1, 2, 1]])
    h_harm_jac = sp.Matrix([[sp.cos(x1),
                            -sp.sin(x2),1/sp.sqrt(x3**2 - 1)]])
    h_exp_jac = sp.Matrix([[-sp.exp(-x1), 2*sp.exp(2*x2), 3*sp.exp(3*x3)]])

    f_jac = sp.Matrix([x2,x1,0,0,0,1,0,1,0]).reshape(3,3)
    f_harm_jac = sp.Matrix([sp.cos(x1),0,0,
                            0,-sp.sin(x2),0,
                            0,0,1/sp.sqrt(x3**2 - 1)]).reshape(3,3)
    f_exp_jac = sp.Matrix([-sp.exp(-x1),0,0,
                           0, 2*sp.exp(2*x2),0,
                           0,0, 3*sp.exp(3*x3)]).reshape(3,3)
    g_jac = sp.Matrix([0,2,0,2*x1,0,0,0,0,2]).reshape(3,3)

    jac_results = (h_jac, h_harm_jac, h_exp_jac,
                   f_jac, f_harm_jac, f_exp_jac,
                   g_jac)

    # for function lie_deriv()
    hf_deriv = 2*x2*x1**2 + x2 + 2*x3
    hf_deriv2 = 2*x2 + 4*x1**2*x2**2 + x3*(2*x1**2 + 1)
    hg_deriv = 2*x1**2 + 4*x2*x1 + 2*x3
    hg_deriv2 = 4*x3 + 4*x1**3 + 2*x2*(4*x1 + 4*x2)
    h_harm_f_deriv = x2/sp.sqrt(x3**2 - 1) - x3*sp.sin(x2) + x1*x2*sp.cos(x1)
    h_exp_f_deriv = 2*x3*sp.exp(2*x2) + 3*x2*sp.exp(3*x3) - x1*x2*sp.exp(-x1)
    deriv_results = [(f, hf_deriv, hf_deriv2), (g, hg_deriv, hg_deriv2)]

    # for function lie_bracket()
    fg_bracket = sp.Matrix([ -x1**3 - 2*x2**2 + 2*x3,
                            2*x2*x1**2 - 2*x3,
                            -x1**2 + 2*x2])
    # symmetry [f,g] = - [g,f]
    gf_bracket = sp.Matrix([x1**3 + 2*x2**2 - 2*x3,
                            -2*x2*x1**2 + 2*x3,
                            x1**2 - 2*x2])
    # second order bracket
    ffg_bracket = sp.Matrix([ 2*x2 + x2*(x1**3 + 2*x2**2 - 2*x3) - 4*x2*x3
                                - 3*x1**3*x2 + x1*(- 2*x2*x1**2 + 2*x3),
                             4*x1**2*x2**2 - 4*x2 + 2*x1**2*x3 + x1**2,
                             - 4*x2*x1**2 + 4*x3])
    # with harmonic functions
    f_harm_g_bracket = sp.Matrix([2*sp.cos(x2) - 2*x2*sp.cos(x1),
                                  x1**2*sp.sin(x2) + 2*x1*sp.sin(x1),
                                  2*sp.acosh(x3) - (2*x3)/sp.sqrt(x3**2 - 1)])
    # with exponential functions
    f_exp_g_bracket = sp.Matrix([2*sp.exp(2*x2) + 2*x2*sp.exp(-x1),
                                 2*x1*sp.exp(-x1) - 2*x1**2*sp.exp(2*x2),
                                 2*sp.exp(3*x3) - 6*x3*sp.exp(3*x3)])

    # for function lie_deriv_covf()
    wf_deriv_covf = sp.Matrix([[2*x2**2 + 2*x3,
                                2*x2*x1**2 + 2*x2*x1 + 2*x3,
                                x1**2 + 2*x2]])
    wf_deriv_covf_notrans = sp.Matrix([[2*x2**2 + 2*x1*x3,
                                        2*x3 + 4*x1*x2,
                                        x1**2 + 2*x2]])
    wg_deriv_covf = sp.Matrix([[2*x1**3 + 2*x1**2,
                                4*x2 + 4*x1*x2,
                                8*x3]])
    wf_deriv_covf2 = sp.Matrix([[2*x2 + x2*(2*x2**2 + 2*x3) + 4*x2*x3,
                                 4*x2 + x3*(2*x1**2 + 2*x1) +
                                 x1*(2*x2**2 + 2*x3) + x1**2 +
                                 x1*x2*(2*x2 + 4*x1*x2),
                                 4*x2*x1**2 + 2*x2*x1 + 4*x3]])
    w_harm_g_deriv_covf = sp.Matrix([[2*x1*sp.cos(x2) + 2*x2*sp.cos(x1),
                                      2*sp.sin(x1) - x1**2*sp.sin(x2),
                                      2*sp.acosh(x3)+(2*x3)/sp.sqrt(x3**2-1)]])
    w_exp_g_deriv_covf = sp.Matrix([[2*x1*sp.exp(2*x2) - 2*x2*sp.exp(-x1),
                                     2*sp.exp(-x1) + 2*x1**2*sp.exp(2*x2),
                                     2*sp.exp(3*x3) + 6*x3*sp.exp(3*x3)]])

    def testKnownValuesJac(self):
        """ jac should give known result with known input"""
        for field, result in zip(self.fields, self.jac_results):
            py_result = lie.jac(field, self.x)
            self.assertEqual(result, py_result)

    def testKnownValuesLie_Deriv(self):
        """ lie_deriv should give known result with known input
        (also for higher order derivatives)"""
        for field, result1, result2 in self.deriv_results:

            py_result1 = lie.lie_deriv(self.h, field, self.x)
            self.assertEqual(result1, py_result1)

            py_result2 = lie.lie_deriv(self.h, field, self.x, n=2)
            self.assertEqual(result2, py_result2) # second order

        py_result_harm = lie.lie_deriv(self.h_harm, self.f, self.x)
        self.assertEqual(self.h_harm_f_deriv, py_result_harm)

        py_result_exp = lie.lie_deriv(self.h_exp, self.f, self.x)
        self.assertEqual(self.h_exp_f_deriv, py_result_exp)

    def testKnownValuesLie_Bracket(self):
        """lie_bracket should give known result with known input
        (also for higher order derivatives)
        and check for symmetry [f,g] = -[g,f]
        """
        py_result_fg = lie.lie_bracket(self.f, self.g, self.x)
        self.assertEqual(self.fg_bracket, py_result_fg)

        py_result_gf = lie.lie_bracket(self.g, self.f, self.x)
        self.assertEqual(self.gf_bracket, py_result_gf)

        py_result_ffg = lie.lie_bracket(self.f, self.g, self.x, n=2)
        self.assertEqual(sp.expand(self.ffg_bracket),
                         sp.expand(py_result_ffg)) # second order

        self.assertEqual(py_result_fg, -1*py_result_gf) # test for symmetry

        py_result_f_harm_g = lie.lie_bracket(self.f_harm, self.g, self.x)
        self.assertEqual(self.f_harm_g_bracket, py_result_f_harm_g)

        py_result_f_exp_g = lie.lie_bracket(self.f_exp, self.g, self.x)
        self.assertEqual(self.f_exp_g_bracket, py_result_f_exp_g)

    def testKnownValuesLie_Deriv_Covf(self):
        """lie_deriv_covf should give known result with known input
        and check for functionality of kwarg: transpose_jac"""
        py_result_wf = lie.lie_deriv_covf(self.w, self.f, self.x)
        self.assertEqual(self.wf_deriv_covf, py_result_wf)

        py_result_wf_notrans = lie.lie_deriv_covf(self.w, self.f, self.x,
                                                  transpose_jac = False)
        self.assertEqual(self.wf_deriv_covf_notrans, py_result_wf_notrans)

        py_result_wg = lie.lie_deriv_covf(self.w, self.g, self.x)
        self.assertEqual(self.wg_deriv_covf, py_result_wg)

        py_result_wf2 = lie.lie_deriv_covf(self.w, self.f, self.x,
                                           n = 2)
        self.assertEqual(self.wf_deriv_covf2, py_result_wf2) # second order

        py_result_w_harm_g = lie.lie_deriv_covf(self.w_harm, self.g, self.x)
        self.assertEqual(self.w_harm_g_deriv_covf, py_result_w_harm_g)

        py_result_w_exp_g = lie.lie_deriv_covf(self.w_exp, self.g, self.x)
        self.assertEqual(self.w_exp_g_deriv_covf, py_result_w_exp_g)

if __name__ == "__main__":
	unittest.main()