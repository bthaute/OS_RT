# -*- coding: utf-8 -*-
"""
Created on Mon Nov 24 15:18:13 2014

@author: herrmann
"""

import sympy as sp

q11, q12, q21, q22 = sp.symbols(["q11", "q12", "q21", "q22"])
#q11_d, q12_d,q21_d, q22_d =sp.symbols("q11_d, q12_d,q21_d, q22_d")
parameter = sp.symbols("I11, I12, I21, I22, m11, m12, m21, m22, l11, l12, l21, l22, a11, a12, a21, a21, k1, k2, d1, d2, g")
I11, I12, I21, I22, m11, m12, m21, m22, l11, l12, l21, l22, a11, a12, a21, a21, k1, k2, d1, d2, g = parameter

A=sp.Matrix([[I11 + I12 + I21 + I22 + m11*(2*l11**2*sp.sin(q11)**2 + 2*l11**2*sp.cos(q11)**2)/2 + m12*((-2*a11*sp.sin(q11) - 2*l12*sp.sin(q11 + q12))*(-a11*sp.sin(q11) - l12*sp.sin(q11 + q12)) + (a11*sp.cos(q11) + l12*sp.cos(q11 + q12))*(2*a11*sp.cos(q11) + 2*l12*sp.cos(q11 + q12)))/2 + m21*((-2*a11*sp.sin(q11) - 2*a12*sp.sin(q11 + q12) - 2*l21*sp.sin(q11 + 2*q12))*(-a11*sp.sin(q11) - a12*sp.sin(q11 + q12) - l21*sp.sin(q11 + 2*q12)) + (a11*sp.cos(q11) + a12*sp.cos(q11 + q12) + l21*sp.cos(q11 + 2*q12))*(2*a11*sp.cos(q11) + 2*a12*sp.cos(q11 + q12) + 2*l21*sp.cos(q11 + 2*q12)))/2 + m22*((-2*a11*sp.sin(q11) - 2*a12*sp.sin(q11 + q12) - 2*a21*sp.sin(q11 + 2*q12) - 2*l22*sp.sin(q11 + 2*q12 + q22))*(-a11*sp.sin(q11) - a12*sp.sin(q11 + q12) - a21*sp.sin(q11 + 2*q12) - l22*sp.sin(q11 + 2*q12 + q22)) + (a11*sp.cos(q11) + a12*sp.cos(q11 + q12) + a21*sp.cos(q11 + 2*q12) + l22*sp.cos(q11 + 2*q12 + q22))*(2*a11*sp.cos(q11) + 2*a12*sp.cos(q11 + q12) + 2*a21*sp.cos(q11 + 2*q12) + 2*l22*sp.cos(q11 + 2*q12 + q22)))/2, I12 + 2*I21 + 2*I22 + m12*(-l12*(-2*a11*sp.sin(q11) - 2*l12*sp.sin(q11 + q12))*sp.sin(q11 + q12) + l12*(2*a11*sp.cos(q11) + 2*l12*sp.cos(q11 + q12))*sp.cos(q11 + q12))/2 + m21*((-a12*sp.sin(q11 + q12) - 2*l21*sp.sin(q11 + 2*q12))*(-2*a11*sp.sin(q11) - 2*a12*sp.sin(q11 + q12) - 2*l21*sp.sin(q11 + 2*q12)) + (a12*sp.cos(q11 + q12) + 2*l21*sp.cos(q11 + 2*q12))*(2*a11*sp.cos(q11) + 2*a12*sp.cos(q11 + q12) + 2*l21*sp.cos(q11 + 2*q12)))/2 + m22*((-a12*sp.sin(q11 + q12) - 2*a21*sp.sin(q11 + 2*q12) - 2*l22*sp.sin(q11 + 2*q12 + q22))*(-2*a11*sp.sin(q11) - 2*a12*sp.sin(q11 + q12) - 2*a21*sp.sin(q11 + 2*q12) - 2*l22*sp.sin(q11 + 2*q12 + q22)) + (a12*sp.cos(q11 + q12) + 2*a21*sp.cos(q11 + 2*q12) + 2*l22*sp.cos(q11 + 2*q12 + q22))*(2*a11*sp.cos(q11) + 2*a12*sp.cos(q11 + q12) + 2*a21*sp.cos(q11 + 2*q12) + 2*l22*sp.cos(q11 + 2*q12 + q22)))/2],
[                                                                                                                                                                I12 + 2*I21 + 2*I22 + m12*(-2*l12*(-a11*sp.sin(q11) - l12*sp.sin(q11 + q12))*sp.sin(q11 + q12) + 2*l12*(a11*sp.cos(q11) + l12*sp.cos(q11 + q12))*sp.cos(q11 + q12))/2 + m21*((-2*a12*sp.sin(q11 + q12) - 4*l21*sp.sin(q11 + 2*q12))*(-a11*sp.sin(q11) - a12*sp.sin(q11 + q12) - l21*sp.sin(q11 + 2*q12)) + (2*a12*sp.cos(q11 + q12) + 4*l21*sp.cos(q11 + 2*q12))*(a11*sp.cos(q11) + a12*sp.cos(q11 + q12) + l21*sp.cos(q11 + 2*q12)))/2 + m22*((-2*a12*sp.sin(q11 + q12) - 4*a21*sp.sin(q11 + 2*q12) - 4*l22*sp.sin(q11 + 2*q12 + q22))*(-a11*sp.sin(q11) - a12*sp.sin(q11 + q12) - a21*sp.sin(q11 + 2*q12) - l22*sp.sin(q11 + 2*q12 + q22)) + (2*a12*sp.cos(q11 + q12) + 4*a21*sp.cos(q11 + 2*q12) + 4*l22*sp.cos(q11 + 2*q12 + q22))*(a11*sp.cos(q11) + a12*sp.cos(q11 + q12) + a21*sp.cos(q11 + 2*q12) + l22*sp.cos(q11 + 2*q12 + q22)))/2,                                                                                                                                       I12 + 4*I21 + 4*I22 + m12*(2*l12**2*sp.sin(q11 + q12)**2 + 2*l12**2*sp.cos(q11 + q12)**2)/2 + m21*((-2*a12*sp.sin(q11 + q12) - 4*l21*sp.sin(q11 + 2*q12))*(-a12*sp.sin(q11 + q12) - 2*l21*sp.sin(q11 + 2*q12)) + (a12*sp.cos(q11 + q12) + 2*l21*sp.cos(q11 + 2*q12))*(2*a12*sp.cos(q11 + q12) + 4*l21*sp.cos(q11 + 2*q12)))/2 + m22*((-2*a12*sp.sin(q11 + q12) - 4*a21*sp.sin(q11 + 2*q12) - 4*l22*sp.sin(q11 + 2*q12 + q22))*(-a12*sp.sin(q11 + q12) - 2*a21*sp.sin(q11 + 2*q12) - 2*l22*sp.sin(q11 + 2*q12 + q22)) + (a12*sp.cos(q11 + q12) + 2*a21*sp.cos(q11 + 2*q12) + 2*l22*sp.cos(q11 + 2*q12 + q22))*(2*a12*sp.cos(q11 + q12) + 4*a21*sp.cos(q11 + 2*q12) + 4*l22*sp.cos(q11 + 2*q12 + q22)))/2]])
print "..."
n=A.shape[0]
L=sp.zeros(n)

for i in range(0,n):
    for j in range(0,i+1):
        Summe = A[i,j]
        for k in range(0,j+1):
            Summe = Summe - A[i,k] * A[j,k]
            if i > j:
                L[i,j] = Summe / L[j,j]
            else:
                L[i,i] = sp.sqrt(Summe)
            