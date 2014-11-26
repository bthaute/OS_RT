# -*- coding: utf-8 -*-
"""
Created on Mon Nov 24 15:18:13 2014

@author: herrmann
"""

import sympy as sp
import pickle
import numpy as np
from parameter_springs import para_g, para_m, para_l, para_a, para_k, para_d, para_I
import symb_tools as st
from IPython import embed as IPS
#
#
#q11, q12, q21, q22 = sp.symbols(["q11", "q12", "q21", "q22"])
##q11_d, q12_d,q21_d, q22_d =sp.symbols("q11_d, q12_d,q21_d, q22_d")
#parameter = sp.symbols("I11, I12, I21, I22, m11, m12, m21, m22, l11, l12, l21, l22, a11, a12, a21, a21, k1, k2, d1, d2, g")
#I11, I12, I21, I22, m11, m12, m21, m22, l11, l12, l21, l22, a11, a12, a21, a21, k1, k2, d1, d2, g = parameter
#params_values = {"m11":para_m[0,0], "m12":para_m[0,1], "m21":para_m[1,0], "m22":para_m[1,1], "I11":para_I[0,0] ,"I12":para_I[0,1], "I21":para_I[1,0] ,"I22":para_I[1,1], "a11":para_a[0,0], "a12":para_a[0,1], "a21":para_a[1,0], "a22":para_a[1,1],
#                 "l11": para_l[0,0], "l12": para_l[0,1], "l21": para_l[1,0], "l22": para_l[1,1], "k1":para_k[0], "k2":para_k[1], "d1":para_d[0], "d2":para_d[1], "g":para_g }
#para_q={"q11":0.1,"q12":0.2,"q21":0.3,"q22":0.4,"q11_d":0.5,"q12_d":0.6,"q21_d":0.7,"q22_d":0.8}
##A=sp.Matrix([[I11 + I12 + I21 + I22 + m11*(2*l11**2*sp.sin(q11)**2 + 2*l11**2*sp.cos(q11)**2)/2 + m12*((-2*a11*sp.sin(q11) - 2*l12*sp.sin(q11 + q12))*(-a11*sp.sin(q11) - l12*sp.sin(q11 + q12)) + (a11*sp.cos(q11) + l12*sp.cos(q11 + q12))*(2*a11*sp.cos(q11) + 2*l12*sp.cos(q11 + q12)))/2 + m21*((-2*a11*sp.sin(q11) - 2*a12*sp.sin(q11 + q12) - 2*l21*sp.sin(q11 + 2*q12))*(-a11*sp.sin(q11) - a12*sp.sin(q11 + q12) - l21*sp.sin(q11 + 2*q12)) + (a11*sp.cos(q11) + a12*sp.cos(q11 + q12) + l21*sp.cos(q11 + 2*q12))*(2*a11*sp.cos(q11) + 2*a12*sp.cos(q11 + q12) + 2*l21*sp.cos(q11 + 2*q12)))/2 + m22*((-2*a11*sp.sin(q11) - 2*a12*sp.sin(q11 + q12) - 2*a21*sp.sin(q11 + 2*q12) - 2*l22*sp.sin(q11 + 2*q12 + q22))*(-a11*sp.sin(q11) - a12*sp.sin(q11 + q12) - a21*sp.sin(q11 + 2*q12) - l22*sp.sin(q11 + 2*q12 + q22)) + (a11*sp.cos(q11) + a12*sp.cos(q11 + q12) + a21*sp.cos(q11 + 2*q12) + l22*sp.cos(q11 + 2*q12 + q22))*(2*a11*sp.cos(q11) + 2*a12*sp.cos(q11 + q12) + 2*a21*sp.cos(q11 + 2*q12) + 2*l22*sp.cos(q11 + 2*q12 + q22)))/2, I12 + 2*I21 + 2*I22 + m12*(-l12*(-2*a11*sp.sin(q11) - 2*l12*sp.sin(q11 + q12))*sp.sin(q11 + q12) + l12*(2*a11*sp.cos(q11) + 2*l12*sp.cos(q11 + q12))*sp.cos(q11 + q12))/2 + m21*((-a12*sp.sin(q11 + q12) - 2*l21*sp.sin(q11 + 2*q12))*(-2*a11*sp.sin(q11) - 2*a12*sp.sin(q11 + q12) - 2*l21*sp.sin(q11 + 2*q12)) + (a12*sp.cos(q11 + q12) + 2*l21*sp.cos(q11 + 2*q12))*(2*a11*sp.cos(q11) + 2*a12*sp.cos(q11 + q12) + 2*l21*sp.cos(q11 + 2*q12)))/2 + m22*((-a12*sp.sin(q11 + q12) - 2*a21*sp.sin(q11 + 2*q12) - 2*l22*sp.sin(q11 + 2*q12 + q22))*(-2*a11*sp.sin(q11) - 2*a12*sp.sin(q11 + q12) - 2*a21*sp.sin(q11 + 2*q12) - 2*l22*sp.sin(q11 + 2*q12 + q22)) + (a12*sp.cos(q11 + q12) + 2*a21*sp.cos(q11 + 2*q12) + 2*l22*sp.cos(q11 + 2*q12 + q22))*(2*a11*sp.cos(q11) + 2*a12*sp.cos(q11 + q12) + 2*a21*sp.cos(q11 + 2*q12) + 2*l22*sp.cos(q11 + 2*q12 + q22)))/2],
#[                                                                                                                                                                I12 + 2*I21 + 2*I22 + m12*(-2*l12*(-a11*sp.sin(q11) - l12*sp.sin(q11 + q12))*sp.sin(q11 + q12) + 2*l12*(a11*sp.cos(q11) + l12*sp.cos(q11 + q12))*sp.cos(q11 + q12))/2 + m21*((-2*a12*sp.sin(q11 + q12) - 4*l21*sp.sin(q11 + 2*q12))*(-a11*sp.sin(q11) - a12*sp.sin(q11 + q12) - l21*sp.sin(q11 + 2*q12)) + (2*a12*sp.cos(q11 + q12) + 4*l21*sp.cos(q11 + 2*q12))*(a11*sp.cos(q11) + a12*sp.cos(q11 + q12) + l21*sp.cos(q11 + 2*q12)))/2 + m22*((-2*a12*sp.sin(q11 + q12) - 4*a21*sp.sin(q11 + 2*q12) - 4*l22*sp.sin(q11 + 2*q12 + q22))*(-a11*sp.sin(q11) - a12*sp.sin(q11 + q12) - a21*sp.sin(q11 + 2*q12) - l22*sp.sin(q11 + 2*q12 + q22)) + (2*a12*sp.cos(q11 + q12) + 4*a21*sp.cos(q11 + 2*q12) + 4*l22*sp.cos(q11 + 2*q12 + q22))*(a11*sp.cos(q11) + a12*sp.cos(q11 + q12) + a21*sp.cos(q11 + 2*q12) + l22*sp.cos(q11 + 2*q12 + q22)))/2,                                                                                                                                       I12 + 4*I21 + 4*I22 + m12*(2*l12**2*sp.sin(q11 + q12)**2 + 2*l12**2*sp.cos(q11 + q12)**2)/2 + m21*((-2*a12*sp.sin(q11 + q12) - 4*l21*sp.sin(q11 + 2*q12))*(-a12*sp.sin(q11 + q12) - 2*l21*sp.sin(q11 + 2*q12)) + (a12*sp.cos(q11 + q12) + 2*l21*sp.cos(q11 + 2*q12))*(2*a12*sp.cos(q11 + q12) + 4*l21*sp.cos(q11 + 2*q12)))/2 + m22*((-2*a12*sp.sin(q11 + q12) - 4*a21*sp.sin(q11 + 2*q12) - 4*l22*sp.sin(q11 + 2*q12 + q22))*(-a12*sp.sin(q11 + q12) - 2*a21*sp.sin(q11 + 2*q12) - 2*l22*sp.sin(q11 + 2*q12 + q22)) + (a12*sp.cos(q11 + q12) + 2*a21*sp.cos(q11 + 2*q12) + 2*l22*sp.cos(q11 + 2*q12 + q22))*(2*a12*sp.cos(q11 + q12) + 4*a21*sp.cos(q11 + 2*q12) + 4*l22*sp.cos(q11 + 2*q12 + q22)))/2]])
#def inv(A):
with open('Matrix.pcl','r') as matr:
    A=pickle.load(matr)

#A=sp.Matrix([[9,2,4],[2,9,6],[4,6,9]])
#L_calc=sp.cholesky(A)


    print "..."
    #B=sp.zeros(4)
    #B[0:4,0:4]=A[0:4,0:4]
    #A=B
    n=A.shape[0]
    L=sp.zeros(n)
    
    print "compute lower triangular matrix (Cholesky)"
    
    for i in range(0,n):
        #print "i= ",i
        for j in range(0,i+1):
            #print "j= ",j
            Summe = A[i,j]
            for k in range(0,j):
                #print "k= ",k
                Summe = Summe - L[i,k] * L [j,k]
            if i > j:
                L[i,j] = Summe / L[j,j]
            else:
                L[i,i] = sp.sqrt(Summe)
    det=1
    
    print "calculate determinant and inverse"
    
    for z in range(n):          # ist nur von den Diagonalelementen abhängig
        det=det*L[z,z]
    
    A_inv=(L*L.T).adjugate()/(det**2)     # A^{-1}=\fraq{1}{\det{LL^{T}}}(LL^{T}_{adj}
    #print "compare ..."
    #B=A.subs(params_values)
    #print "1"
    #H=B.subs(para_q)
    #print "2"
    #C=A_inv.subs(params_values)
    #print "3"
    #L=C.subs(para_q)
    #print "4"
    #D=H.inv()
    #print "5"
    #E=D-L
    #print E    
    print "done ... thanks"
    
    IPS()
    #return A_inv      