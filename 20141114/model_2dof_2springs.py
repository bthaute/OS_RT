# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

import sympy as sp
#import numpy as np
import model_tools as mt
import pickle
from parameter_springs import para_g, para_m, para_l, para_a, para_k, para_d, para_I
from IPython import embed as IPS
import cholesky as chol

pfilepath = "model2dof2springs.pcl"

flag_new_model_generation = False
params_values = {"m11":para_m[0,0], "m12":para_m[0,1], "m21":para_m[1,0], "m22":para_m[1,1], "I11":para_I[0,0] ,"I12":para_I[0,1], "I21":para_I[1,0] ,"I22":para_I[1,1], "a11":para_a[0,0], "a12":para_a[0,1], "a21":para_a[1,0], "a22":para_a[1,1],
                 "l11": para_l[0,0], "l12": para_l[0,1], "l21": para_l[1,0], "l22": para_l[1,1], "k1":para_k[0], "k2":para_k[1], "d1":para_d[0], "d2":para_d[1], "g":para_g }

if flag_new_model_generation:
    print "create new model"
    # Variablen definieren
    t = sp.Symbol("t")
    x1 = sp.Symbol("x1")
    q11 = sp.Function("q11")(t)
    q12 = sp.Function("q12")(t)
    q21 = sp.Function("q21")(t)
    q22 = sp.Function("q22")(t)
    qq = sp.Matrix([q11, q12, q21, q22])
    FF = sp.Matrix(sp.symbols("F11, F12, F21, F22"))
    parameter = sp.symbols("I11, I12, I21, I22, m11, m12, m21, m22, l11, l12, l21, l22, a11, a12, a21, a22, k1, k2, d1, d2, g")
    I11, I12, I21, I22, m11, m12, m21, m22, l11, l12, l21, l22, a11, a12, a21, a22, k1, k2, d1, d2, g = parameter

    # Geometrie (karthesische Koordinaten der Schwerpunkte)
    ex = sp.Matrix([1, 0])
    ey = sp.Matrix([0, 1])
    
    # Schwerpunkte 
    s11 = l11*mt.Rz(q11)*ex
    s12 = a11*mt.Rz(q11)*ex + l12*mt.Rz(q11+q12)*ex
    s21 = a11*mt.Rz(q11)*ex + a12*mt.Rz(q11+q12)*ex + l21*mt.Rz(q11+q12+q21)*ex
    s22 = a11*mt.Rz(q11)*ex + a12*mt.Rz(q11+q12)*ex + a21*mt.Rz(q11+q12+q21)*ex + l22*mt.Rz(q11+q12+q21+q22)*ex
    
    # definiere neue Variablen für Ableitungen
    q11d = q11.diff(t)
    q12d = q12.diff(t)
    q21d = q21.diff(t)
    q22d = q22.diff(t)
    s11d = s11.diff(t)
    s12d = s12.diff(t)
    s21d = s21.diff(t)
    s22d = s22.diff(t)

    # kinetische und potentielle Energie
    T = (I11*q11d**2 + I12*(q11d+q12d)**2 + I21*(q11d+q12d+q21d)**2 + I22*(q11d+q12d+q21d+q22d)**2 + m11*(s11d.T*s11d)[0] + m12*(s12d.T*s12d)[0] + m21*(s21d.T*s21d)[0] + m22*(s22d.T*s22d)[0])/2
    #T = (I1*q1d**2 + I2*(q2d)**2 + m1*(s1d.T*s1d)[0] + m2*(s2d.T*s2d)[0] )/2
    V = (k1*q12**2 + k2*q22**2)/2 + g*m11*s11[1] + g*m12*s12[1] + g*m21*s21[1] + g*m22*s22[1]
    
    # Erzeuge Modell nach Lagrange
    print "create model by model_tools"
    mod1 = mt.generate_model(T, V, qq, FF)
#    mod1.eq_list.simplify()

    
    # Dissipative Kräfte einbeziehen
    q11_d, q12_d,q21_d, q22_d =sp.symbols("q11_d, q12_d,q21_d, q22_d")
    mod1.eq_list=mod1.eq_list+sp.Matrix([[d1*q11_d],[d2*q12_d],[d1*q21_d],[d2*q22_d]])

    #Loese nach q_dot_dot auf
    eq_list_temp = mod1.eq_list.subs(params_values)
    MM_temp = mod1.MM.subs(params_values)
    #MM_temp.simplify()
    print "solving equation system"
    subslist = zip(mod1.qdds, [0,0,0,0])
    temp = eq_list_temp.subs(subslist)
    sol = chol.inv(MM_temp)*(-temp) 

#    mt.solve_eq(mod1)      
#    q11_dd_expr = mod1.solved_eq[0].rhs
#    q12_dd_expr = mod1.solved_eq[1].rhs
#    q21_dd_expr = mod1.solved_eq[2].rhs
#    q22_dd_expr = mod1.solved_eq[3].rhs
#    q11_dd_expr.simplify()
#    q12_dd_expr.simplify()
#    q21_dd_expr.simplify()
#    q22_dd_expr.simplify()
    #dictionary verwenden um Parameter durch Werte zu ersetzen
    print "craete system q_dd_expr"
    q11_dd_expr = sol[0]
    q12_dd_expr = sol[1]
    q21_dd_expr = sol[2]
    q22_dd_expr = sol[3]
#    q11_dd_expr.simplify()
#    q12_dd_expr.simplify()
#    q21_dd_expr.simplify()
#    q22_dd_expr.simplify()
    # Speichere Modell in Datei ab.
    print "saving model"
    pdict = {"q11_dd_expr": q11_dd_expr, "q12_dd_expr": q12_dd_expr, "q21_dd_expr": q21_dd_expr, "q22_dd_expr": q22_dd_expr, "mod1": mod1}

    with  open(pfilepath, "w") as pfile:
        pickle.dump(pdict, pfile)

else:  # flag_new_model_generation == False

    print "just open exsiting model"
    with  open(pfilepath, "r") as pfile:
        pdict = pickle.load(pfile)

    q11_dd_expr = pdict["q11_dd_expr"]
    q12_dd_expr = pdict["q12_dd_expr"]
    q21_dd_expr = pdict["q21_dd_expr"]
    q22_dd_expr = pdict["q22_dd_expr"]
    mod1 = pdict["mod1"]

print "read model"
#IPS()