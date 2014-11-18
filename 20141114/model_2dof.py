# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

import sympy as sp
import model_tools as mt
import pickle
#import symb_tools as st
#import numpy as np
#from numpy import r_,pi
#from scipy.integrate import odeint
#import pylab as pl
#import matplotlib.pyplot as plt


pfilepath = "model2dof.pcl"

flag_new_model_generation = True

params_values = {"m1":2550, "m2":1700, "I1":53.125 ,"I2":17.354, "a1":9,
                 "l1": 4.5, "l2": 8, "k1":0, "k2":0, "d1":100, "d2":100, "g":9.81 }


if flag_new_model_generation:
    print "create new model"
    # Variablen definieren
    t = sp.Symbol("t")
    x1 = sp.Symbol("x1")
    q1 = sp.Function("q1")(t)
    q2 = sp.Function("q2")(t)
    qq = sp.Matrix([q1, q2])
    FF = sp.Matrix(sp.symbols("F1, F2"))
    parameter = sp.symbols("I1, I2, m1, m2, l1, l2, a1, k1, k2, d1, d2, g")
    I1, I2, m1, m2, l1, l2, a1, k1, k2, d1, d2, g = parameter

    # Geometrie (karthesische Koordinaten der Schwerpunkte)
    ex = sp.Matrix([1, 0])
    ey = sp.Matrix([0, 1])
    
    # Rotationsmatrix
    mt.Rz(q1)
    # Zwangsbedingungen
    s1 = l1*mt.Rz(q1)*ex
    s2 = a1*mt.Rz(q1)*ex + l2*mt.Rz(q1+q2)*ex
    # definiere neue Variablen für Ableitungen
    q1d = q1.diff(t)
    q2d = q2.diff(t)
    s1d = s1.diff(t)
    s2d = s2.diff(t)

    # kinetische und potentielle Energie
    T = (I1*q1d**2 + I2*(q1d+q2d)**2 + m1*(s1d.T*s1d)[0] + m2*(s2d.T*s2d)[0] )/2
    #T = (I1*q1d**2 + I2*(q2d)**2 + m1*(s1d.T*s1d)[0] + m2*(s2d.T*s2d)[0] )/2
    V = (k1*q1**2 + k2*q2**2)/2 + g*m1*s1[1] + g*m2*s2[1]
    
    # Erzeuge Modell nach Lagrange
    print "create model by model_tools"
    mod1 = mt.generate_model(T, V, qq, FF)
    mod1.eq_list.simplify()
    mod1.eq_list
    
    # Dissipative Kräfte einbeziehen
    q1_d, q2_d =sp.symbols("q1_d,q2_d")
    mod1.eq_list=mod1.eq_list+sp.Matrix([[d1*q1_d],[d2*q2_d]])
    # Enzelne Modell-Gleichungen auslesen
    eq1 = mod1.eq_list[0]#.subs(sublistF)
    eq2 = mod1.eq_list[1]#.subs(sublistF)
    eq1
    #Loese nach q_dot_dot auf
    print "solving equation system"
    M_temp = mod1.eq_list.jacobian(mod1.qdds)
    sol = M_temp.inv()*(-mod1.eq_list + M_temp*mod1.qdds)    
#    sol.simplify();
    #dictionary verwenden um Parameter durch Werte zu ersetzen
    q1_dd_expr = sol[0].subs(params_values)
    q2_dd_expr = sol[1].subs(params_values)
    q1_dd_expr.simplify()
    q2_dd_expr.simplify()
    # Speichere Modell in Datei ab.
    print "saving model"
    pdict = {"q1_dd_expr": q1_dd_expr, "q2_dd_expr": q2_dd_expr, "mod1": mod1}

    with  open(pfilepath, "w") as pfile:
        pickle.dump(pdict, pfile)

else:  # flag_new_model_generation == False

    print "just open exsiting model"
    with  open(pfilepath, "r") as pfile:
        pdict = pickle.load(pfile)

    q1_dd_expr = pdict["q1_dd_expr"]
    q2_dd_expr = pdict["q2_dd_expr"]
    mod1 = pdict["mod1"]

print "done"