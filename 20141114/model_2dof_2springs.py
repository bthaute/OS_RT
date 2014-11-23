# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

import sympy as sp
import model_tools as mt
import pickle
from parameter_springs import para_g, para_m, para_l, para_a, para_k, para_d, para_I

pfilepath = "model2dof2springs.pcl"

flag_new_model_generation = True
params_values = {"m11":para_m[0,0], "m12":para_m[0,1], "m21":para_m[1,0], "m22":para_m[1,1], "I11":para_I[0,0] ,"I12":para_m[0,1], "I21":para_I[1,0] ,"I22":para_m[1,1], "a11":para_a[0,0], "a12":para_a[0,1], "a21":para_a[1,0], "a22":para_a[1,1],
                 "l11": para_l[0,0], "l12": para_l[0,1], "l21": para_l[1,0], "l22": para_l[1,1], "k1":para_k[0], "k2":para_k[1], "d1":para_d[0], "d2":para_d[1], "g":para_g }

if flag_new_model_generation:
    print "create new model"
    # Variablen definieren
    t = sp.Symbol("t")
    x1 = sp.Symbol("x1")
    q11 = sp.Function("q11")(t)
    q12 = sp.Function("q12")(t)
    q21 = sp.Function("q12")(t)
    q22 = sp.Function("q22")(t)
    qq = sp.Matrix([q11, q12, q21, q22])
    FF = sp.Matrix(sp.symbols("F11, F12, F21, F22"))
    parameter = sp.symbols("I11, I12, I21, I22, m11, m12, m21, m22, l11, l12, l21, l22, a11, a12, a21, a21, k1, k2, d1, d2, g")
    I11, I12, I21, I22, m11, m12, m21, m22, l11, l12, l21, l22, a11, a12, a21, a21, k1, k2, d1, d2, g = parameter

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
    V = (k1*q11**2 + k2*q21**2)/2 + g*m11*s11[1] + g*m12*s12[1] + g*m21*s21[1] + g*m22*s22[1]
    
    # Erzeuge Modell nach Lagrange
    print "create model by model_tools"
    mod1 = mt.generate_model(T, V, qq, FF)
#    mod1.eq_list.simplify()

    
    # Dissipative Kräfte einbeziehen
    q11_d, q12_d,q21_d, q22_d =sp.symbols("q11_d, q12_d,q21_d, q22_d")
    mod1.eq_list=mod1.eq_list+sp.Matrix([[d1*q11_d],[d1*q12_d],[d2*q21_d],[d1*q22_d]])
    # Enzelne Modell-Gleichungen auslesen
    eq1 = mod1.eq_list[0]#.subs(sublistF)
    eq2 = mod1.eq_list[1]#.subs(sublistF)

    #Loese nach q_dot_dot auf
    print "solving equation system"
    M_temp = mod1.eq_list.jacobian(mod1.qdds)
    temp = -mod1.eq_list + M_temp*mod1.qdds
#    temp.simplify()
    sol = M_temp.inv()*(temp)       
    #dictionary verwenden um Parameter durch Werte zu ersetzen
    q1_dd_expr = sol[0].subs(params_values)
    q2_dd_expr = sol[1].subs(params_values)
#    q1_dd_expr.simplify()
#    q2_dd_expr.simplify()
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

print "read model"