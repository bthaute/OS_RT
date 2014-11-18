# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import sympy as sp
import symb_tools as st
import model_tools as mt
import numpy as np
from numpy import r_,pi
from scipy.integrate import odeint
import pylab as pl
import matplotlib.pyplot as plt

import pickle

pfilepath = "model2dof.pcl"

flag_new_model_generation = True

params_values = {"m1":2550, "m2":1700, "I1":53.125 ,"I2":17.354, "a1":9,
                 "l1": 4.5, "l2": 8, "k1":0, "k2":0, "d1":100, "d2":100, "g":9.81 }


if flag_new_model_generation:

    # <codecell>

    t = sp.Symbol("t")
    x1 = sp.Symbol("x1")

    # <codecell>

    q1 = sp.Function("q1")(t)
    q2 = sp.Function("q2")(t)
    qq = sp.Matrix([q1, q2])
    FF = sp.Matrix(sp.symbols("F1, F2"))
    parameter = sp.symbols("I1, I2, m1, m2, l1, l2, a1, k1, k2, d1, d2, g")
    #Bill: was bewirkt diese Zeile? Quasi um jedes Symbol, was oben angelegt wurde als Listenelemente ansprechen zu können
    I1, I2, m1, m2, l1, l2, a1, k1, k2, d1, d2, g = parameter

    # <markdowncell>

    # Geometrie (karthesische Koordinaten der Schwerpunkte)

    # <codecell>

    ex = sp.Matrix([1, 0])
    ey = sp.Matrix([0, 1])

    # <codecell>

    # Rotationsmatrix
    mt.Rz(q1)

    # <codecell>

    s1 = l1*mt.Rz(q1)*ex
    s1

    # <codecell>

    s2 = a1*mt.Rz(q1)*ex + l2*mt.Rz(q1+q2)*ex
    s2

    # <markdowncell>

    # kinetische und potentielle Energie

    # <codecell>

    q1d = q1.diff(t)
    q2d = q2.diff(t)
    s1d = s1.diff(t)
    s2d = s2.diff(t)

    # <codecell>

    T = (I1*q1d**2 + I2*(q1d+q2d)**2 + m1*(s1d.T*s1d)[0] + m2*(s2d.T*s2d)[0] )/2
    #T = (I1*q1d**2 + I2*(q2d)**2 + m1*(s1d.T*s1d)[0] + m2*(s2d.T*s2d)[0] )/2

    # <codecell>

    #V = (k1*q1**2 + k2*q2**2)/2
    V = (k1*q1**2 + k2*q2**2)/2 + g*m1*s1[1] + g*m2*s2[1]
    # <codecell>
    m2*((s2.T*ey)[0] + a1 + l2)
    # <codecell>

    mod1 = mt.generate_model(T, V, qq, FF)
    mod1.eq_list.simplify()
    mod1.eq_list
    # <codecell>
    # Dissipative Kräfte einbeziehen
    q1_d, q2_d =sp.symbols("q1_d,q2_d")
    mod1.eq_list=mod1.eq_list+sp.Matrix([[d1*q1_d],[d2*q2_d]])
    #mod1.eq_list
    # <codecell>

    Mq = mod1.MM.expand()

    # <codecell>

    #Mq.simplify()
    Mq

    # <markdowncell>

    # ---
    # Alternativer Zugang (Koordinatendefinition wie Prof. Janschek)

    # <codecell>

    z1 = sp.Function("z1")(t)
    z2 = sp.Function("z2")(t)
    zz = sp.Matrix([z1, z2])
    # Überschreiben von s1, s2

    s1 = l1*mt.Rz(z1)*ex
    s2 = a1*mt.Rz(z1)*ex + l2*mt.Rz(z2)*ex

    z1d = z1.diff(t)
    z2d = z2.diff(t)
    s1d = s1.diff(t)
    s2d = s2.diff(t)

    T = (I1*z1d**2 + I2*z2d**2 + m1*st.norm2(s1d) + m2*st.norm2(s2d) )/2
    V = (k1*z1**2 + k2*z2**2)/2

    # <codecell>

    mod2 = mt.generate_model(T, V, zz, FF)
    mod2.eq_list.simplify()
    mod2.eq_list

    # <codecell>

    s1d

    # <markdowncell>

    # Das ist die Massenmatrix, die bei Prof. Janschek rauskam

    # <codecell>

    Mz = mod2.MM
    Mz

    # <markdowncell>

    # ---

    # <headingcell level=3>

    # Frage: wie kann man die beiden Massenmatrizen ineinander umrechnen?

    # <markdowncell>

    # **Variante 1**: Schritt für Schritt: wir nehmen unsere Bewegungsgleichungen (mod1), substituieren die alten $q$ Variablen durch die neuen $z$ Variablen (auch für alle Ableitungen), bilden dann die Jacobi-Matrix bezüglich $\ddot z$ und bilden schiließlich eine Linearkombination der Gleichungen (Zeilen der Matrix), so dass die resultierende Matrix wieder symmetrisch ist.

    # <codecell>

    # Transformationsmatrix z = S * q:
    S = sp.Matrix([[1, 0], [1,1]])
    S

    # <markdowncell>

    # $q$-Koordinaten durch $z$-Koordinaten ersetzen

    # <codecell>

    q_new = S.inv()*mod2.qs
    qd_new = S.inv()*mod2.qds
    qdd_new = S.inv()*mod2.qdds

    # Kontrollanzeige:
    st.col_stack(q_new, qd_new, qdd_new)

    # <codecell>

    subslist = zip(mod1.qdds, qdd_new) + zip(mod1.qds, qd_new) + zip(mod1.qs, q_new)
    subslist

    # <markdowncell>

    # temporäre Zwischenschritt-Matrix (sie ist noch nicht symmetrisch)

    # <codecell>

    M_tmp = mod1.eq_list.subs(subslist).jacobian(mod2.qdds)
    M_tmp

    # <markdowncell>

    # Multiplikation von links (=linearkombination der Zeile (also der Gleichungen))

    # <codecell>

    C = sp.Matrix([[1, -1], [0, 1]])

    # <codecell>

    res = C*M_tmp
    res.simplify()
    res

    # <codecell>

    res == Mz

    # <markdowncell>

    # ---
    # Das klappt offensichtlich
    #
    # **Variante 2**: (so wollte ich es ursprünglich machen):
    # Wir wissen dass folgende Zusammenhange gelten:
    # $$ z = S q$$
    # $$ \ddot z = S \ddot q$$
    # $$ \ddot q = S^{-1} \ddot z.$$
    #
    # Deswegen können wir auch schreiben:
    #
    # $$M \ddot q + ... = M S^{-1} \ddot z + ... \enspace.$$
    #
    # Diese Matrix können wir dann auch von links mit einener Matrix multiplizieren (Linearkombination der Zeilen, also Gleichungen)
    # Dazu wählt man die *Transponierte* der Matrix auf der rechten Seite. **Hier lag mein Fehler**, weil ich die *Inverse* von $S$ genommen hatte. Ma bekommt
    #
    # $$\underbrace{ (S^{-1})^T M S^{-1} }_{\hat M} \ddot z + ... \enspace .$$
    #
    # Das ist zwar etwas länger zu erklären aber einfacher zu rechnen.

    # <codecell>

    res2 = S.inv().T*Mq*S.inv()
    res2 == Mz

    # <markdowncell>

    # logisch. Die Massenmatrix hängt ja noch von $q$ bzw. $z$ ab.

    # <codecell>

    res2.subs(subslist) == Mz

    # <markdowncell>

    # Auflösen der Bewegungsgleichungen nach den zweifachen Winkelableitungen

    # <codecell>

    #q1,q2=sp.symbols("q1,q2") # erstmal keine Zeitabhängigkeit

    # <codecell>

    #externe Kräfte einbeziehen ohne Dissipationskräfte D*q_d
    #F1 = (m1*l1+m2*a1)*sp.cos(q1)*9.81
    #F2 = m2*l2*sp.cos(q1+q2)*9.81

    # <codecell>

    #sublistF=zip(FF,[(m1*l1+m2*a1)*sp.cos(q1)*9.81,m2*l2*sp.cos(q1+q2)*9.81])

    # <codecell>

    eq1 = mod1.eq_list[0]#.subs(sublistF)
    eq2 = mod1.eq_list[1]#.subs(sublistF)
    eq1

    # <codecell>

    print "solving equation system"
    sol = sp.solve([eq1,eq2],mod1.qdds)
    sol

    # <codecell>

    #dictionary verwenden um Parameter durch Werte zu ersetzen


    # <codecell>

    q1_dd_expr = sol[mod1.qdds[0]].subs(params_values)
    q2_dd_expr = sol[mod1.qdds[1]].subs(params_values)


    pdict = {"q1_dd_expr": q1_dd_expr, "q2_dd_expr": q2_dd_expr, "mod1": mod1}



    with  open(pfilepath, "w") as pfile:
        pickle.dump(pdict, pfile)


else:  # flag_new_model_generation == False


    with  open(pfilepath, "r") as pfile:
        pdict = pickle.load(pfile)

    q1_dd_expr = pdict["q1_dd_expr"]
    q2_dd_expr = pdict["q2_dd_expr"]
    mod1 = pdict["mod1"]

