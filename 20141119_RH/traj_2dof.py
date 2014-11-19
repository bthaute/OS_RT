# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import sympy as sp
import symb_tools as st
from sympy import pi
import numpy as np
import pylab as pl

# <codecell>

def calc_traj(T0,T1,q1_t0,q1_t1,q2_t0,q2_t1):
    

    t = sp.Symbol('t')
    
    #T0 = 0
    #T1 = 10
    #q1_t0=-pi/2
    #q1_t1=pi/4
    #q2_t0=0
    #q2_t1=-pi/4
    
    # Übergangspolynom bestimmen
    q1_poly = st.trans_poly(t, cn=2, left=(T0, q1_t0, 0, 0), right=(T1, q1_t1, 0, 0))
    q1_poly
    
    # <codecell>
    
    # Stückweise definierten Ausdruck anlegen
    q1_piecewise = sp.Piecewise( (q1_poly.subs(t, T0), t < T0), (q1_poly, t < T1), (q1_poly.subs(t, T1), True))
    q1_piecewise
    
    # <codecell>
    
    q1_piecewise_d=q1_piecewise.diff('t')
    q1_piecewise_dd=q1_piecewise_d.diff('t')
    
    # <markdowncell>
    
    # → Siehe auch Doku (Eintippen von `sp.Piecewise?`)
    
    # <codecell>
    
    # Das gleiche für q2:
    q2_poly = st.trans_poly(t, cn=2, left=(T0, q2_t0, 0, 0),
                                     right=(T1, q2_t1, 0, 0))
    
    q2_piecewise = sp.Piecewise( (q2_poly.subs(t, T0), t < T0), (q2_poly, t < T1), (q2_poly.subs(t, T1), True))
    
    # <codecell>
    
    q2_piecewise_d=q2_piecewise.diff('t')
    q2_piecewise_dd=q2_piecewise_d.diff('t')
    
    # <codecell>
    
    # als Spaltenvektor:
    qq_traj = sp.Matrix([[q1_piecewise, q2_piecewise],[q1_piecewise_d, q2_piecewise_d],[q1_piecewise_dd, q2_piecewise_dd]])
    qq_traj
    
    # <codecell>
    
    # Ausführbare Python-Funktion erstellen, die man auch mit numpy-Vektoren auswerten kann
    qq_func = st.expr_to_func(t, list(qq_traj), eltw_vectorize=True)
    
    #qq = qq_func(tt)
    
    return qq_func

# <markdowncell>

# Eine solche Funktion kann man jetzt entweder an bestimmten Zeitpunkten auswerten und plotten oder in der RHS-Funktion der Simulation für eine Vorsteuerung verwenden

# <codecell>
def example():
    tt = np.linspace(0, 10, 1e3)
    
    
    T0 = 0
    T1 = 10
    q1_t0=-pi/2
    q1_t1=pi/4
    q2_t0=0
    q2_t1=-pi/4
    
    
    # <codecell>
    qq_traj=calc_traj(T0,T1,q1_t0,q1_t1,q2_t0,q2_t1)
    qq = qq_traj(tt)
    #np.save('traj_01',qq)
    
    # <codecell>
    
    # Einbetten von Graphen direkt in den Browser
    #matplotlib inline 
    
    # <codecell>
    
    pl.plot(tt, qq[:, 0], label='$q_1$')
    pl.plot(tt, qq[:, 1], label='$q_1_d$')
    #pl.plot(tt, qq[:, 4], label='$q_1_dd$')
    #pl.plot(tt, qq[:, 1], label='$q_2$')
    pl.show()
    

