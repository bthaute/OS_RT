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



from model_2dof import q1_dd_expr, q2_dd_expr, mod1

# <codecell>

#Umformen zu einer lambda function
q1_dd_fnc = sp.lambdify([mod1.qs[0],mod1.qs[1], mod1.qds[0],mod1.qds[1],\
                         mod1.extforce_list[0],mod1.extforce_list[1]],q1_dd_expr,'numpy')
q2_dd_fnc = sp.lambdify([mod1.qs[0],mod1.qs[1], mod1.qds[0],mod1.qds[1],\
                         mod1.extforce_list[0],mod1.extforce_list[1]],q2_dd_expr,'numpy')

# <codecell>

def pd_controller(z):
    q1, q2, q1_d, q2_d = z


    k1 = 100000
    k2 = 100000
    f1 = -k1*(q1+pi/2) - k1*q1_d
    f2 = -k2*q2 - k2*q2_d

    return f1, f2

def get_zd(z,t):
    q1, q2, q1_d, q2_d = z
    f1 = 0 #force1(q1)
    f2 = 0 #force2(q1,q2)
    f1, f2 = pd_controller(z)
    q1_dd = q1_dd_fnc(q1, q2, q1_d, q2_d,f1,f2)
    q2_dd = q2_dd_fnc(q1, q2, q1_d, q2_d,f1,f2)
    #print q1_dd
    #print q2_dd
    return r_[q1_d,q2_d,q1_dd,q2_dd]

# <codecell>

tt = np.linspace(1,30,10000)

z0 = r_[-0.45*pi, 0, 0, 0]

# <codecell>

lsg = odeint(get_zd,z0,tt)
#get_zd(z0,0)


lsg2 = np.concatenate( (tt.reshape(-1, 1) , lsg), axis=1)

from IPython import embed as IPS
#IPS()


# <codecell>

fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.set_title("odeint")
ax.plot(tt,(lsg[:,1]*180/pi), label= r"$q_1$")
ax.plot(tt, (lsg[:,0]*180/pi), label= r"$q_2$")
ax.set_xlabel("$t$ in s")
ax.set_ylabel("$\phi_2$")
ax.legend()
plt.show()
np.save('lsg_outfile',lsg2)

# <codecell>


