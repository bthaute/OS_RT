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
from parameter_springs import para_g, para_m, para_l, para_a, para_k, para_d, para_I
traj=np.load('traj_01.npy')


params_values = {"m11":para_m[0,0], "m12":para_m[0,1], "m21":para_m[1,0], "m22":para_m[1,1], "I11":para_I[0,0] ,"I12":para_I[0,1], "I21":para_I[1,0] ,"I22":para_I[1,1], "a11":para_a[0,0], "a12":para_a[0,1], "a21":para_a[1,0], "a22":para_a[1,1],
                 "l11": para_l[0,0], "l12": para_l[0,1], "l21": para_l[1,0], "l22": para_l[1,1], "k1":para_k[0], "k2":para_k[1], "d1":para_d[0], "d2":para_d[1], "g":para_g }

# <codecell>                 
# Stelle Modell nach externen Groessen um

subslist=zip(mod1.extforce_list,[0,0,0,0])
sol = mod1.eq_list.subs(subslist)
F11 = sol[0].subs(params_values)
F12 = sol[1].subs(params_values)
F21 = sol[2].subs(params_values)
F22 = sol[3].subs(params_values)
# <codecell>

#Umformen zu einer lambda function
q11_dd_fnc = sp.lambdify([mod1.qs[0],mod1.qs[1], mod1.qs[2],mod1.qs[3], mod1.qds[0],mod1.qds[1],\
                         mod1.qds[2],mod1.qds[3],mod1.extforce_list[0],mod1.extforce_list[1],\
                         mod1.extforce_list[2],mod1.extforce_list[3]],q11_dd_expr,'numpy')
q12_dd_fnc = sp.lambdify([mod1.qs[0],mod1.qs[1], mod1.qs[2],mod1.qs[3], mod1.qds[0],mod1.qds[1],\
                         mod1.qds[2],mod1.qds[3],mod1.extforce_list[0],mod1.extforce_list[1],\
                         mod1.extforce_list[2],mod1.extforce_list[3]],q12_dd_expr,'numpy')
q21_dd_fnc = sp.lambdify([mod1.qs[0],mod1.qs[1], mod1.qs[2],mod1.qs[3], mod1.qds[0],mod1.qds[1],\
                         mod1.qds[2],mod1.qds[3],mod1.extforce_list[0],mod1.extforce_list[1],\
                         mod1.extforce_list[2],mod1.extforce_list[3]],q21_dd_expr,'numpy')
q22_dd_fnc = sp.lambdify([mod1.qs[0],mod1.qs[1], mod1.qs[2],mod1.qs[3], mod1.qds[0],mod1.qds[1],\
                         mod1.qds[2],mod1.qds[3],mod1.extforce_list[0],mod1.extforce_list[1],\
                         mod1.extforce_list[2],mod1.extforce_list[3]],q22_dd_expr,'numpy')
# <codecell>
def calc_F_traj(t):
       
    sub=zip(['q1'],[traj[t,0]])+zip(['q2'],[traj[t,1]])+zip(['q1_d'],[traj[t,2]])+zip(['q2_d'],\
    [traj[t,3]])+zip(['q1_dd'],[traj[t,4]])+zip(['q2_dd'],[traj[t,5]])
    
    f1=F1.subs(sub)
    f2=F2.subs(sub)
    
    return f1, f2

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
#    f1, f2 = calc_F_traj(t)
    q1_dd = q1_dd_fnc(q1, q2, q1_d, q2_d,f1,f2)
    q2_dd = q2_dd_fnc(q1, q2, q1_d, q2_d,f1,f2)
    #print q1_dd
    #print q2_dd
    return r_[q1_d,q2_d,q1_dd,q2_dd]

# <codecell>

tt = np.linspace(1,30,10000)

z0 = r_[-0.0*pi, 0, 0, 0]

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


