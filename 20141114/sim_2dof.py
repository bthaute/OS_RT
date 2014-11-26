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
from parameter import para_g, para_m, para_l, para_a, para_k, para_d, para_I



from model_2dof import q1_dd_expr, q2_dd_expr, mod1
import traj_2dof as traj

print "import finished ..."

params_values = {"m1":para_m[0], "m2":para_m[1], "I1":para_I[0] ,"I2":para_m[1], "a1":para_a[0],
                 "l1": para_l[0], "l2": para_l[1], "k1":para_k[0], "k2":para_k[1], "d1":para_d[0], "d2":para_d[1], "g":para_g }

# <codecell>                 
#sol = sp.solve(mod1.eq_list,mod1.extforce_list)
subslist=zip(mod1.extforce_list,[0,0])
sol = mod1.eq_list.subs(subslist)
F1 = sol[0].subs(params_values)
F2 = sol[1].subs(params_values)
#mod_temp = mod1.eq_list - mod1.eq_list.jacobian(mod1.extforce_list)*mod1.extforce_list
#mod_temp.simplify()
#F1 = mod_temp[0].subs(params_values)
#F2 = mod_temp[1].subs(params_values)

print "force equations created"

#F1 = sol[mod1.extforce_list[0]].subs(params_values)
#F2 = sol[mod1.extforce_list[1]].subs(params_values)
print "build lambda-functions"
F1_fnc = sp.lambdify([mod1.qs[0],mod1.qs[1], mod1.qds[0],mod1.qds[1],mod1.qdds[0],mod1.qdds[1]],F1,'numpy')
F2_fnc = sp.lambdify([mod1.qs[0],mod1.qs[1], mod1.qds[0],mod1.qds[1],mod1.qdds[0],mod1.qdds[1]],F2,'numpy')

# <codecell>

#Umformen zu einer lambda function
q1_dd_fnc = sp.lambdify([mod1.qs[0],mod1.qs[1], mod1.qds[0],mod1.qds[1],\
                         mod1.extforce_list[0],mod1.extforce_list[1]],q1_dd_expr,'numpy')
q2_dd_fnc = sp.lambdify([mod1.qs[0],mod1.qs[1], mod1.qds[0],mod1.qds[1],\
                         mod1.extforce_list[0],mod1.extforce_list[1]],q2_dd_expr,'numpy')

# <codecell>
def calc_F_traj(time):
    t = time
    states=qq_traj(t)  
  
    f1 = F1_fnc(states[0],states[1],states[2],states[3],states[4],states[5])
    f2 = F2_fnc(states[0],states[1],states[2],states[3],states[4],states[5])
        
    return f1, f2

def pd_controller(z,t):
    q1, q2, q1_d, q2_d = z
    #
    states=qq_traj(t)
    k1 = 1e7
    k2 = 1e7
    f1 = k1*(states[0]-q1) + k1*(states[2]-q1_d) + F1_fnc(states[0],states[1],states[2],states[3],states[4],states[5])
    f2 = k2*(states[1]-q2) + k2*(states[3]-q2_d) + F2_fnc(states[0],states[1],states[2],states[3],states[4],states[5])
    #
    return f1, f2

def get_zd(z,t):
    q1, q2, q1_d, q2_d = z
    ti = t
    #f1 = 0
    #f2 = 0
    f1, f2 = pd_controller(z,ti)
    #f1, f2 = calc_F_traj(ti)
    #if ti>25:
    #    f1 -= 100000
    q1_dd = q1_dd_fnc(q1, q2, q1_d, q2_d, f1, f2)
    q2_dd = q2_dd_fnc(q1, q2, q1_d, q2_d, f1, f2)
    #print "q1_dd: ",q1_dd
    #print q2_dd
    return r_[q1_d,q2_d,q1_dd,q2_dd]

# <codecell>

tt = np.linspace(0,40,10000)

z0 = r_[-pi/4, pi/4, 0, 0]
q1_end =3* pi/4
q2_end = pi/2
qq_traj=traj.calc_traj(tt[0],tt[-1]/2,z0[0],q1_end,z0[1],q2_end)

# <codecell>
print "simulate"
lsg = odeint(get_zd,z0,tt)
#get_zd(z0,0)


lsg2 = np.concatenate( (tt.reshape(-1, 1) , lsg), axis=1)

#from IPython import embed as IPS
#IPS()


# <codecell>
deg=180/pi
fig = pl.figure()
ax = fig.add_subplot(1,1,1)
ax.set_title("odeint")
ax.plot(tt,lsg[:,0]*deg, label= r"$q_1$")
ax.plot(tt, lsg[:,1]*deg, label= r"$q_2$")
ax.set_xlabel("$t$ in s")
ax.set_ylabel("$\phi_2$")
ax.legend()
pl.show()
np.save('lsg_outfile',lsg2)

print "thanks"
# <codecell>


