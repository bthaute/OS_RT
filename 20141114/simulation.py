# -*- coding: utf-8 -*-
"""
Created on Tue Dec  9 16:58:04 2014

@author: herrmann
"""
# <codecell>

import sympy as sp
import symb_tools as st
import model_tools as mt
import numpy as np
from numpy import r_,pi
from scipy.integrate import odeint
import pylab as pl
import matplotlib.pyplot as plt
from IPython import embed as IPS

import traj_2dof as traj
from parameter_springs import para_g, para_m, para_a, para_k, para_d, para_I

from model import mod1

# Number of the actuated joints -> uses .simplify()
nr_aj = 2  

II = sp.Matrix([sp.symbols("I1")])
mm = sp.Matrix([sp.symbols("m1")])
aa = sp.Matrix([sp.symbols("a1")])
kk = sp.Matrix([sp.symbols("k1")])
dd = sp.Matrix([sp.symbols("d1")])
qq_0 = sp.Matrix([sp.symbols("q1_0")])
g = sp.Matrix([sp.Symbol("g")])
for index in range(2,2*(nr_aj)+1):
    II = st.row_stack(II,sp.Matrix([sp.symbols("I"+np.str(index))]))
    mm = st.row_stack(mm,sp.Matrix([sp.symbols("m"+np.str(index))]))
    aa = st.row_stack(aa,sp.Matrix([sp.symbols("a"+np.str(index))]))
    kk = st.row_stack(kk,sp.Matrix([sp.symbols("k"+np.str(index))]))
    dd = st.row_stack(dd,sp.Matrix([sp.symbols("d"+np.str(index))]))
    qq_0 = st.row_stack(qq_0,sp.Matrix([sp.symbols("q"+np.str(index)+"_0")]))

# Modellfehler in der Parametrierung
model_error_on = False
if model_error_on:
    model_error=0.2*np.random.randn(10,6)+1
else:
    model_error=np.ones((10,6))
    
params_values_error = zip(mm,para_m*model_error[:,0]) + zip(II,para_I*model_error[:,1]) + zip(aa,para_a*model_error[:,2]) + zip(kk,para_k*model_error[:,3]) + zip(dd,para_d*model_error[:,4]) + zip(g,para_g*model_error[0,5])
params_values = zip(mm,para_m) + zip(II,para_I) + zip(aa,para_a) + zip(kk,para_k) + zip(dd,para_d) + zip(g,para_g)

# <codecell>                 
# Stelle Modell nach externen Groessen um

subslist_force=zip(mod1.extforce_list,sp.zeros(mod1.qdds.shape[0],1))
sol = mod1.eq_list.subs(subslist_force)
FF = sp.Matrix([sol[0].subs(params_values_error)])
for index in range(1,nr_aj):
    FF = st.row_stack(FF,sp.Matrix([sol[index*2].subs(params_values_error)]))

# <codecell>
print "build lambda-functions"

FF_fnc =    sp.lambdify([mod1.qs[0],mod1.qs[1], mod1.qs[2],mod1.qs[3], mod1.qds[0],mod1.qds[1],\
                         mod1.qds[2],mod1.qds[3],mod1.qdds[0],mod1.qdds[1],mod1.qdds[2],mod1.qdds[3]],FF,'numpy')

# <codecell>

# alternativer Zugang
subslist_M = zip(mod1.qdds, sp.zeros(mod1.qdds.shape[0],1))
temp = -mod1.eq_list.subs(subslist_M).subs(params_values)
Mq_dd_func = sp.lambdify([mod1.qs[0],mod1.qs[1], mod1.qs[2],mod1.qs[3], mod1.qds[0],mod1.qds[1],\
                         mod1.qds[2],mod1.qds[3],mod1.extforce_list[0],mod1.extforce_list[1],\
                         mod1.extforce_list[2],mod1.extforce_list[3]],temp,'numpy')
M=mod1.MM.subs(params_values)
M_func = sp.lambdify([mod1.qs[0],mod1.qs[1], mod1.qs[2],mod1.qs[3], mod1.qds[0],mod1.qds[1],\
                         mod1.qds[2],mod1.qds[3],mod1.extforce_list[0],mod1.extforce_list[1],\
                         mod1.extforce_list[2],mod1.extforce_list[3]],M,'numpy')
#IPS()
# <codecell>
def calc_F_traj(t):
       
    states=qq_traj(t)  
  
    f1 = FF_fnc(states[0],0,states[1],0,states[2],0,states[3],0,states[4],0,states[5],0)[0]
    f2 = FF_fnc(states[0],0,states[1],0,states[2],0,states[3],0,states[4],0,states[5],0)[1]
     
    return f1, f2
    
Bk_opt=np.array([[ -3850742.68840943,  -1168179.30944712,    572685.92341211,   -829342.73665955,  -2321800.26669223,   -958161.40449268,    -34473.98775092,   -248108.98384053],
       [ -1896612.40241954,   -732166.89093391,   -181480.50007226,   -537116.708522  ,  -1138256.78085381,  19388676.83651922,   -188894.41163585,   -169773.20668387],
       [  -385962.37262645,   -354504.81738105,   -753379.34580358,   -329553.21279944,   -240036.62640163,   -222917.85184836,   -305371.60841489,   -114109.99513054],
       [  -117495.28088422,   -104823.25767166,   -227355.66425041,   -104686.26100281,    -74138.88126445,    -68152.74985214,    -93554.54750739,   9968492.99219859]])
# <codecell>
def pd_controller(z,t):
    q1, q2, q3, q4, q1_d, q2_d, q3_d, q4_d = z
    
    k=1e7*np.ones((nr_aj,1))
    states=qq_traj(t)
    zstern=[states[0],q2,states[1],q4,states[2],q2_d,states[3],q4_d]
    states2=states[:4:]
    states2=states2.reshape((2,-1))
    #states2=states2[::,:nr_aj:]
    # Fehler
    e=states2-z[::2].reshape((2,-1))
    e=e.T
    man=np.ones((nr_aj,1))

    #f = (k*e).dot(man)+FF_fnc(states[0],0,states[1],0,states[2],0,states[3],0,states[4],0,states[5],0)
    f = -np.dot(Bk_opt,z-zstern)[0::2]+FF_fnc(states[0],0,states[1],0,states[2],0,states[3],0,states[4],0,states[5],0)
    return f[0,0],f[1,0]
# <codecell>
def get_zd(z,t):
    q1, q2, q3, q4, q1_d, q2_d, q3_d, q4_d = z
    f1 = 0 #force1(q1)
    f2 = 0 #force2(q1,q2)
    f1, f2 = pd_controller(z,t)
    
    # Stellgrößenbegrenzung
    #f1 = np.clip(f1,-3e5,3e5)
    #f2 = np.clip(f2,-3e5,3e5)
    #global f,temp_zd
    #if temp_zd:
    #    f=np.array([[f1,f2]])
    #    temp_zd = False
    #else: f = np.concatenate((f,np.array([[f1,f2]])),axis=0)
    #q0=r_[pi/3, 0, -pi/2, 0, 0, 0.0, 0,0]
    #qq_dd=mod1.A_para.subs(sub_qq_0)*(z-q0).reshape((-1,1))+mod1.B_para.subs(sub_qq_0)*np.array([[f1],[0],[f2],[0]])
    A_schlange=A#A_fnc(q1,q2,q3,q4)
    B_schlange=B#B_fnc(q1,q2,q3,q4)
    #qq_dd=np.dot(A,(z-q0))+np.dot(B_schlange,np.array([f1,0,f2,0]))
    qq_dd=np.dot(A_schlange,(z-z1))+np.dot(B_schlange,np.array([f1,0,f2,0]))
    #rhs_eq = Mq_dd_func(q1, q2, q3,q4, q1_d, q2_d, q3_d, q4_d,f1,0,f2,0)
    #M_eq = M_func(q1, q2, q3,q4, q1_d, q2_d, q3_d, q4_d,f1,0,f2,0)
    #qq_dd = M_eq**-1*rhs_eq
        
    #print q1_dd
    #print q2_dd
    #print t
    #return r_[q1_d,q2_d,q3_d,q4_d,q1_dd,q2_dd,q3_dd,q4_dd]
    #return r_[q1_d,q2_d,q3_d,q4_d,qq_dd[0,0],qq_dd[1,0],qq_dd[2,0],qq_dd[3,0]]
    return qq_dd
    
# <codecell>

tt = np.linspace(0,10,1000)

error_q11 = 0
error_q21 = 0
q11_t0 = pi/3 # -pi/4
q21_t0 = -pi/2 # pi/4
z0 = r_[q11_t0 + error_q11, 0, q21_t0 + error_q21, 0, 0.08, 0.0, 0,0.056]
z1 = r_[q11_t0 + error_q11, 0, q21_t0 + error_q21, 0, 0, 0.0, 0,0]

q1_end = q11_t0 # 3*pi/4
q2_end = q21_t0 # pi/2
#qq_traj=traj.calc_traj(tt[0],tt[-1]/2,z0[0],q1_end,z0[2],q2_end)
qq_traj=traj.calc_traj(tt[0],tt[-1]/2,q11_t0,q1_end,q21_t0,q2_end)

sub_qq_0=zip(['q1_0'],[q11_t0])+zip(['q2_0'],[0])+zip(['q3_0'],[q21_t0])+zip(['q4_0'],[0])
#A=st.to_np(mod1.A.subs(params_values).subs(sub_qq_0))
A=mod1.A.subs(params_values)
A=A.subs(sub_qq_0)
A=st.to_np(A)
#A_fnc = sp.lambdify(['q1_0','q2_0','q3_0','q4_0'],A,'numpy')
#B=st.to_np(mod1.B.subs(params_values).subs(sub_qq_0))
B=mod1.B.subs(params_values)
B=B.subs(sub_qq_0)
B=st.to_np(B)
#B_fnc = sp.lambdify(['q1_0','q2_0','q3_0','q4_0'],B,'numpy')
q0=r_[pi/3, 0, -pi/2, 0, 0, 0.0, 0,0]
# <codecell>
print "simulate"
lsg = odeint(get_zd,z0,tt)
#get_zd(z0,0)


lsg2 = np.concatenate( (tt.reshape(-1, 1) , lsg), axis=1)

#Regelfehler
qq_soll = qq_traj(tt)
e1 = qq_soll[:,0]-lsg2[:,1]
e1_d = qq_soll[:,2]-lsg2[:,5]
e2 = qq_soll[:,1]-lsg2[:,3]
e2_d = qq_soll[:,3]-lsg2[:,7]
#IPS()

# später für das linearisierte Modell
sub_qq_0=zip(['q1_0'],[q11_t0])+zip(['q2_0'],[0])+zip(['q3_0'],[q21_t0])+zip(['q4_0'],[0])

# <codecell>
deg=180/pi
fig = plt.figure(1)
ax = fig.add_subplot(1,1,1)
ax.set_title("odeint")
ax.plot(tt,(lsg[:,2]*deg), label= r"$q_3$")
ax.plot(tt, (lsg[:,0]*deg), label= r"$q_1$")
ax.plot(tt,(lsg[:,1]*deg), label= r"$q_2$")
ax.plot(tt,(lsg[:,3]*deg), label= r"$q_4$")
ax.set_xlabel("$t$ in s")
ax.set_ylabel("$\phi_2$")
ax.legend()
plt.show()
np.save('lsg_outfile',lsg2)

# <codecell>
# Regelabweichungen veranschaulicht
fig2 = plt.figure(2)
a1=fig2.add_subplot(2,2,1)
a1.set_title('e1')
a1.plot(tt,e1*deg)
a2=fig2.add_subplot(2,2,2)
a2.set_title('e1_d')
a2.plot(tt,e1_d*deg)
a3=fig2.add_subplot(2,2,3)
a3.set_title('e2')
a3.plot(tt,e2*deg)
a4=fig2.add_subplot(2,2,4)
a4.set_title('e2_d')
a4.plot(tt,e2_d*deg)
plt.show()

#fig3 = plt.figure(3)
#ay = fig3.add_subplot(1,1,1)
#ay.set_title("control variable")
#ay.plot(np.linspace(0,40,f.shape[0]),f[:,0])
#ay.set_xlabel("$t$ in s")
#ay.set_ylabel("$r$")
#ay.legend()
#plt.show()
# <codecell>

def control_effort(lsg2):
    f_contr=[0]
    f_feedforw=[0]
    k=1e7*np.ones((nr_aj,1))
    man=np.ones((nr_aj,1))
    ti=lsg2[:,0]
    
    for z in range(len(ti)):
        states=qq_traj(ti[z])
        states2=states[:4:]
        states2=states2.reshape((2,-1))
        #states2=states2[::,:nr_aj:]
        # Fehler
        e=states2-lsg2[z,1::2].reshape((2,-1))
        e=e.T
        if z==0:
            f_contr=np.array((k*e).dot(man))
            f_feedforw=np.array(FF_fnc(states[0],0,states[1],0,states[2],0,states[3],0,states[4],0,states[5],0))
            #f2_feedforw=[F21_fnc(states[0],0,states[1],0,states[2],0,states[3],0,states[4],0,states[5],0)]
        else:    
            f_contr = np.concatenate((f_contr,(k*e).dot(man)),axis=1)  
            f_feedforw = np.concatenate((f_feedforw,FF_fnc(states[0],0,states[1],0,states[2],0,states[3],0,states[4],0,states[5],0)),axis=1)
            #f2_feedforw = np.concatenate((f2_feedforw,[F21_fnc(states[0],0,states[1],0,states[2],0,states[3],0,states[4],0,states[5],0)]),axis=1)
        
    f_feedforw=np.array(f_feedforw)
    fig3 = plt.figure(3)
    ay = fig3.add_subplot(1,2,1)
    ay.set_title("Regler, Vorsteueuerung: q1")
    ay.plot(ti,f_contr[0,:], label= r"f_contr q1" )
    ay.plot(ti,f_feedforw[0,:], label= r"f_feedfor q1")
    ay.set_xlabel("$t$ in s")
    ay.set_ylabel("$M$")
    ay.legend()
    az = fig3.add_subplot(1,2,2)
    az.set_title("Regler, Vorsteueuerung: q2")
    az.plot(ti,f_contr[1,:], label= r"f_contr q2" )
    az.plot(ti,f_feedforw[1,:], label= r"f_feedfor q2")
    az.set_xlabel("$t$ in s")
    az.set_ylabel("$M$")        
    az.legend()
    plt.show()
        #f_contr,f_feedforw=control_effort(lsg2)
    return f_contr,f_feedforw
f_contr,f_feedforw=control_effort(lsg2)
