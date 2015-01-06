# -*- coding: utf-8 -*-
"""
Created on Thu Apr 10 11:02:56 2014

@author: chris2
"""
import sympy as sp
import numpy as np
from scipy.integrate import odeint
import scipy.signal

#import matplotlib.patheffects
import matplotlib.pyplot as pl
from matplotlib import rcParams
#import svgutils.transform as sg
#from pylab import *

import symb_tools as st
import robust_poleplacement as opt

#==============================================================================
    # setting rcParams for figures
    #=========================
rcParams['legend.numpoints'] = 1    # sets number of markers
                                        # (in one row) in legend
rcParams['text.usetex'] = True # sets if LaTeX is used for figs
    # sets which fonts to use
rcParams['font.family'] = 'serif'
rcParams['font.serif'] = 'lmodern'
rcParams['font.sans-serif'] = 'lmodern'
rcParams['font.cursive'] = 'lmodern'
rcParams['font.fantasy'] = 'lmodern'
rcParams['font.monospace'] = 'lmodern'
    
rcParams['text.latex.unicode'] = True # sets unicode using
rcParams['legend.fontsize'] = 'large' # sets fontsize
rcParams['figure.facecolor'] = '1' # sets color behind subplots
rcParams['lines.markersize'] = '8' # sets size of markers
rcParams['lines.linewidth'] = '1.0' # sets linewidth of functions

rcParams['xtick.labelsize'] = '15' # sets fontsize of the tick labels
rcParams['ytick.labelsize'] = '15' # sets fontsize of the tick labels
rcParams['figure.figsize'] = '12,8' # sets figure size in inches
rcParams['axes.titlesize'] = '18' # sets fontsize of the axes title
rcParams['axes.labelsize'] = '16' # sets fontsize of the x any y labels
#==============================================================================


# SIMULATION OF THE OPTIMAL POLEPLACEMENT BY REINSCHKE

# example 8.8
n = 5
m = 3


A = np.array([[1, 1, 0, 0, 0],
              [0, 1, 1, 0, 1],
              [-1, 1, 1, 0, -1],
              [0, 0, 1, 1, 1],
              [1, -1, 0, 0, 2],])

B = np.array([[1, 2, 1, 1, 1],
              [-1, -1, -1, 0, 1],
              [-2, -1, -1, -1, 1]]).T
              
eigvals = [-5, -4, -3, -2, -1]

# step 1: Bk is designed (not robust)
# design process (see: chapter 8.5.3)
####################################

kk = sp.symbols("k1:6")
s = sp.symbols("s")
I = sp.eye(5)

Bk = sp.Matrix([[-1, -2, -1, -2, -1],
              [2, 1, 2, 1, 2],
              kk])
A2 = A+ B*Bk

Lambda = np.diag(eigvals)


d1 = (s*I-A2).det()
d2 = (s*I-Lambda).det()

eqns = sp.Matrix (st.coeffs(d1, s)) - sp.Matrix(st.coeffs(d2, s))

sol = sp.solve(eqns, kk)

A2_bad = st.to_np(A2.subs(sol))
Bk_bad = st.to_np(Bk.subs(sol))
foo, V = np.linalg.eig(A2_bad)

# step 2:   Bk_opt is designed with optimal poleplacement
#           see module: opt_polplatzierung
#########################################################

Bk_opt = opt.opt_place_MI(A, B, eigvals)
A2_opt = A + B*sp.Matrix(Bk_opt)

## step 3:  simulation of the state variables
#           with the two different control matrices
#############################################################

tt = np.linspace(0, 10, 500)

# Initializing symbols for x(t) and derivative(x(t), t)
x_s = sp.symbols('x1:%i'%(n+1))
xd_s = sp.symbols('xd1:%i'%(n+1))

xt_vec = sp.Matrix(x_s)
xdt_vec = sp.Matrix(xd_s)

# Solving the linear equation system
sol_bad = sp.solve(xdt_vec - A2_bad*xt_vec, xdt_vec)
sol_opt = sp.solve(xdt_vec - A2_opt*xt_vec, xdt_vec)

# DEBUGGING
if 0:
    for i in range(len(xd_s)):
        print xd_s[i] , ':', sol_bad[xd_s[i]]
        print xd_s[i] , ':', sol_opt[xd_s[i]]

# List of expressions for every xdt
xd_bad_expr = [ sol_bad[xd_s[i]] for i in range(len(xd_s)) ]
xd_opt_expr = [ sol_opt[xd_s[i]] for i in range(len(xd_s)) ]

# lambdify gives function for fast calculation of xd_expr
xd_bad_fnc = [  sp.lambdify(x_s, xd_bad_expr[i], 'numpy')
                for i in range(len(xd_bad_expr)) ]
xd_opt_fnc = [  sp.lambdify(x_s, xd_opt_expr[i], 'numpy')
                for i in range(len(xd_opt_expr)) ]


def rhs_bad(z, t):
    x1, x2, x3, x4, x5 = z # Entpacken
    
    xd1 = xd_bad_fnc[0](x1, x2, x3, x4, x5)
    xd2 = xd_bad_fnc[1](x1, x2, x3, x4, x5)
    xd3 = xd_bad_fnc[2](x1, x2, x3, x4, x5)
    xd4 = xd_bad_fnc[3](x1, x2, x3, x4, x5)
    xd5 = xd_bad_fnc[4](x1, x2, x3, x4, x5)

    
    return np.array([xd1, xd2, xd3, xd4, xd5])
    
def rhs_opt(z, t):
    x1, x2, x3, x4, x5 = z # Entpacken
    
    xd1 = xd_opt_fnc[0](x1, x2, x3, x4, x5)
    xd2 = xd_opt_fnc[1](x1, x2, x3, x4, x5)
    xd3 = xd_opt_fnc[2](x1, x2, x3, x4, x5)
    xd4 = xd_opt_fnc[3](x1, x2, x3, x4, x5)
    xd5 = xd_opt_fnc[4](x1, x2, x3, x4, x5)

    
    return np.array([xd1, xd2, xd3, xd4, xd5])

x0 = [1,1,1,1,1] # Anfangszustand

res_bad = odeint(rhs_bad, x0, tt)
res_opt = odeint(rhs_opt, x0, tt)

#if 0:
#    for i in range(100):
#        print res[i][0], "#", res[i][1], "#", res[i][2], "#", res[i][3]
#
#if 0:
#    for i in range(500):
#        print "res:", res[i][0], "##", "res_opt:", res_opt[i][0]
        
        
# step 4: Visualization of the state variables
##############################################

fig1 = pl.figure()
fig1.canvas.set_window_title('Zustandsgrößen')
#fig1.suptitle(r'Zustandsgroessen', fontsize=20)




# subplot for the stable, but not optimal, control system
ax1 = fig1.add_subplot(2,1,1)

# Figure Settings
#==============================================================================
ax1.set_title(r'\bf{LAG-Regler}',
                    loc = 'center')
# hide right and top spine
ax1.spines['right'].set_color('none')
ax1.spines['top'].set_color('none')
ax1.set_ylim([-10, 10])
ax1.set_xlim([0,4])
#==============================================================================


ax1.plot(tt, res_bad[:,0], 'k')
ax1.plot(tt, res_bad[:,1], 'k')
ax1.plot(tt, res_bad[:,2], 'k')
ax1.plot(tt, res_bad[:,3], 'k')
ax1.plot(tt, res_bad[:,4], 'k')
pl.ylabel(u'Zustandsgrößen')


pl.grid()

ax2 = fig1.add_subplot(2,1,2)
# Figure Settings
#==============================================================================
ax2.set_title(r'\bf{Parameterrobuster Regler}',
                    loc = 'center')
# hide right and top spine
ax2.spines['right'].set_color('none')
ax2.spines['top'].set_color('none')
ax2.set_ylim([-1, 1])
ax2.set_xlim([0,4])
#==============================================================================

ax2.plot(tt, res_opt[:,0], 'k')
ax2.plot(tt, res_opt[:,1], 'k')
ax2.plot(tt, res_opt[:,2], 'k')
ax2.plot(tt, res_opt[:,3], 'k')
ax2.plot(tt, res_opt[:,4], 'k')
ax2.set_xlabel(r'Zeit[s]')
pl.ylabel(u'Zustandsgrößen')

pl.grid()

# saving the plot
pl.savefig('xsim_opt_place_MI.svg')

# step 5: simulation of the steering signals
#         with the two different control matrices
#############################################################

u_s = sp.symbols('u1:%i'%(m+1))

ut_vec = sp.Matrix(u_s)

sol_bad = sp.solve(ut_vec - Bk_bad*xt_vec, ut_vec)
sol_opt = sp.solve(ut_vec - Bk_opt*xt_vec, ut_vec)

u_bad_expr = [ sol_bad[u_s[i]] for i in range(len(u_s)) ]
u_opt_expr = [ sol_opt[u_s[i]] for i in range(len(u_s)) ]

u_bad_fnc = [  sp.lambdify(x_s, u_bad_expr[i], 'numpy')
                for i in range(len(u_bad_expr)) ]
u_opt_fnc = [  sp.lambdify(x_s, u_opt_expr[i], 'numpy')
                for i in range(len(u_opt_expr)) ]
                    
def u_rhs_bad(res_bad, t):
    x1, x2, x3, x4, x5 = res_bad # Entpacken
    
    u1 = u_bad_fnc[0](x1, x2, x3, x4, x5)
    u2 = u_bad_fnc[1](x1, x2, x3, x4, x5)
    u3 = u_bad_fnc[2](x1, x2, x3, x4, x5)
    
    return np.array([u1, u2, u3])
    
def u_rhs_opt(res_opt, t):
    x1, x2, x3, x4, x5 = res_opt # Entpacken
    
    u1 = u_opt_fnc[0](x1, x2, x3, x4, x5)
    u2 = u_opt_fnc[1](x1, x2, x3, x4, x5)
    u3 = u_opt_fnc[2](x1, x2, x3, x4, x5)
    
    return np.array([u1, u2, u3])
    
resu_bad = []
for i in range(len(tt)):
    resu_bad.append(u_rhs_bad(res_bad[i], i))
    
resu_opt = []
for i in range(len(tt)):
    resu_opt.append(u_rhs_opt(res_opt[i], i))


resu_bad = st.to_np(resu_bad)
resu_opt = st.to_np(resu_opt)

# step 6: Visualization of the steering variables
##############################################

fig2 = pl.figure()
fig2.canvas.set_window_title('Steuersignale')
#fig2.suptitle(r'Steuersignale', fontsize=20)



# subplot for the stable, but not optimal, control system
ax3 = fig2.add_subplot(2,1,1)

# Figure Settings
#==============================================================================
ax3.set_title(r'\bf{LAG-Regler}',
                    loc = 'center')
# hide right and top spine
ax3.spines['right'].set_color('none')
ax3.spines['top'].set_color('none')
ax3.set_ylim([-20, 40])
ax3.set_xlim([0,4])
pl.ylabel(r'Steuersignale')
#==============================================================================


ax3.plot(tt, resu_bad[:,0], 'k')
ax3.plot(tt, resu_bad[:,1], 'k')
ax3.plot(tt, resu_bad[:,2], 'k')

pl.grid()

ax4 = fig2.add_subplot(2,1,2)
# Figure Settings
#==============================================================================
ax4.set_title(r'\bf{Parameterrobuster Regler}',
                    loc = 'center')
# hide right and top spine
ax4.spines['right'].set_color('none')
ax4.spines['top'].set_color('none')
ax4.set_ylim([-4, 2])
ax4.set_xlim([0,4])
pl.ylabel(u'Steuersignale')
#==============================================================================

ax4.plot(tt, resu_opt[:,0], 'k')
ax4.plot(tt, resu_opt[:,1], 'k')
ax4.plot(tt, resu_opt[:,2], 'k')
ax4.set_xlabel(r'Zeit[s]')

pl.grid()

# saving the plot
pl.savefig('usim_opt_place_MI.svg')

pl.show()
