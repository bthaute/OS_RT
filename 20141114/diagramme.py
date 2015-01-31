# -*- coding: utf-8 -*-
"""
Created on Tue Jan  6 12:57:20 2015

@author: richard
"""

import sympy as sp
import numpy as np
from numpy import pi, exp, r_

#import matplotlib.patheffects
import matplotlib.pyplot as pl
from matplotlib import rcParams
#import svgutils.transform as sg
#from pylab import *

anzahl_subplot = 2

lsg=np.load('lsg_outfile_unvollvoll.npy')
#lsg=np.load('Trajektorie_Bsp.npy')

tt=lsg[:,0]
#tt = np.linspace(0,30,1000)
#lsg = [np.e**(-np.log(10)/30*tt)*(np.sin(2*pi*tt)),1*np.e**(-np.log(10)/30*tt),-1*np.e**(-np.log(10)/30*tt)]

#if pl.fignum_exists(1):
#    pl.close()
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

rcParams['xtick.labelsize'] = '24'#15 # sets fontsize of the tick labels
rcParams['ytick.labelsize'] = '24'#15 # sets fontsize of the tick labels
rcParams['figure.figsize'] = '12,8' # sets figure size in inches
rcParams['axes.titlesize'] = '24' #18# sets fontsize of the axes title
rcParams['axes.labelsize'] = '24' #16# sets fontsize of the x any y labels
#==============================================================================

fig1 = pl.figure()
#fig1.canvas.set_window_title(r'Zustandsgrößen')
#fig1.suptitle(r'Zustandsgroessen', fontsize=20)




# subplot for the stable, but not optimal, control system
ax1 = fig1.add_subplot(anzahl_subplot,1,1)

# Figure Settings
#==============================================================================
#ax1.set_title(r'\bf{Gelenkwinkel}',
#                    loc = 'center')
# hide right and top spine
ax1.spines['right'].set_color('none')
ax1.spines['top'].set_color('none')
#ax1.set_ylim([-10, 10])
#ax1.set_xlim([0,4])
#==============================================================================

t_end=250
#ax1.plot(tt, res_bad[:,0], 'k')
#ax1.plot(tt[:t_end], lsg[:t_end,1]*180/pi, label= r'$\theta_{11}$')#, 'k')
#ax1.plot(tt[:t_end], lsg[:t_end,2]*180/pi, label= r'$\theta_{12}$')#, 'k')
ax1.plot(tt[:t_end], lsg[:t_end,1]*180/pi, label= r'$\theta_{11}$')#, 'k')
#ax1.plot(tt[:t_end], lsg[:t_end,4]*180/pi, label= r'$\theta_{22}$')#, 'k')
ax1.set_ylabel(r'Winkel in $^\circ$')
ax1.set_xlabel(r'Zeit in s')
ax1.legend()
pl.grid()

if anzahl_subplot>1:
    
    ax2 = fig1.add_subplot(2,1,2)
    #ax2.spines['right'].set_color('none')
    #ax2.spines['top'].set_color('none')
    ax2.plot(tt[:t_end], lsg[:t_end,3]*180/pi, label= r'$\theta_{21}$')#, 'k')
    ax2.set_ylabel(r'Winkel in $^\circ$')
    ax2.set_xlabel(r'Zeit in s')
    ax2.legend()
    pl.grid()
    #
    pl.show()
    #pl.close(fig1)
#pl.savefig('trajtest.png')

#
#
#l2 = 0.5
#l3 = 0.5
#l4 = 0.5
#
#t=200
#
## Grafik-Fenster einrichten:
#
#pl.rcParams['figure.subplot.bottom']=.1
#pl.rcParams['figure.subplot.left']=.05
#pl.rcParams['figure.subplot.top']=.98
#pl.rcParams['figure.subplot.right']=.98
#
#mm = 1./25.4 #mm to inch
#scale = 3
#fs = (85*mm*scale, 65*mm*scale)
#fig = pl.figure(figsize = fs)
#ax = fig.add_subplot(1, 1, 1)#
#pl.axis('equal')
#pl.axis([-3, 3, -3, 3])
#my_axis = pl.axis() # zoom merken (wegen 'equal' nicht identisch mit Vorgaben)
#ax.axis('off')
#
#q1=lsg[t,1]
#q2=lsg[t,2]
#q3=lsg[t,3]
#q4=lsg[t,4]
#
#J0= 0+0j 
#J1= J0 + l1*exp(1j*q1)
#J2= J1 + l2*exp(1j*(q1+q2))
#J3= J2 + l3*exp(1j*(q1+q2+q3))
#J4= J3 + l4*exp(1j*(q1+q2+q3+q4))
#pl.plot(r_[J0,].real, r_[J0,].imag, 'ks', ms = 8)
#pl.plot(r_[J0, J1].real, r_[J0, J1].imag, 'k-', lw=3)
#pl.plot(r_[J2, J1].real, r_[J2, J1].imag, 'ko-', lw=2)
#pl.plot(r_[J2, J3].real, r_[J2, J3].imag, 'ko-', lw=2)
#pl.plot(r_[J4, J3].real, r_[J4, J3].imag, 'ko-')
#pl.xticks= []
#pl.yticks= []
#pl.axis('equal')
#l1 = 0.5
#pl.grid()