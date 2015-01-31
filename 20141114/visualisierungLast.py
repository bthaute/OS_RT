# -*- coding: utf-8 -*-

import pylab as pl # = matplotlib
import numpy as np
from numpy import pi, exp, r_
import time

## <codecell>
#from main01_20141110 import lsg
lsg=np.load('lsg_outfile_Last.npy')
#lsg=np.load('traj_01.npy')

## <codecell>

def live_animation():
    dt = 0.3  # in s # hat irgendwie keinen Einfluss -> Egal.
    interval=int(dt * 1000)  # interval in ms
    timer = fig.canvas.new_timer(interval=interval)
    
    timer.add_callback(update_plot,ax)
    timer.start()
    
    pl.show()



def update_plot(axes):
    """
    Das Plot-Fenster mit neuen Daten ausstatten.
    Diese Funktion wird vom matplotlib-timer-Objekt aufgerufen.
    """
    axes.clear()

    i = C.i
    C.i += di  # globale Zählvariable erhöhen
    if C.i >= len(tt):
        time.sleep(2)
        C.i = 0

    t = tt[i]
    q1 = qq1[i]
    q2 = qq2[i]
    q3 = qq3[i]
    q4 = qq4[i]
    CCframe(q1, q2, q3, q4, t)

    # Ausgabe der aktuellen Zeit
    pl.text(0.06, 0.05, "t = %3.2fs" % t, transform = axes.transAxes)
    pl.axis([-1, 2, -1, 1.5])
    axes.figure.canvas.draw()




def CCframe(q1, q2, q3, q4, t, xy = 0):
    """
    creates a cartesian coordinate frame
    """

    # Postitions of the joints
    # Komplexe Zahlen für Geometrie in der Ebene nutzen
    # (Rotationsmatrizen gingen genauso)

    J0= 0+0j + xy # offset
    J1= J0 + l1*exp(1j*q1)
    J2= J1 + l2*exp(1j*(q1+q2))
    J3= J2 + l3*exp(1j*(q1+q2+q3))
    J4= J3 + l4*exp(1j*(q1+q2+q3+q4))
    F = J4 + 0.5*exp(1j*(pi/2))
    pfeil_rechts = J4 + 0.1*exp(1j*(pi/3))
    pfeil_links = J4 + 0.1*exp(1j*(2*pi/3))
    pl.plot(r_[J0,].real, r_[J0,].imag, 'ks', ms = 8)
    pl.plot(r_[J0, J1].real, r_[J0, J1].imag, 'k-', lw=3)
    pl.plot(r_[J2, J1].real, r_[J2, J1].imag, 'ko-', lw=2)
    pl.plot(r_[J2, J3].real, r_[J2, J3].imag, 'ko-', lw=2)
    pl.plot(r_[J4, J3].real, r_[J4, J3].imag, 'ko-')
    if (t%2)<=0.5:
        pl.plot(r_[J4, J4].real, r_[J4, F].imag, 'r-', lw=2)
        pl.plot(r_[J4, pfeil_rechts].real, r_[J4, pfeil_rechts].imag, 'r-', lw=2)
        pl.plot(r_[J4, pfeil_links].real, r_[J4, pfeil_links].imag, 'r-', lw=2)
        
    pl.xticks= []
    pl.yticks= []
    pl.axis('equal')

#
#    xmin = np.min(r_[J0, J1, J2, J3].real)
#    xmax = np.max(r_[J0, J1, J2, J3].real)
#
#    ymin = np.min(r_[J0, J1, J2, J3].imag)
#    ymax = np.max(r_[J0, J1, J2, J3].imag)



class Container:
    """
    Hilfslasse um global schreibend auf Variablen zugreifen zu können (als Attribute)
    """
    pass

C = Container()
C.i = 0 # globaler counter
di = 20


l1 = 0.5
l2 = 0.5
l3 = 0.5
l4 = 0.5



# Grafik-Fenster einrichten:

pl.rcParams['figure.subplot.bottom']=.1
pl.rcParams['figure.subplot.left']=.05
pl.rcParams['figure.subplot.top']=.98
pl.rcParams['figure.subplot.right']=.98

mm = 1./25.4 #mm to inch
scale = 3
fs = (85*mm*scale, 65*mm*scale)
fig = pl.figure(figsize = fs)
ax = fig.add_subplot(1, 1, 1)#
pl.axis('equal')
pl.axis([-1, 3, -1, 3])
my_axis = pl.axis() # zoom merken (wegen 'equal' nicht identisch mit Vorgaben)
ax.axis('off')


#L = 1000
#L=500
#tt = np.linspace(0, 10, L)
#qq1 = np.linspace(0, 3*pi, L)
#qq2 = np.linspace(0, 5*pi, L)
#qq3 = np.sin(np.linspace(0, 2*pi, L)) * pi
#tt=np.linspace(1, 500, L)

#tt = np.linspace(-5, 15, 1e3)
#qq1=lsg[:,0]
#qq2=lsg[:,1]
#qq3=tt*0

tt=lsg[:,0]
qq1 = lsg[:,1]
qq2 = 10 * lsg[:,2]
qq3 = lsg[:,3]
qq4 = 30 * lsg[:,4]

if 0:
    from IPython import embed as IPS
    IPS()
    pl.figure()
    pl.plot(tt, qq1)
    pl.plot(tt, qq2)
    pl.plot(tt, qq3)
    pl.plot(tt, qq4)
    
    pl.show()

    raise SystemExit


#C.i = 0
#C.i = 50

#update_plot(ax)
#pl.show()

live_animation()