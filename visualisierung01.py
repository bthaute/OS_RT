# -*- coding: utf-8 -*-

import pylab as pl # = matplotlib
import numpy as np
from numpy import pi, exp, r_

# <codecell>
#from main01_20141110 import lsg
lsg=np.load('lsg_outfile.npy')

# <codecell>

def live_animation():
    dt = 0.3  # in s # hat irgendwie keinen Einfluss -> Egal.
    interval=int(dt * 1000)  # interval in ms
    timer = fig.canvas.new_timer(interval=interval)

    timer.add_callback(update_plot, ax)
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
    CCframe(q1, q2, q3)

    # Ausgabe der aktuellen Zeit
    pl.text(0.06, 0.05, "t = %3.2fs" % t, transform = axes.transAxes)
    pl.axis([-3, 3, -3, 3])
    axes.figure.canvas.draw()




def CCframe(q1, q2, q3, xy = 0):
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
    pl.plot(r_[J0,].real, r_[J0,].imag, 'ks', ms = 8)
    pl.plot(r_[J0, J1].real, r_[J0, J1].imag, 'k-', lw=3)
    pl.plot(r_[J2, J1].real, r_[J2, J1].imag, 'ko-', lw=2)
    pl.plot(r_[J2, J3].real, r_[J2, J3].imag, 'ko-')
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
di = 2


l1 = 1.0
l2 = 1.0
l3 = 1.0



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
pl.axis([-3, 3, -3, 3])
my_axis = pl.axis() # zoom merken (wegen 'equal' nicht identisch mit Vorgaben)
ax.axis('off')


#L = 1000
L=500
#tt = np.linspace(0, 10, L)
#qq1 = np.linspace(0, 3*pi, L)
#qq2 = np.linspace(0, 5*pi, L)
#qq3 = np.sin(np.linspace(0, 2*pi, L)) * pi
tt=np.linspace(1, 500, L)
qq1=lsg[0:L,0]
qq2=lsg[0:L,1]
qq3=np.linspace(0, 0, L)


C.i = 100
#C.i = 50
#
#update_plot(ax)
#pl.show()

live_animation()
