# coding: utf-8

# Vorgabe reale Werte Betonpumpe

# Berechnung Trägheitsmomente

# Hohlrechteckträger

# H...Höhe
# 
# B...Breite
# 
# L...Länge
# 
# D...Wandstärke
# 
# M...Masse

import sympy as sp
import numpy as np
import pickle

pfilepath = "para5dof_springs.pcl"

flag_new_parameter = True
if flag_new_parameter:
    
    parameter = sp.symbols("H, B, L, D, M")
    H, B, L, D, M = parameter
    
    V1 = H*B*L
    V2 = (B - 2*D)*(H - 2*D)*L
    M1 = M*(V1/(V1-V2))
    M2 = M*(V2/(V1-V2))
    
    Jx1 = (M1*(H**2 + L**2))/12
    Jy1 = (M1*(B**2 + L**2))/12
    Jz1 = (M1*(B**2 + H**2))/12
    Jx2 = (M2*((H - 2*D)**2 + L**2))/12
    Jy2 = (M2*((B - 2*D)**2 + L**2))/12
    Jz2 = (M2*((B - 2*D)**2 + (H - 2*D)**2))/12
    
    Jx = Jx1 - Jx2 + M*(L/2)**2
    Jy = Jy1 - Jy2 + M*(L/2)**2
    Jz = Jz1 - Jz2
    Jx.simplify()
    Jy.simplify()
    Jz.simplify()
    Jx_fnc = sp.lambdify([H, B, L, D, M],Jx,'numpy')
    Jy_fnc = sp.lambdify([H, B, L, D, M],Jy,'numpy')
    Jz_fnc = sp.lambdify([H, B, L, D, M],Jz,'numpy')
    
    parameter = sp.symbols("J, g, omega_0")
    J, g, omega_0 = parameter
    
    k=(J + M*(L))*(2*np.pi*omega_0)**2 - M*g*(L)
    k_fnc = sp.lambdify([J, M, L, g, omega_0],k,'numpy')
    c=-(2*J + M*L)*(np.log(0.5)/10)
    c_fnc = sp.lambdify([J, M, L],c,'numpy')
    # Vorgabe 
    
    # Konstanten
    
    # g in [m/s^2]
    para_g =  np.array([9.80665])
    
    # Massen der Arme inkl. Gelenke in [kg]
    para_m = np.dot(np.array([(1.6 + 0.95),(1.6 + 0.95),(1.1 + 0.6),(1.1 + 0.6),\
    (0.8 + 0.55),(0.8 + 0.55),(0.5 + 0.4),(0.5 + 0.4),(0.28 + 0.2),(0.28 + 0.2)]),500)
        
    # Längen der Arme in [m]
    para_a = np.dot(np.array([9,9,8,8,7,7,7,7,6,6]),0.5)
   
    # Höhen, Breiten und Dicken der Arme bei Annahme eines rechteckigen Hohträgers
    # Höhen in [m]
    para_h = np.array([0.5,0.5,0.35,0.35,0.3,0.3,0.25,0.25,0.2,0.2])
    
    # Breiten in [m]
    para_b = para_h
    
    # Dicken in [m]
    para_di = np.array([0.02,0.02,0.015,0.015,0.015,0.015,0.01,0.01,0.01,0.01])
    
    
     #Federkonstanten -> rad/[Nm] bei 100kg belastung soll es sich um 10° biegen
#    bend = 10*np.pi/180 
#    load = 100*para_g
#    para_k = np.array([0,1,0,1,0,1,0,1,0,1])*para_a*(load/bend)
#    para_k = np.array([0,2e7,0,1e7,0,0,0,0,0,0])
    # Dämpfung-> noch keine Ahnung
#    para_d = np.array([0,1,0,1,0,1,0,1,0,1])* 1e1  
    
    # Trägheitsmomente Beachte: Bei uns wird momentan nur Trägheit in x benötigt in[kg*m^2]
    para_Jx = np.ones(10)
    para_Jy = np.ones(10)
    para_Jz = np.ones(10)
    for index in range(0,para_Jx.shape[0]):
            para_Jx[index] = Jx_fnc(para_h[index],para_b[index],para_a[index],para_di[index],para_m[index])
            para_Jy[index] = Jx_fnc(para_h[index],para_b[index],para_a[index],para_di[index],para_m[index])
            para_Jz[index] = Jx_fnc(para_h[index],para_b[index],para_a[index],para_di[index],para_m[index])
    
    para_I = para_Jx
    
    # Federkonstanten -> rad/[Nm] bei 1Hz Schwingung
    para_k = np.zeros(10,)
    m_temp = 0
    a_temp = 0
    I_temp = 0
    omega_0 = np.array([1,1,1.2,1.2,1.7,1.7,2.5,2.5,5,5])
    for index in range(9,-1,-1): 
        m_temp += para_m[index]
        a_temp += para_a[index]
        I_temp += para_I[index]
        for index2 in range(8,index-1,-1):
            I_temp += para_m[index2]*para_a[index]
        if (((index+1)%2)==0):
            para_k[index] = k_fnc(I_temp, m_temp, a_temp-(para_a[index]/2), para_g, omega_0[index])
      
    # Dämpfung-> nach 10s auf 10% gedämpft
    para_d = np.zeros(10,)
    for index in range(9,-1,-1): 
        m_temp += para_m[index]
        a_temp += para_a[index]
        I_temp += para_I[index]
        for index2 in range(8,index-1,-1):
            I_temp += para_m[index2]*para_a[index]
        para_d[index] = c_fnc(I_temp, m_temp, a_temp-(para_a[index]/2))
        
        
    print "saving parameter"
    pdict = {"para_g": para_g, "para_m": para_m,  "para_a": para_a, "para_k": para_k, "para_d": para_d, "para_I": para_I}

    with  open(pfilepath, "w") as pfile:
        pickle.dump(pdict, pfile)

else:  # flag_new_parameter == False

    with  open(pfilepath, "r") as pfile:
        pdict = pickle.load(pfile)

    para_g = pdict["para_g"]
    para_m = pdict["para_m"]
    para_a = pdict["para_a"]
    para_k = pdict["para_k"]
    para_d = pdict["para_d"]
    para_I = pdict["para_I"]

print "read parameters"