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

# In[55]:

import sympy as sp
import numpy as np
import pickle

pfilepath = "para5dof.pcl"

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
    Jx2 = (M1*((H - 2*D)**2 + L**2))/12
    Jy2 = (M1*((B - 2*D)**2 + L**2))/12
    Jz2 = (M1*((B - 2*D)**2 + (H - 2*D)**2))/12
    
    Jx = Jx1 - Jx2
    Jy = Jy1 - Jy2
    Jz = Jz1 - Jz2
    Jx.simplify()
    Jy.simplify()
    Jz.simplify()
    Jx_fnc = sp.lambdify([H, B, L, D, M],Jx,'numpy')
    Jy_fnc = sp.lambdify([H, B, L, D, M],Jy,'numpy')
    Jz_fnc = sp.lambdify([H, B, L, D, M],Jz,'numpy')
    # Vorgabe 
    
    # Konstanten
    
    # g in [m/s^2]
    para_g =  9.80665
    
    # Massen der Arme inkl. Gelenke in [kg]
    para_m = np.dot([(1.6 + 0.95),(1.1 + 0.6),(0.8 + 0.55),(0.5 + 0.4),(0.28 + 0.2)],1000)
        
    # Längen der Arme in [m]
    para_l = [9,8,7,7,6]
   
    # Schwerpunkte Annahme bei der halben Länge
    para_a = np.dot(para_l,0.5)
        
    # Höhen, Breiten und Dicken der Arme bei Annahme eines rechteckigen Hohträgers
    
    # Höhen in [m]
    para_h = [0.5,0.35,0.3,0.25,0.2]
    
    # Breiten in [m]
    para_b = [0.5,0.35,0.3,0.25,0.2]
    
    # Dicken in [m]
    para_di = [0.02,0.015,0.015,0.01,0.01]
    
    
    # Federkonstanten -> noch keine Ahnung
    para_k = [0,0,0,0,0]
    
#    Federkonstanten -> noch keine Ahnung
    para_d = [100,100,100,100,100]   
    
    # Trägheitsmomente Beachte: Bei uns wird momentan nur Trägheit in x benötigt in[kg*m^2]
    para_Jx = [0,0,0,0,0]
    para_Jy = [0,0,0,0,0]
    para_Jz = [0,0,0,0,0]
    for index in range(0,5):
        para_Jx[index] = Jx_fnc(para_h[index],para_b[index],para_l[index],para_di[index],para_m[index])
        para_Jy[index] = Jx_fnc(para_h[index],para_b[index],para_l[index],para_di[index],para_m[index])
        para_Jz[index] = Jx_fnc(para_h[index],para_b[index],para_l[index],para_di[index],para_m[index])
    
    para_I = para_Jx
        
    print "saving parameter"
    pdict = {"para_g": para_g, "para_m": para_m, "para_l": para_l, "para_a": para_a, "para_k": para_k, "para_d": para_d, "para_I": para_I}

    with  open(pfilepath, "w") as pfile:
        pickle.dump(pdict, pfile)

else:  # flag_new_parameter == False

    with  open(pfilepath, "r") as pfile:
        pdict = pickle.load(pfile)

    para_g = pdict["para_g"]
    para_m = pdict["para_m"]
    para_l = pdict["para_l"]
    para_a = pdict["para_a"]
    para_k = pdict["para_k"]
    para_d = pdict["para_d"]
    para_I = pdict["para_I"]

print "read parameters"