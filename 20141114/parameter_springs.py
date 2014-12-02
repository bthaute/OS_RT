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
    para_m = np.dot(np.matrix([[(1.6 + 0.95),(1.6 + 0.95)],[(1.1 + 0.6),(1.1 + 0.6)],[(0.8 + 0.55),(0.8 + 0.55)],[(0.5 + 0.4),(0.5 + 0.4)],[(0.28 + 0.2),(0.28 + 0.2)]]),500)
        
    # Längen der Arme in [m]
    para_a = np.dot(np.matrix([[9,9],[8,8],[7,7],[7,7],[6,6]]),0.5)
   
    # Schwerpunkte Annahme bei der halben Länge
    para_l = np.dot(para_a,0.5)
        
    # Höhen, Breiten und Dicken der Arme bei Annahme eines rechteckigen Hohträgers
    
    # Höhen in [m]
    para_h = np.dot(np.matrix([[0.5,0.5],[0.35,0.35],[0.3,0.3],[0.25,0.25],[0.2,0.2]]),0.5)
    
    # Breiten in [m]
    para_b = np.dot(np.matrix([[0.5,0.5],[0.35,0.35],[0.3,0.3],[0.25,0.25],[0.2,0.2]]),0.5)
    
    # Dicken in [m]
    para_di = [0.02,0.015,0.015,0.01,0.01]
    
    
    # Federkonstanten -> noch keine Ahnung
    para_k = [1e6,1e6,0,0,0]
    
#    Federkonstanten -> noch keine Ahnung
    para_d = [100,1e4,100,100,100]   
    
    # Trägheitsmomente Beachte: Bei uns wird momentan nur Trägheit in x benötigt in[kg*m^2]
    para_Jx = np.matrix([[0.1,0.1],[0.1,0.1],[0.1,0.1],[0.1,0.1],[0.1,0.1]])
    para_Jy = np.matrix([[0.1,0.1],[0.1,0.1],[0.1,0.1],[0.1,0.1],[0.1,0.1]])
    para_Jz = np.matrix([[0.1,0.1],[0.1,0.1],[0.1,0.1],[0.1,0.1],[0.1,0.1]])
    for index1 in range(0,5):
        for index2 in range(0,2):
            para_Jx[index1,index2] = Jx_fnc(para_h[index1,index2],para_b[index1,index2],para_a[index1,index2],para_di[index1],para_m[index1,index2])
            para_Jy[index1,index2] = Jx_fnc(para_h[index1,index2],para_b[index1,index2],para_a[index1,index2],para_di[index1],para_m[index1,index2])
            para_Jz[index1,index2] = Jx_fnc(para_h[index1,index2],para_b[index1,index2],para_a[index1,index2],para_di[index1],para_m[index1,index2])
    
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