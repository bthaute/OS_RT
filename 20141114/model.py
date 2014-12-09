# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

import sympy as sp
import numpy as np
import model_tools as mt
import symb_tools as st
import pickle
import cholesky as chol
from parameter_springs import para_g, para_m, para_a, para_k, para_d, para_I
from IPython import embed as IPS
# Number of the actuated joints
nr_aj = 2
# Number of the actuated joints -> uses .simplify()
nr_simplify = 0
# Do you want to create an new model?
flag_new_model_generation = False


pfilepath = "model_"+np.str(nr_aj)+"_actuatedJoints_"+np.str(nr_aj)+"_nonactuatedJoints.pcl"
if flag_new_model_generation:  

    print "create new model"
    # Variablen definieren
    t = sp.Symbol("t")
    q1 = sp.Function("q1")(t)
    qq = sp.Matrix([q1])
    FF = sp.Matrix([sp.symbols("F1")])
    for index in range(2,2*(nr_aj)+1):
        locals()["q"+np.str(index)] = sp.Function("q"+np.str(index))(t)
        qq = st.row_stack(qq,sp.Matrix([locals()["q"+np.str(index)]]))
        FF = st.row_stack(FF,sp.Matrix([sp.symbols("F"+np.str(index))]))
    #Parameter
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
        
    params_values = zip(mm,para_m) + zip(II,para_I) + zip(aa,para_a) + zip(kk,para_k) + zip(dd,para_d) + zip(g,para_g)

    # Geometrie (karthesische Koordinaten der Schwerpunkte)
    ex = sp.Matrix([1, 0])
    ey = sp.Matrix([0, 1])
    # Schwerpunkte 
    q_sum = q1
    s1 = (aa[0]/2)*mt.Rz(q1)*ex
    ss = sp.Matrix([s1])
    for index in range(2,2*nr_aj+1):
        locals()["s"+np.str(index)] =  locals()["s"+np.str(index-1)] + (aa[index-2]/2)*mt.Rz(q_sum)*ex + (aa[index-1]/2)*mt.Rz(q_sum + locals()["q"+np.str(index)])*ex
        locals()["s"+np.str(index)].simplify()
        ss = st.col_stack(ss,sp.Matrix([locals()["s"+np.str(index)]]))
        q_sum = q_sum + locals()["q"+np.str(index)]
        
    # kinetische und potentielle Energie
    T = 0
    V = 0
    q_dd_sum = 0
    for index in range(0,2*nr_aj):
        q_dd_sum = q_dd_sum + qq[index].diff(t)
        s_diff_temp =ss[:,index].diff(t)
        T = T + (II[index]*q_dd_sum**2 + mm[index]*(s_diff_temp.T*s_diff_temp)[0])/2
        V = V + (kk[index]*qq[index]**2)/2 + g[0]*mm[index]*ss[1,index]
        
    # Erzeuge Modell nach Lagrange
    print "create model by model_tools"
    mod1 = mt.generate_model(T, V, qq, FF)
    if (nr_aj <= nr_simplify):
        for index in range(0,(mod1.eq_list.shape[0])):
            print ("vereinfache eq list",index+1)
            mod1.eq_list[index].simplify();
    
    # Dissipative Kräfte einbeziehen
    for index in range(0,mod1.eq_list.shape[0]):
        mod1.eq_list[index]=mod1.eq_list[index]+dd[index]*mod1.qds[index]
    
    #Linearisiere
    subslist_lin = zip(mod1.qdds, sp.zeros(mod1.qdds.shape[0]))+zip(mod1.qds,\
    sp.zeros(mod1.qds.shape[0]))+zip(mod1.qs, qq_0)+zip(mod1.extforce_list, sp.zeros(mod1.extforce_list.shape[0]))
    vector_exp = mod1.qdds
    vector_exp = st.row_stack(vector_exp,mod1.qds)
    vector_exp = st.row_stack(vector_exp,mod1.qs)
    vector_exp = st.row_stack(vector_exp,mod1.extforce_list)
    mod1.eq_list_lin = mod1.eq_list.jacobian([mod1.qdds,mod1.qds,mod1.qs,mod1.extforce_list]).subs(subslist_lin)*vector_exp    
    if (nr_aj <= nr_simplify):
        for index in range(0,(mod1.eq_list_lin.shape[0])):
            print ("vereinfache eq_list_lin",index+1)
            mod1.eq_list_lin[index].simplify();
        
    print "solving lin equation system"
    MM_lin = mod1.eq_list_lin.jacobian(mod1.qdds)

    subslist = zip(mod1.qdds, sp.zeros(mod1.qdds.shape[0]))
    temp_lin = mod1.eq_list_lin.subs(subslist)
    MM_temp_lin_inv = chol.inv(MM_lin) 
    if (nr_aj <= nr_simplify):
        for index1 in range(0,MM_temp_lin_inv.shape[0]):
            for index2 in range(0,MM_temp_lin_inv.shape[1]):
                print ("vereinfache MM_inv",index1+1,index2+1)
                MM_temp_lin_inv[index1,index2].simplify()
    
    q_dd_expr_lin =  MM_temp_lin_inv*(-temp_lin) 
    if (nr_aj <= nr_simplify):
        for index in range(0,q_dd_expr_lin.shape[0]):
            print ("vereinfache sol",index+1)
            q_dd_expr_lin[index].simplify()
    A_ = q_dd_expr_lin.jacobian([mod1.qds,mod1.qs])
    mod1.A = st.col_stack(sp.zeros(q_dd_expr_lin.shape[0]),sp.eye(q_dd_expr_lin.shape[0]))
    mod1.A = st.row_stack(mod1.A,A_)
    B_ = q_dd_expr_lin.jacobian(mod1.extforce_list)
    mod1.B = st.row_stack(sp.zeros(mod1.qs.shape[0],mod1.extforce_list.shape[0]),B_)
    print "create state space description"
    A_para_temp = mod1.A.subs(params_values)
    mod1.A_para = A_para_temp
    # Achtung hier gibt es
#    mod1.A_para.simpify()
    B_para_temp = mod1.B.subs(params_values)
    mod1.B_para = B_para_temp
#    mod1.B_para.simpify()
    IPS()
    # Speichere Modell in Datei ab.
    print "saving model"
    pdict = {"mod1": mod1}

    with  open(pfilepath, "w") as pfile:
        pickle.dump(pdict, pfile)

else:  # flag_new_model_generation == False

    print "just open exsiting model"
    with  open(pfilepath, "r") as pfile:
        pdict = pickle.load(pfile)

    mod1 = pdict["mod1"]
   
print "read model"