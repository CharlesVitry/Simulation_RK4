import pandas as pd
import numpy as np
import plotly.express as px

def transform_model_sir(t, model):
    S, I, R = model.T
    df = pd.DataFrame ({"t" : t, "S": S, "I": I, "R": R})
    fig = px.line(df, x="t", y=["S", "I", "R"], 
    title='Simulation du modèle SIR'
    )
    return df, fig  

def transform_model_sircvd_global(t, model):
    Sn, Sr, Vn, Vr, Csn, Csr, Cvn, Cvr, Isn, Ivn, Isr, Ivr, Rn, Rr, Dn, Dr = model.T
    df_global = pd.DataFrame({"t" : t, "Sn": Sn, "Sr": Sr, "Vn": Vn, "Vr": Vr, "Csn": Csn, "Csr": Csr, "Cvn": Cvn, "Cvr": Cvr, "Isn": Isn, "Ivn": Ivn, "Isr": Isr, "Ivr": Ivr, "Rn": Rn, "Rr": Rr, "Dn": Dn, "Dr": Dr})

    I = Isn + Ivn + Isr + Ivr
    R = Rn + Rr
    D = Dn + Dr
    S = Sn + Sr 
    V = Vn + Vr
    C = Csn + Csr + Cvn + Cvr

    df_global_merge = pd.DataFrame ({"t" : t, "S": S, "I": I, "R": R, "D": D, "V": V, "C": C})
    fig_global_merge = px.line(df_global_merge, x="t", y=["S", "I", "R", "D", "V", "C"],
    title='Simulation du modèle SIRCVD global'
    )
    fig_global_merge.update_traces(line=dict(color='green', width=2), selector=dict(name='S'))
    fig_global_merge.update_traces(line=dict(color='red', width=2), selector=dict(name='I'))
    fig_global_merge.update_traces(line=dict(color='lightblue', width=2), selector=dict(name='R'))
    fig_global_merge.update_traces(line=dict(color='black', width=2), selector=dict(name='D'))
    fig_global_merge.update_traces(line=dict(color='purple', width=2), selector=dict(name='V'))
    fig_global_merge.update_traces(line=dict(color='orange', width=2), selector=dict(name='C'))


    df_global_sans_risque = pd.DataFrame ({"t" : t, "S": Sn, "V" : Vn, "C" : Csn + Cvn, "I": Isn + Ivn, "R": Rn, "D": Dn})
    fig_global_sans_risque = px.line(df_global_sans_risque, x="t", y=["S", "I", "R", "D", "V", "C"],
    title='Simulation du modèle SIRCVD évolution des personnes sans risque'
    )
    fig_global_sans_risque.update_traces(line=dict(color='green', width=2), selector=dict(name='S'))
    fig_global_sans_risque.update_traces(line=dict(color='red', width=2), selector=dict(name='I'))
    fig_global_sans_risque.update_traces(line=dict(color='lightblue', width=2), selector=dict(name='R'))
    fig_global_sans_risque.update_traces(line=dict(color='black', width=2), selector=dict(name='D'))
    fig_global_sans_risque.update_traces(line=dict(color='purple', width=2), selector=dict(name='V'))
    fig_global_sans_risque.update_traces(line=dict(color='orange', width=2), selector=dict(name='C'))


    df_global_avec_risque = pd.DataFrame ({"t" : t, "S": Sr, "V" : Vr, "C" : Csr + Cvr, "I": Isr + Ivr, "R": Rr, "D": Dr})
    fig_global_avec_risque = px.line(df_global_avec_risque, x="t", y=["S", "I", "R", "D", "V", "C"],
    title='Simulation du modèle SIRCVD évolution des personnes avec risque'
    )
    fig_global_avec_risque.update_traces(line=dict(color='green', width=2), selector=dict(name='S'))
    fig_global_avec_risque.update_traces(line=dict(color='red', width=2), selector=dict(name='I'))
    fig_global_avec_risque.update_traces(line=dict(color='lightblue', width=2), selector=dict(name='R'))
    fig_global_avec_risque.update_traces(line=dict(color='black', width=2), selector=dict(name='D'))
    fig_global_avec_risque.update_traces(line=dict(color='purple', width=2), selector=dict(name='V'))
    fig_global_avec_risque.update_traces(line=dict(color='orange', width=2), selector=dict(name='C'))

    return df_global, df_global_merge, fig_global_merge, df_global_sans_risque, fig_global_sans_risque, df_global_avec_risque, fig_global_avec_risque

def transform_model_sircvd_echange(t, model):
    Sn, Sr, Vn, Vr, Csn, Csr, Cvn, Cvr, Isn, Ivn, Isr, Ivr, Rn, Rr, Dn, Dr, Sn_P, Sr_P, Vn_P, Vr_P, Csn_P, Csr_P, Cvn_P, Cvr_P, Isn_P, Ivn_P, Isr_P, Ivr_P, Rn_P, Rr_P, Dn_P, Dr_P = model.T
    df_global = pd.DataFrame({"t" : t, "Sn": Sn, "Sr": Sr, "Vn": Vn, "Vr": Vr, "Csn": Csn, "Csr": Csr, "Cvn": Cvn, "Cvr": Cvr, "Isn": Isn, "Ivn": Ivn, "Isr": Isr, "Ivr": Ivr, "Rn": Rn, "Rr": Rr, "Dn": Dn, "Dr": Dr, "Sn_P": Sn_P, "Sr_P": Sr_P, "Vn_P": Vn_P, "Vr_P": Vr_P, "Csn_P": Csn_P, "Csr_P": Csr_P, "Cvn_P": Cvn_P, "Cvr_P": Cvr_P, "Isn_P": Isn_P, "Ivn_P": Ivn_P, "Isr_P": Isr_P, "Ivr_P": Ivr_P, "Rn_P": Rn_P, "Rr_P": Rr_P, "Dn_P": Dn_P, "Dr_P": Dr_P})

    I = Isn + Ivn + Isr + Ivr 
    R = Rn + Rr
    D = Dn + Dr
    S = Sn + Sr
    V = Vn + Vr
    C = Csn + Csr + Cvn + Cvr

    I_P = Isn_P + Ivn_P + Isr_P + Ivr_P
    R_P = Rn_P + Rr_P
    D_P = Dn_P + Dr_P
    S_P = Sn_P + Sr_P
    V_P = Vn_P + Vr_P
    C_P = Csn_P + Csr_P + Cvn_P + Cvr_P

    df_merge = pd.DataFrame ({"t" : t, "S": S, "I": I, "R": R, "D": D, "V": V, "C": C, "S_P": S_P, "I_P": I_P, "R_P": R_P, "D_P": D_P, "V_P": V_P, "C_P": C_P})    
    fig_merge = px.line(df_merge, x="t", y=["S", "I", "R", "D", "V", "C", "S_P", "I_P", "R_P", "D_P", "V_P", "C_P"],
    title='Simulation du modèle SIRCVD échange'
    )
    fig_merge.update_traces(line=dict(color='green', width=2), selector=dict(name='S'))
    fig_merge.update_traces(line=dict(color='red', width=2), selector=dict(name='I'))
    fig_merge.update_traces(line=dict(color='lightblue', width=2), selector=dict(name='R'))
    fig_merge.update_traces(line=dict(color='black', width=2), selector=dict(name='D'))
    fig_merge.update_traces(line=dict(color='purple', width=2), selector=dict(name='V'))
    fig_merge.update_traces(line=dict(color='orange', width=2), selector=dict(name='C'))
    fig_merge.update_traces(line=dict(color='green', width=2), selector=dict(name='S_P'))
    fig_merge.update_traces(line=dict(color='red', width=2), selector=dict(name='I_P'))
    fig_merge.update_traces(line=dict(color='lightblue', width=2), selector=dict(name='R_P'))
    fig_merge.update_traces(line=dict(color='black', width=2), selector=dict(name='D_P'))
    fig_merge.update_traces(line=dict(color='purple', width=2), selector=dict(name='V_P'))
    fig_merge.update_traces(line=dict(color='orange', width=2), selector=dict(name='C_P'))


    df_classique = pd.DataFrame ({"t" : t, "S": S, "I": I, "R": R, "D": D, "V": V, "C": C})
    df_P = pd.DataFrame ({"t" : t, "S": S_P, "I": I_P, "R": R_P, "D": D_P, "V": V_P, "C": C_P})

    fig_classique = px.line(df_classique, x="t", y=["S", "I", "R", "D", "V", "C"],
    title='Simulation du modèle SIRCVD échange pop initiale'
    )
    fig_classique.update_traces(line=dict(color='green', width=2), selector=dict(name='S'))
    fig_classique.update_traces(line=dict(color='red', width=2), selector=dict(name='I'))
    fig_classique.update_traces(line=dict(color='lightblue', width=2), selector=dict(name='R'))
    fig_classique.update_traces(line=dict(color='black', width=2), selector=dict(name='D'))
    fig_classique.update_traces(line=dict(color='purple', width=2), selector=dict(name='V'))
    fig_classique.update_traces(line=dict(color='orange', width=2), selector=dict(name='C'))

    fig_P = px.line(df_P, x="t", y=["S", "I", "R", "D", "V", "C"],
    title='Simulation du modèle SIRCVD échange pop étrangère'
    )
    fig_P.update_traces(line=dict(color='green', width=2), selector=dict(name='S'))
    fig_P.update_traces(line=dict(color='red', width=2), selector=dict(name='I'))
    fig_P.update_traces(line=dict(color='lightblue', width=2), selector=dict(name='R'))
    fig_P.update_traces(line=dict(color='black', width=2), selector=dict(name='D'))
    fig_P.update_traces(line=dict(color='purple', width=2), selector=dict(name='V'))
    fig_P.update_traces(line=dict(color='orange', width=2), selector=dict(name='C'))


    return df_global, df_merge, fig_merge, df_classique, fig_classique, df_P, fig_P
