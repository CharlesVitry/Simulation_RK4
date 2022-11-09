from IHM import *
from models import *
import numpy as np
import pandas as pd
import plotly.express as px
import math

def pleasework():
    N = 400

    facteur = [0.4, 10]
    I0, R0 = 1, 0
    S0 = N - I0 - R0
    X0 = S0, I0, R0 
    model_rk4 = rk4(SIR, X0, t, N, delta, facteur)
    return model_rk4

def PLEASEWORK():
    # facteur = [B, Nu, mu, lambd, alpha, tau] avec B, mu, lambda vecteur
    N = 400

    B = [0.6, 0.4, 0.7, 0.5]
    Nu = 7
    mu = [0.04, 0.0001, 0.05, 0.001]
    lambd = [10, 10, 10, 10]
    alpha = 100
    tau = 10

    facteur = [B, Nu, mu, lambd, alpha, tau]

    

    # On initialise le nombre d'infectés, de recovered et de sains
    Isn0, Ivn0, Isr0, Ivr0 = 1, 1, 1, 1
    I0 = Isn0 + Ivn0 + Isr0 + Ivr0

    Rn0, Rr0 = 0, 0
    R0 = Rn0 + Rr0

    Sn0, Sr0 = (N - I0 - R0)/2, (N - I0 - R0)/2
    S0 = Sn0 + Sr0

    Vn0 = Vr0 = 0
    V0 = Vn0 + Vr0

    Csn, Csr0, Cvn0, Cvr0 = 0, 0, 0, 0
    C0 = Csn + Csr0 + Cvn0 + Cvr0

    Dn0, Dr0 = 0, 0
    D0 = Dn0 + Dr0

    X0 = Sn0, Sr0, Vn0, Vr0, Csn, Csr0, Cvn0, Cvr0, Isn0, Ivn0, Isr0, Ivr0, Rn0, Rr0, Dn0, Dr0 
    #print('XO', X0)
    #print('facteur', facteur)
    model_rk4 = rk4(SIRCVD, X0, t, N, delta, facteur)

    return model_rk4


N = 400


tmax = 365
Nt = 365
t = np.linspace(0, tmax, Nt+1)

delta = 2


model_rka_simple = pleasework()
S, I, R = model_rka_simple.T



model_rka_dur = PLEASEWORK()
Sn, Sr, Vn, Vr, Csn, Csr, Cvn, Cvr, Isn, Ivn, Isr, Ivr, Rn, Rr, Dn, Dr = model_rka_dur.T

toprint = [Sn, Sr, Vn, Vr, Csn, Csr, Cvn, Cvr, Isn, Ivn, Isr, Ivr, Rn, Rr, Dn, Dr]



I= Isn + Ivn + Isr + Ivr
R = Rn + Rr
S = Sn + Sr
V = Vn + Vr
C = Csn + Csr + Cvn + Cvr
D = Dn + Dr


#df = pd.DataFrame ({"t" : t, "S": S, "V" : V, "C" : C, "I": I, "R": R, "D": D})

#fig = px.line(df, x="t", y=["S", "V", "C", "I", "R", "D"], 
#    title='Simulation du modèle SIRCVD'
#    )
#fig.show()
#print(df.to_markdown())

