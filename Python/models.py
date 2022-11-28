import pandas as pd
import numpy as np
import math

def SIR(X0, t, N, facteur):
    beta, lambd = facteur
    S0, I0, R0 = X0
    dSdt = - beta * I0 * S0 / N
    dIdt = beta * I0 * S0 / N - I0 / lambd
    dRdt = I0 / lambd
    return np.array([dSdt, dIdt, dRdt])


def SIRCVD(X0, t, N, facteur, geste_barriere, confinement, vaccination):
    beta, Nu, mu, lambd, alpha, tau = facteur # facteur = [B, Nu, mu, lambd, alpha, tau] avec B et mu vecteur

    Sn0, Sr0, Vn0, Vr0, Csn0, Csr0, Cvn0, Cvr0, Isn0, Ivn0, Isr0, Ivr0, Rn0, Rr0, Dn0, Dr0 = X0

    if (confinement and (t > 15 and t < 30)):
        beta = [beta[i] * (1 - 0.95) for i in range(len(beta))]
    else :
        #beta = beta
        pass

    ValeurVaccination = [None , None]
    if (vaccination) :
        alpha = (math.sin(t) * 0.5 * 100 + 0.5)  * alpha
        ValeurVaccination = [Sn0 / alpha , Sr0 / alpha]
    else :
        ValeurVaccination = [0, 0]

    I = Isn0 + Ivn0 + Isr0 + Ivr0

    dSndt = - beta[0] * I * Sn0 / N  - ValeurVaccination[0]
    dSrdt = - beta[1] * I * Sr0 / N - ValeurVaccination[1]

    dVndt = Rn0 / Nu + ValeurVaccination[0] - beta[2] * I * Vn0 / N
    dVrdt = Rr0 / Nu + ValeurVaccination[1] - beta[3] * I * Vr0 / N

    dCsn_dt = beta[0] * I * Sn0 / N - Csn0 / tau
    dCsrdt = beta[1] * I * Sr0 / N - Csr0 / tau
    dCvndt = beta[2] * I * Vn0 / N - Cvn0 / tau
    dCvrdt = beta[3] * I * Vr0 / N - Cvr0 / tau

    dIsn_dt = Csn0 / tau - Isn0 / lambd[0] - Isn0 * mu[0]
    dIsrdt = Csr0 / tau - Isr0 / lambd[1] - Isr0 * mu[1]
    dIvndt = Cvn0 / tau - Ivn0 / lambd[2] - Ivn0 * mu[2]
    dIvrdt = Cvr0 / tau - Ivr0 / lambd[3] - Ivr0 * mu[3]

    dRndt = Isn0 / lambd[0] + Ivn0 / lambd[2] - Rn0 / Nu
    dRrdt = Isr0 / lambd[1] + Ivr0 / lambd[3] - Rr0 / Nu

    dDndt = Isn0 * mu[0] + Ivn0 * mu[2]
    dDrdt = Isr0 * mu[1] + Ivr0 * mu[3]
    return np.array([dSndt, dSrdt, dVndt, dVrdt, dCsn_dt, dCsrdt, dCvndt, dCvrdt, dIsn_dt, dIvndt, dIsrdt, dIvrdt, dRndt, dRrdt, dDndt, dDrdt])


def SIRCVD_echange(X0, t, N, facteur, geste_barriere, confinement, vaccination):
    beta, Nu, mu, lambd, alpha, tau, echange = facteur # facteur = [B, Nu, mu, lambd, alpha, tau] avec B et mu vecteur

    Sn0, Sr0, Vn0, Vr0, Csn0, Csr0, Cvn0, Cvr0, Isn0, Ivn0, Isr0, Ivr0, Rn0, Rr0, Dn0, Dr0, Sn_P0, Sr_P0, Vn_P0, Vr_P0, Csn_P0, Csr_P0, Cvn_P0, Cvr_P0, Isn_P0, Ivn_P0, Isr_P0, Ivr_P0, Rn_P0, Rr_P0, Dn_P0, Dr_P0 = X0

    N_I, N_P = N[0], N[1]

    if (confinement and (t > 30 and t < 40)):
        beta = [beta[i] * (1 - 0.95) for i in range(len(beta))]
    else :
        #beta = beta
        pass

    ValeurVaccination = [None , None, None, None]
    if (vaccination) :
        alpha = (math.sin(t) * 0.5 * 100 + 0.5)  * alpha
        ValeurVaccination = [Sn0 / alpha , Sr0 / alpha, Sn_P0 / alpha, Sr_P0 / alpha]
    else :
        ValeurVaccination = [0, 0, 0, 0]

    I = Isn0 + Ivn0 + Isr0 + Ivr0
    I_P = Isn_P0 + Ivn_P0 + Isr_P0 + Ivr_P0

    # first pop
    dSndt = - beta[0] * I * Sn0 / N_I  - ValeurVaccination[0]
    dSrdt = - beta[1] * I * Sr0 / N_I - ValeurVaccination[1]

    dVndt = Rn0 / Nu + ValeurVaccination[0] - beta[2] * I * Vn0 / N_I
    dVrdt = Rr0 / Nu + ValeurVaccination[1] - beta[3] * I * Vr0 / N_I

    dCsndt = beta[0] * I * Sn0 / N_I - Csn0 / tau + echange * (Csn_P0 - Csn0)
    dCsrdt = beta[1] * I * Sr0 / N_I - Csr0 / tau + echange * (Csr_P0 - Csr0)
    dCvndt = beta[2] * I * Vn0 / N_I - Cvn0 / tau + echange * (Cvn_P0 - Cvn0)
    dCvrdt = beta[3] * I * Vr0 / N_I - Cvr0 / tau + echange * (Cvr_P0 - Cvr0)

    dIsndt = Csn0 / tau - Isn0 / lambd[0] - Isn0 * mu[0]
    dIsrdt = Csr0 / tau - Isr0 / lambd[1] - Isr0 * mu[1]
    dIvndt = Cvn0 / tau - Ivn0 / lambd[2] - Ivn0 * mu[2]
    dIvrdt = Cvr0 / tau - Ivr0 / lambd[3] - Ivr0 * mu[3]

    dRndt = Isn0 / lambd[0] + Ivn0 / lambd[2] - Rn0 / Nu
    dRrdt = Isr0 / lambd[1] + Ivr0 / lambd[3] - Rr0 / Nu

    dDndt = Isn0 * mu[0] + Ivn0 * mu[2]
    dDrdt = Isr0 * mu[1] + Ivr0 * mu[3]

    # second pop
    dSn_Pdt = - beta[0] * I_P * Sn_P0 / N_P - ValeurVaccination[2]
    dSr_Pdt = - beta[1] * I_P * Sr_P0 / N_P - ValeurVaccination[3]

    dVn_Pdt = Rn_P0 / Nu + ValeurVaccination[2] - beta[2] * I_P * Vn_P0 / N_P
    dVr_Pdt = Rr_P0 / Nu + ValeurVaccination[3] - beta[3] * I_P * Vr_P0 / N_P

    dCsn_Pdt = beta[0] * I_P * Sn_P0 / N_P - Csn_P0 / tau + echange * (Csn0 - Csn_P0)
    dCsr_Pdt = beta[1] * I_P * Sr_P0 / N_P - Csr_P0 / tau + echange * (Csr0 - Csr_P0)
    dCvn_Pdt = beta[2] * I_P * Vn_P0 / N_P - Cvn_P0 / tau + echange * (Cvn0 - Cvn_P0)
    dCvr_Pdt = beta[3] * I_P * Vr_P0 / N_P - Cvr_P0 / tau + echange * (Cvr0 - Cvr_P0)

    dIsn_Pdt = Csn_P0 / tau - Isn_P0 / lambd[0] - Isn_P0 * mu[0]
    dIsr_Pdt = Csr_P0 / tau - Isr_P0 / lambd[1] - Isr_P0 * mu[1]
    dIvn_Pdt = Cvn_P0 / tau - Ivn_P0 / lambd[2] - Ivn_P0 * mu[2]
    dIvr_Pdt = Cvr_P0 / tau - Ivr_P0 / lambd[3] - Ivr_P0 * mu[3]

    dRn_Pdt = Isn_P0 / lambd[0] + Ivn_P0 / lambd[2] - Rn_P0 / Nu
    dRr_Pdt = Isr_P0 / lambd[1] + Ivr_P0 / lambd[3] - Rr_P0 / Nu

    dDn_Pdt = Isn_P0 * mu[0] + Ivn_P0 * mu[2]
    dDr_Pdt = Isr_P0 * mu[1] + Ivr_P0 * mu[3]


    return np.array([dSndt, dSrdt, dVndt, dVrdt, dCsndt, dCsrdt, dCvndt, dCvrdt, dIsndt, dIvndt, dIsrdt, dIvrdt, dRndt, dRrdt, dDndt, dDrdt, dSn_Pdt, dSr_Pdt, dVn_Pdt, dVr_Pdt, dCsn_Pdt, dCsr_Pdt, dCvn_Pdt, dCvr_Pdt, dIsn_Pdt, dIvn_Pdt, dIsr_Pdt, dIvr_Pdt, dRn_Pdt, dRr_Pdt, dDn_Pdt, dDr_Pdt])


def rk2(f, X0, t, N, dt, facteur):
    nt = len(t)
    x = np.zeros([nt, len(X0)])
    x[0] = X0
    for i in range(nt-1):
        k1 = dt*f(x[i], t[i], N, facteur)
        k2 = dt*f(x[i] + k1, t[i] + dt, N, facteur)
        x[i+1] = x[i] + (k1 + k2)/2
    return x

def rk4(f, X0, t, N, dt, facteur, geste_barriere, confinement, vaccination):
    nt = len(t)
    x = np.zeros([nt, len(X0)])
    x[0] = X0
    for i in range(nt-1):
        k1 = dt*f(x[i], t[i], N, facteur, geste_barriere, confinement, vaccination)
        k2 = dt*f(x[i] + k1/2, t[i] + dt/2, N, facteur, geste_barriere, confinement, vaccination)
        k3 = dt*f(x[i] + k2/2, t[i] + dt/2, N, facteur, geste_barriere, confinement, vaccination)
        k4 = dt*f(x[i] + k3, t[i] + dt, N, facteur, geste_barriere, confinement, vaccination)
        x[i+1] = x[i] + (k1 + 2*k2 + 2*k3 + k4)/6
        #if i < 5:
        #    print("k1", k1)
        #    print("k2", k2)
        #    print("k3", k3)
        #    print("k4", k4)
        #    print("xi+1", x[i+1])
        #    print("")
    return x