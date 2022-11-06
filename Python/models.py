import pandas as pd
import numpy as np

def SIR(X0, t, N, facteur):
    beta, lambd = facteur
    S0, I0, R0 = X0
    dSdt = - beta * I0 * S0 / N
    dIdt = beta * I0 * S0 / N - I0 / lambd
    dRdt = I0 / lambd
    return np.array([dSdt, dIdt, dRdt])


def SIRCVD(X0, t, N, facteur):
    beta, Nu, mu, lambd, alpha, tau = facteur # facteur = [B, Nu, mu, lambd, alpha, tau] avec B et mu vecteur

    Sn0, Sr0, Vn0, Vr0, Csn, Csr0, Cvn0, Cvr0, Isn0, Ivn0, Isr0, Ivr0, Rn0, Rr0, Dn0, Dr0 = X0

    I = Isn0 + Ivn0 + Isr0 + Ivr0

    dSndt = - beta[0] * I * Sn0 / N  - Sn0 / alpha
    dSrdt = - beta[1] * I * Sr0 / N - Sr0 / alpha

    dVndt = Rn0 / Nu + Sn0 / alpha - beta[2] * I * Vn0 / N
    dVrdt = Rr0 / Nu + Sr0 / alpha - beta[3] * I * Vr0 / N

    dCsn_dt = beta[0] * I * Sn0 / N - Csn / tau
    dCsrdt = beta[1] * I * Sr0 / N - Csr0 / tau
    dCvndt = beta[2] * I * Vn0 / N - Cvn0 / tau
    dCvrdt = beta[3] * I * Vr0 / N - Cvr0 / tau

    dIsn_dt = Csn / tau - Isn0 / lambd[0] - Isn0 * mu[0]
    dIsrdt = Csr0 / tau - Isr0 / lambd[1] - Isr0 * mu[1]
    dIvndt = Cvn0 / tau - Ivn0 / lambd[2] - Ivn0 * mu[2]
    dIvrdt = Cvr0 / tau - Ivr0 / lambd[3] - Ivr0 * mu[3]

    dRndt = Isn0 / lambd[0] + Ivn0 / lambd[2] - Rn0 / Nu
    dRrdt = Isr0 / lambd[1] + Ivr0 / lambd[3] - Rr0 / Nu

    dDndt = Isn0 * mu[0] + Ivn0 * mu[2]
    dDrdt = Isr0 * mu[1] + Ivr0 * mu[3]
    return np.array([dSndt, dSrdt, dVndt, dVrdt, dCsn_dt, dCsrdt, dCvndt, dCvrdt, dIsn_dt, dIvndt, dIsrdt, dIvrdt, dRndt, dRrdt, dDndt, dDrdt])








def rk2(f, X0, t, N, dt, facteur):
    nt = len(t)
    x = np.zeros([nt, len(X0)])
    x[0] = X0
    for i in range(nt-1):
        k1 = dt*f(x[i], t[i], N, facteur)
        k2 = dt*f(x[i] + k1, t[i] + dt, N, facteur)
        x[i+1] = x[i] + (k1 + k2)/2
    return x

def rk4(f, X0, t, N, dt, facteur):
    nt = len(t)
    x = np.zeros([nt, len(X0)])
    x[0] = X0
    for i in range(nt-1):
        k1 = dt*f(x[i], t[i], N, facteur)
        k2 = dt*f(x[i] + k1/2, t[i] + dt/2, N, facteur)
        k3 = dt*f(x[i] + k2/2, t[i] + dt/2, N, facteur)
        k4 = dt*f(x[i] + k3, t[i] + dt, N, facteur)
        x[i+1] = x[i] + (k1 + 2*k2 + 2*k3 + k4)/6
        #if i < 5:
        #    print("k1", k1)
        #    print("k2", k2)
        #    print("k3", k3)
        #    print("k4", k4)
        #    print("xi+1", x[i+1])
        #    print("")
    return x