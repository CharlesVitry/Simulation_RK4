import streamlit as st
import pandas as pd
import numpy as np
import altair as alt
import plotly.express as px
from PIL import Image

from models import *

def IHM():
    
    st.set_page_config(layout="wide")                                           # Pour avoir une page large
    style()                                                                     # Pour changer le style de la page

    N = 500

    # Sidebar
    st.sidebar.title("Menu")
    choix_page = st.sidebar.selectbox("Choisissez une page :", 
    ["Présentation",
    "Modèle SIR classique",
    "Modèle SIRCVD"])

    st.sidebar.markdown("## Paramètres")

    N = st.sidebar.slider("Taille de la population :", min_value=50, max_value=1000, value=500, step=50)

    tmax = st.sidebar.slider("Durée de la simulation (en jours) :", min_value=30, max_value=365, value=100, step=5)

    if choix_page == "Présentation" :
        load_page_accueil()
    elif choix_page == "Modèle SIR classique" :
        load_page_sir_classique(N, tmax)
    elif choix_page == "Modèle SIRCVD" :
        st.sidebar.write("Scénario :")
        geste_barriere = st.sidebar.checkbox("Gestes barrières")
        confinement = st.sidebar.checkbox("Confinement") 
        vaccination =  st.sidebar.checkbox("Vaccination")
        load_page_sir_modif(N, tmax, geste_barriere, confinement, vaccination)


def model_rk4_SIR(N, tmax):
    I0, R0 = 1, 0                                                               # On initialise le nombre d'infectés et de recovered	
    S0 = N - I0 - R0                                                            # On initialise le nombre de sains                       

    beta, lambd = 0.4, 10                                                       # On initialise les paramètres du modèle               
    facteur = [beta, lambd]                                                     # On met les paramètres dans un vecteur                      
    delta = 2                                                                   # On initialise le pas de temps                        

    Nt = tmax                                                                   # On initialise le nombre de points de la grille                       
    t = np.linspace(0, tmax, Nt+1)                                              # On initialise la grille de temps        

    X0 = S0, I0, R0                                                             # On initialise le vecteur d'état initial            

    model_rk4 = rk4(SIR, X0, t, N, delta, facteur)                              # On résout le système d'équations différentielles avec la méthode rk4      

    return t, model_rk4

def model_rk4_SIRCVD(N, tmax, geste_barriere, confinement, vaccination):
    # facteur = [B, Nu, mu, lambd, alpha, tau] avec B, mu, lambda vecteur
    B = [0.5, 0.6, 0.2, 0.3]
    Nu = 7
    mu = [0.0001, 0.05, 0.005, 0.01]
    lambd = [10, 10, 10, 10]
    alpha = 100
    tau = 10

    if geste_barriere:                                                          # Si on a coché la case geste barrière > on diminue le taux de transmission                   
        B = [B[i]*(1 - 0.2) for i in range(len(B))]
    else :                                                                      # Sinon on garde les taux de transmission initiaux              
        #B = B
        pass

    facteur = [B, Nu, mu, lambd, alpha, tau]                                    # On met les paramètres dans un vecteur

    Isn0, Ivn0, Isr0, Ivr0 = 1, 1, 1, 1                                         # On initialise le nombre d'infectés    
    I0 = Isn0 + Ivn0 + Isr0 + Ivr0                                         

    Rn0, Rr0 = 0, 0                                                             # On initialise le nombre de recovered               
    R0 = Rn0 + Rr0                                         

    Sn0, Sr0 = (N - I0 - R0)/2, (N - I0 - R0)/2                                 # On initialise le nombre de sains
    S0 = Sn0 + Sr0                                      

    Vn0 = Vr0 = 0                                                               # On initialise le nombre de vaccinés          
    V0 = Vn0 + Vr0                         

    Csn, Csr0, Cvn0, Cvr0 = 0, 0, 0, 0                                          # On initialise le nombre de cas confirmés   
    C0 = Csn + Csr0 + Cvn0 + Cvr0           

    Dn0, Dr0 = 0, 0                                                             # On initialise le nombre de décès 
    D0 = Dn0 + Dr0

    X0 = Sn0, Sr0, Vn0, Vr0, Csn, Csr0, Cvn0, Cvr0, Isn0, Ivn0, Isr0, Ivr0, Rn0, Rr0, Dn0, Dr0 # On initialise le vecteur d'état initial
    
    Nt = tmax                                                                   # On initialise le nombre de points de la grille    
    t = np.linspace(0, tmax, Nt+1)                                              # On initialise la grille de temps

    delta = 2                                                                   # On initialise le pas de temps                        

    model_rk4 = rk4(SIRCVD, X0, t, N, delta, facteur, geste_barriere, confinement, vaccination) # On résout le système d'équations différentielles avec la méthode rk4

    return t,model_rk4

def affichage(t, S, I, R):
    df = pd.DataFrame ({"t" : t, "S": S, "I": I, "R": R})

    fig = px.line(df, x="t", y=["S", "I", "R"], 
    title='Simulation du modèle SIR'
    )
    
    st.plotly_chart(fig, use_container_width=True)
    
def style():
    st.markdown(   
        f"""   
        <style>
        .stApp {{
             background-color : white
             text-justify: auto
         }}
         </style>
         """,
         unsafe_allow_html=True
     )

def load_page_accueil():
    st.header("Projet Simulation IMA 5")
    st.markdown("> Clovis Delêtre & Charles Vitry")
    
    st.write("L'ambition de ce projet est de simuler par une approche Systèmes Dynamiques l'évolution d'une épidémie. "
             "Cette méthode se base sur les équations différentielles déterminites établies sur un modèle SIR donnée")
    
    st.write("Dans le cadre de ce projet nous avons utilisés deux modélès : ")
    st.write("> - Un SIR classique,")
    st.write("> - Un SIR que nous avons modifiés.")

    st.markdown(" là y aura une image askip (model1)")
    st.markdown(" là y aura une image askip (model2)")

def espace_entre_parties():
    st.write("")
    st.markdown("""<hr style="height:2px;width:90%;border:none;color:#333;background-color:#333;" /> """, unsafe_allow_html=True) 
    st.write("")

def load_page_sir_classique(N, tmax):
    st.header("Simulation du modèle SIR")

    st.write("Le modèle SIR se base sur une notion de compartiments. "
             "Pour une population (de taille N donnée), on étudie la taille des trois sous populations au cours du temps (t). ")

    st.write("> - S(t) représente les personnes saines")
    st.write("> - I(t) représente aux personnes infectées")
    st.write("> - R(t) représente les personnes rétablies")

    espace_entre_parties()
    st.write("Dans notre cas on va considérer que :")
    st.write("> - La maladie / virus n'est pas mortel,")
    st.write("> - Les personnes rétablies sont immunisées,")
    st.write("> - Il n'y a pas de naissances pendant la période 't'.")

    espace_entre_parties()
    st.write("On peut alors définir les équations différentielles suivantes : ")
    st.write("> - dS/dt = -beta * S * I")
    st.write("> - dI/dt = beta * S * I - lambda * I")
    st.write("> - dR/dt = lambda * I")
    st.write("Avec : ")
    st.write("> - beta : taux d'infection")
    st.write("> - lambda : taux de guérison")

    st.write("A l'aide de la méthode Runge-Kutta de degré 4 (RK4) on observer la simulation suivante :")

    t, model_rk4_SIR_value = model_rk4_SIR(N, tmax) 

    S, I, R = model_rk4_SIR_value.T
    affichage(t, S, I, R)


def load_page_sir_modif(N, tmax, geste_barriere, confinement, vaccination): 

    t, model_rka_SIRCVD_value = model_rk4_SIRCVD(N, tmax, geste_barriere, confinement, vaccination)
    Sn, Sr, Vn, Vr, Csn, Csr, Cvn, Cvr, Isn, Ivn, Isr, Ivr, Rn, Rr, Dn, Dr = model_rka_SIRCVD_value.T

    df_global = pd.DataFrame({"t" : t, "Sn": Sn, "Sr": Sr, "Vn": Vn, "Vr": Vr, "Csn": Csn, "Csr": Csr, "Cvn": Cvn, "Cvr": Cvr, "Isn": Isn, "Ivn": Ivn, "Isr": Isr, "Ivr": Ivr, "Rn": Rn, "Rr": Rr, "Dn": Dn, "Dr": Dr})

    t = t.astype(int)
    I = (Isn + Ivn + Isr + Ivr).astype(int)
    R = (Rn + Rr).astype(int)
    S = (Sn + Sr).astype(int)
    V = (Vn + Vr).astype(int)
    C = (Csn + Csr + Cvn + Cvr).astype(int)
    D = (Dn + Dr).astype(int)

    st.header("Simulation du modèle SIR adapté : SIRCVD")

    st.markdown("<div align='center'><br>"
                "<img src='model_SIRCVD.png' width='500' height='500' alt='Model-SIRCVD image problem' border='0'>"
               ,unsafe_allow_html=True)

    espace_entre_parties()

    col1, col2 = st.columns(2)

    col1.write("Le modèle SIRCVD se base également sur une notion de compartiments. ")
    col1.write("> - S(t) représente les personnes saines")
    col1.write("> - V(t) représente les personnes vaccinées")
    col1.write("> - C(t) représente les personnes contaminées")
    col1.write("> - I(t) représente aux personnes infectées")
    col1.write("> - R(t) représente les personnes rétablies")
    col1.write("> - D(t) représente les personnes décédées")

    col2.write("Et les paramètres suivant :")
    col2.write("> - B le vecteur des taux de transmission")
    col2.write("> - α le taux de vaccination")
    col2.write("> - v le taux d'incubation")
    col2.write("> - λ le vecteur des taux de guérison")
    col2.write("> - µ le vecteur des taux de mortalité")
    col2.write("> - τ le taux de rétablissement")

    espace_entre_parties()
    df_global = pd.DataFrame ({"t" : t, "S": S, "V" : V, "C" : C, "I": I, "R": R, "D": D})
    fig_global = px.line(df_global, x="t", 
                    y=["S", "V", "C", "I", "R", "D"], 
                    title='Simulation du modèle SIRCVD',
                    color_discrete_map={
                        "S": "green",
                        "V": "lightgreen",
                        "C": "orange",
                        "I": "red",
                        "R": "darkgreen",
                        "D": "black"
                    },
                    #markers = True
    
    )    
    st.plotly_chart(fig_global, use_container_width=True)
    
    espace_entre_parties()
    st.subheader("Comparaison à risque / sans risque")
    col1, col2 = st.columns(2)
    
    df_global_sans_risque = pd.DataFrame ({"t" : t, "S": Sn, "V" : Vn, "C" : Csn + Cvn, "I": Isn + Ivn, "R": Rn, "D": Dn})
    fig_global_sans_risque = px.line(df_global_sans_risque, x="t", 
                    y=["S", "V", "C", "I", "R", "D"], 
                    color_discrete_map={
                        "S": "green",
                        "V": "lightgreen",
                        "C": "orange",
                        "I": "red",
                        "R": "darkgreen",
                        "D": "black"
                    },
    )    
    col1.plotly_chart(fig_global_sans_risque, use_container_width=True)
    col1.write(df_global_sans_risque)

    df_global_a_risque = pd.DataFrame ({"t" : t, "S": Sr, "V" : Vr, "C" : Csr + Cvr, "I": Isr + Ivr, "R": Rr, "D": Dr})
    fig_global_a_risque = px.line(df_global_a_risque, x="t", 
                    y=["S", "V", "C", "I", "R", "D"],
                    color_discrete_map={
                        "S": "green",
                        "V": "lightgreen",
                        "C": "orange",
                        "I": "red",
                        "R": "darkgreen",
                        "D": "black"
                    },
    )    
    col2.plotly_chart(fig_global_a_risque, use_container_width=True)
    col2.write(df_global_a_risque)
    
