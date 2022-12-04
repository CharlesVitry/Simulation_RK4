import streamlit as st
import pandas as pd
import numpy as np
import altair as alt
import plotly.express as px
import matplotlib.pyplot as plt
import math

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
    "Modèle SIRCVD",
    "Modèle SIRCVD avec échange"])

    radio_markdown_gestes_barrieres = '''
        Reduit le taux de transmission par **0.4** par instant t. 
        '''.strip()
    radio_markdown_confinement = '''
        Reduit le taux de transmission par **0.05** entre les instants 15 et 30. 
        '''.strip()
    radio_markdown_vaccination = '''
        Modifie le taux de vaccination selon une fonction sinusoïdale.  
        '''.strip()

    if choix_page != "Présentation":
        st.sidebar.markdown("## Paramètres")

    if choix_page == "Présentation" :
        load_page_accueil()

    elif choix_page == "Modèle SIR classique" :
        N = st.sidebar.slider("Taille de la population :", min_value=50, max_value=1000, value=500, step=50)
        tmax = st.sidebar.slider("Durée de la simulation (en jours) :", min_value=30, max_value=365, value=100, step=5)
        geste_barriere = st.sidebar.checkbox("Gestes barrières", help=radio_markdown_gestes_barrieres)
        confinement = st.sidebar.checkbox("Confinement", help=radio_markdown_confinement) 
        
        load_page_sir_classique(N, tmax, geste_barriere, confinement, None)

    elif choix_page == "Modèle SIRCVD" :
        N = st.sidebar.slider("Taille de la population :", min_value=50, max_value=1000, value=500, step=50)
        tmax = st.sidebar.slider("Durée de la simulation (en jours) :", min_value=30, max_value=365, value=100, step=5)
        
        st.sidebar.write("Scénario :")
        geste_barriere = st.sidebar.checkbox("Gestes barrières", help=radio_markdown_gestes_barrieres)
        confinement = st.sidebar.checkbox("Confinement", help=radio_markdown_confinement)
        vaccination =  st.sidebar.checkbox("Vaccination", help=radio_markdown_vaccination)
        load_page_sir_modif(N, tmax, geste_barriere, confinement, vaccination)
    
    elif choix_page == "Modèle SIRCVD avec échange" :
        N = st.sidebar.slider("Taille de la population :", min_value=50, max_value=1000, value=500, step=50)
        N_P = st.sidebar.slider("Taille de la population étrangère :", min_value=50, max_value=1000, value=500, step=50)
        tmax = st.sidebar.slider("Durée de la simulation (en jours) :", min_value=30, max_value=365, value=100, step=5)
        
        N = [N, N_P]
        st.sidebar.write("Scénario :")
        geste_barriere = st.sidebar.checkbox("Gestes barrières", help=radio_markdown_gestes_barrieres)
        confinement = st.sidebar.checkbox("Confinement", help=radio_markdown_confinement)
        vaccination =  st.sidebar.checkbox("Vaccination", help=radio_markdown_vaccination)
        load_page_sir_modif_echange(N, tmax, geste_barriere, confinement, vaccination)

def model_rk4_SIR(N, tmax, geste_barriere, confinement, vaccination):
    I0, R0 = 1, 0                                                               # On initialise le nombre d'infectés et de recovered	
    S0 = N - I0 - R0                                                            # On initialise le nombre de sains                       

    B, lambd = 0.4, 10                                                       # On initialise les paramètres du modèle                     
    delta = 2                                                                   # On initialise le pas de temps                        

    if geste_barriere:                                                          # Si on a coché la case geste barrière > on diminue le taux de transmission                   
        B = B*(1 - 0.4) 
    else :                                                                      # Sinon on garde les taux de transmission initiaux              
        #B = B
        pass

    Nt = tmax                                                                   # On initialise le nombre de points de la grille                       
    t = np.linspace(0, tmax, Nt+1)                                              # On initialise la grille de temps        

    X0 = S0, I0, R0                                                             # On initialise le vecteur d'état initial            
                  
    facteur = [B, lambd]                                                        # On met les paramètres dans un vecteur  
    model_rk4 = rk4(SIR, X0, t, N, delta, facteur, geste_barriere, confinement, vaccination)                              # On résout le système d'équations différentielles avec la méthode rk4      

    return t, model_rk4

def model_rk4_SIRCVD(N, tmax, geste_barriere, confinement, vaccination):
    # facteur = [B, Nu, mu, lambd, alpha, tau] avec B, mu, lambda vecteur
    # [Sn, Sr, Vn, Vr]
    B = [0.5, 0.6, 0.2, 0.3]
    Nu = 5
    mu = [0.002, 0.003, 0.001, 0.001]
    #mu = [0, 0, 0, 0]
    eta = 0.5
    lambd = [10, 10, 10, 10]
    alpha = 500
    tau = 10

    if geste_barriere:                                                          # Si on a coché la case geste barrière > on diminue le taux de transmission                   
        B = [B[i]*(1 - 0.4) for i in range(len(B))]
    else :                                                                      # Sinon on garde les taux de transmission initiaux              
        #B = B
        pass

    facteur = [B, Nu, mu, eta, lambd, alpha, tau]                                    # On met les paramètres dans un vecteur
    
    Nt = tmax                                                                   # On initialise le nombre de points de la grille    
    t = np.linspace(0, tmax, Nt+1)                                              # On initialise la grille de temps

    delta = 2                                                                   # On initialise le pas de temps                        

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


    model_rk4 = rk4(SIRCVD, X0, t, N, delta, facteur, geste_barriere, confinement, vaccination) # On résout le système d'équations différentielles avec la méthode rk4
    
    
    return t,model_rk4


def model_rk4_SIRCVD_echange(N, tmax, geste_barriere, confinement, vaccination):

    # [Sn, Sr, Vn, Vr]
    B = [0.5, 0.6, 0.2, 0.3]
    Nu = 5
    mu = [0.002, 0.003, 0.001, 0.001]
    #mu = [0, 0, 0, 0]
    eta = 0.5
    lambd = [10, 10, 10, 10]
    alpha = 500
    tau = 10
    echange = 0.2

    if geste_barriere:                                                          # Si on a coché la case geste barrière > on diminue le taux de transmission                   
        B = [B[i]*(1 - 0.2) for i in range(len(B))]
    else :                                                                      # Sinon on garde les taux de transmission initiaux              
        #B = B
        pass

    facteur = [B, Nu, mu, eta, lambd, alpha, tau, echange]                                    # On met les paramètres dans un vecteur
    
    Nt = tmax                                                                   # On initialise le nombre de points de la grille    
    t = np.linspace(0, tmax, Nt+1)                                              # On initialise la grille de temps

    delta = 2                                                                   # On initialise le pas de temps 
        
    # On initie la premiere population
    Isn0, Ivn0, Isr0, Ivr0 = 1, 1, 1, 1                                         # On initialise le nombre d'infectés    
    I0 = Isn0 + Ivn0 + Isr0 + Ivr0                                         

    Rn0, Rr0 = 0, 0                                                             # On initialise le nombre de recovered               
    R0 = Rn0 + Rr0                                         

    Sn0, Sr0 = (N[0] - I0 - R0)/2, (N[0] - I0 - R0)/2                           # On initialise le nombre de sains
    S0 = Sn0 + Sr0                                      

    Vn0 = Vr0 = 0                                                               # On initialise le nombre de vaccinés          
    V0 = Vn0 + Vr0                         

    Csn, Csr0, Cvn0, Cvr0 = 0, 0, 0, 0                                          # On initialise le nombre de cas confirmés   
    C0 = Csn + Csr0 + Cvn0 + Cvr0           

    Dn0, Dr0 = 0, 0                                                             # On initialise le nombre de décès 
    D0 = Dn0 + Dr0


    # On initie la deuxième population
    Isn_P0, Ivn_P0, Isr_P0, Ivr_P0 = 0, 0, 0, 0                                 # On initialise le nombre d'infectés    
    I_P0 = Isn_P0 + Ivn_P0 + Isr_P0 + Ivr_P0                                         

    Rn_P0, Rr_P0 = 0, 0                                                         # On initialise le nombre de recovered               
    R_P0 = Rn_P0 + Rr_P0                                         

    Sn_P0, Sr_P0 = (N[1] - I_P0 - R_P0)/2, (N[1] - I_P0 - R_P0)/2               # On initialise le nombre de sains
    S_P0 = Sn_P0 + Sr_P0                                      

    Vn0 = Vr0 = 0                                                               # On initialise le nombre de vaccinés          
    V_P0 = Vn0 + Vr0                         

    Csn, Csr0, Cvn0, Cvr0 = 0, 0, 0, 0                                          # On initialise le nombre de cas confirmés   
    C_P0 = Csn + Csr0 + Cvn0 + Cvr0           

    Dn0, Dr0 = 0, 0                                                             # On initialise le nombre de décès 
    D_P0 = Dn0 + Dr0

    X0 = Sn0, Sr0, Vn0, Vr0, Csn, Csr0, Cvn0, Cvr0, Isn0, Ivn0, Isr0, Ivr0, Rn0, Rr0, Dn0, Dr0, Sn_P0, Sr_P0, Vn0, Vr0, Csn, Csr0, Cvn0, Cvr0, Isn_P0, Ivn_P0, Isr_P0, Ivr_P0, Rn_P0, Rr_P0, Dn0, Dr0 # On initialise le vecteur d'état initial
        
    model_rk4 = rk4(SIRCVD_echange, X0, t, N, delta, facteur, geste_barriere, confinement, vaccination)

    return t, model_rk4

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
    
    st.write("Dans le cadre de ce projet nous avons utilisés trois modélès : ")
    st.write("> - Un SIR classique,")
    st.write("> - Un SIR que nous avons modifiés (SIRCVD),")
    st.write("> - Un SIRCVD avec deux populations.")

    _, col2, _ = st.columns([2,8,1])
    model_image = "https://raw.githubusercontent.com/ClovisDel/SIR_Simulation/main/model_global.png?token=GHSAT0AAAAAAB2JO4ZGBERTLFDK3STOP74EY4MQ35Q"
    col2.image(model_image, width = 700)

    st.write("De plus, nous avons ajouté des paramètres supplémentaires pour rendre le modèle plus réaliste : ")
    st.write("> - Des geste barrière, qui consiste à diminuer le taux de transmission de l'épidémie à chaque instant t")
    st.write("> - Un confinement, qui consiste à diminuer fortement le taux de transmission de l'épidémie entre deux instants t non consécutifs")
    st.write("> - Une vaccination, l'idée est de simuler des vagues de vaccinations comme on peut en voir dans la réalité")

    st.write("Pour la vaccination, nous avons transformé le taux de vaccination de tel façon à ce qu'ils suivent une courbe sinusoïdale, pour ce faire on le multiple par un facteur θ tel que : ")
    st.latex(r'''
            \theta(t)=\left| sin(\frac{1}{15}*(t+15)) \right|
            ''')
    
    t = np.linspace(0, 100, 101)
    y = abs( 1/15 * np.sin(1 / 15 * (t-2)))
    fig_courbe_vaccination = px.line(x=t, y=y)
    fig_courbe_vaccination.update_layout(
    title="Evolution du taux θ modifiant la vaccination en fonction de l'instant t :", xaxis_title="t", yaxis_title="taux")   
    st.plotly_chart(fig_courbe_vaccination, use_container_width=True)

def espace_entre_parties():
    st.write("")
    st.markdown("""<hr style="height:2px;width:100%;border:none;color:#333;background-color:#333;" /> """, unsafe_allow_html=True) 
    st.write("")

def load_page_sir_classique(N, tmax, geste_barriere, confinement, vaccination):
    st.header("Simulation du modèle SIR")

    _, col2, _ = st.columns([2,6,1])
    image_SIR = "https://raw.githubusercontent.com/ClovisDel/SIR_Simulation/main/model_SIR.png?token=GHSAT0AAAAAAB2JO4ZGOETFTXLIYQCFGUM6Y4MP42Q"
    col2.image(image_SIR, width = 700)

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

    t, model_rk4_SIR_value = model_rk4_SIR(N, tmax, geste_barriere, confinement, vaccination) 

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

    _, col2, _ = st.columns([2,6,1])
    image_SIRCVD = "https://raw.githubusercontent.com/ClovisDel/SIR_Simulation/main/model_SIRCVD_simple.png?token=GHSAT0AAAAAAB2JO4ZHFPBW7I6LDLFPHO74Y4MP6OA"
    col2.image(image_SIRCVD, width = 700)

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

    st.write("L'ambition de ce modèle est de prendre en compte différents types de personnes.")
    st.write("Dans un premier temps on veut considérer la différence entre les personnes dites 'à risque' et les personnes 'non à risque'. Ces deux populations vont se voir attribuer des taux de contaminations, de guérisons et de mortalité différents.")

    st.write("Dans un second temps on veut prendre en compte les personnes vaccinées, on va donc considérer que les personnes vaccinées peuvent être contaminées mais avec un taux réduit. Leur taux de décès est également réduit et le taux de guérison est augmenté.")
    
    st.write("On découle sur un modèle qui superpose 4 groupes de personnes :")
    st.write("> - Sn les sains non à rique")
    st.write("> - Sr les sains à rique")
    st.write("> - Vr les vaccinés non à rique")
    st.write("> - Vn les vaccinés à rique")

    st.write("D'une part, les sains peuvent se faire vaccinés. D'autres part un infecté quand il est retablis passe à vacciné.")
    st.write("Dans ce modèle, on considère que les personnes décédées restent dans le modèle et ne sont plus comptabilisées dans les autres groupes.")

    espace_entre_parties()

    st.write("A l'aide de la méthode Runge-Kutta de degré 4 (RK4) on observer la simulation suivante :")

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
    df_global_sans_risque = df_global_sans_risque.astype(int)
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
    
    df_global_a_risque = pd.DataFrame ({"t" : t, "S": Sr, "V" : Vr, "C" : Csr + Cvr, "I": Isr + Ivr, "R": Rr, "D": Dr})
    df_global_a_risque = df_global_a_risque.astype(int)
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

    st.button(label="Afficher donnnées", key="btn_scrape")
    col1.plotly_chart(fig_global_sans_risque, use_container_width=True)
    col2.plotly_chart(fig_global_a_risque, use_container_width=True)

    if st.session_state.get("btn_scrape"):
        _, col1, _, col2, _ = st.columns([1,3,3,3,1])
        col1.dataframe(df_global_sans_risque.round(0))
        col2.dataframe(df_global_a_risque.round(0))
    
def load_page_sir_modif_echange(N, tmax, geste_barriere, confinement, vaccination):
    t, model_rka_SIRCVD_echange_value = model_rk4_SIRCVD_echange(N, tmax, geste_barriere, confinement, vaccination) 
    
    Sn, Sr, Vn, Vr, Csn, Csr, Cvn, Cvr, Isn, Ivn, Isr, Ivr, Rn, Rr, Dn, Dr, Sn_P, Sr_P, Vn_P, Vr_P, Csn_P, Csr_P, Cvn_P, Cvr_P, Isn_P, Ivn_P, Isr_P, Ivr_P, Rn_P, Rr_P, Dn_P, Dr_P = model_rka_SIRCVD_echange_value.T

    t = t.astype(int)
    I = (Isn + Ivn + Isr + Ivr).astype(int)
    R = (Rn + Rr).astype(int)
    S = (Sn + Sr).astype(int)
    V = (Vn + Vr).astype(int)
    C = (Csn + Csr + Cvn + Cvr).astype(int)
    D = (Dn + Dr).astype(int)
    df_global = pd.DataFrame ({"t" : t, "S": S, "V" : V, "C" : C, "I": I, "R": R, "D": D})

    I_P = (Isn_P + Ivn_P + Isr_P + Ivr_P).astype(int)
    R_P = (Rn_P + Rr_P).astype(int)
    S_P = (Sn_P + Sr_P).astype(int)
    V_P = (Vn_P + Vr_P).astype(int)
    C_P = (Csn_P + Csr_P + Cvn_P + Cvr_P).astype(int)
    D_P = (Dn_P + Dr_P).astype(int)
    df_global_P = pd.DataFrame ({"t" : t, "S": S_P, "V" : V_P, "C" : C_P, "I": I_P, "R": R_P, "D": D_P})
    
    fig_global = px.line(df_global, x="t", 
                    y=["S", "V", "C", "I", "R", "D"], 
                    title='Simulation du modèle SIRCVD first Pop',
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
    
    fig_global_P = px.line(df_global_P, x="t", 
                    y=["S", "V", "C", "I", "R", "D"], 
                    title='Simulation du modèle SIRCVD Second Pop',
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

    col1, col2 = st.columns(2)

    col1.dataframe(df_global)
    col1.plotly_chart(fig_global, use_container_width=True)
    col2.dataframe(df_global_P)
    col2.plotly_chart(fig_global_P, use_container_width=True)