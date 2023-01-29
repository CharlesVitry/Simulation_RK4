import streamlit as st
import pandas as pd
import numpy as np

from models import *
from utils import *


def espace_entre_parties():
    st.write("")
    st.markdown("""<hr style="height:2px;width:100%;border:none;color:#333;background-color:#333;" /> """, unsafe_allow_html=True) 
    st.write("") 

def sidebar():
    radio_markdown_gestes_barrieres = '''
        Reduit le taux de transmission par **0.4** par instant t. 
        '''.strip()
    radio_markdown_confinement = '''
        Reduit le taux de transmission par **0.05** entre les instants 15 et 30. 
        '''.strip()
    radio_markdown_vaccination = '''
        Modifie le taux de vaccination selon une fonction sinusoïdale.  
        '''.strip()
    st.sidebar.markdown("## Paramètres SIR / SIRCVD")
    N = st.sidebar.slider("Taille de la population :", min_value=50, max_value=1000, value=500, step=50)
    tmax = st.sidebar.slider("Durée de la simulation (en jours) :", min_value=30, max_value=365, value=100, step=5)
        
    st.sidebar.markdown("## Paramètres SIRCVD avec échange")
    N_P = st.sidebar.slider("Taille de la population P :", min_value=50, max_value=1000, value=500, step=50)
    
    st.sidebar.markdown("## Scénarios")
    geste_barriere = st.sidebar.checkbox("Gestes barrières", help=radio_markdown_gestes_barrieres)
    confinement = st.sidebar.checkbox("Confinement", help=radio_markdown_confinement)
    vaccination =  st.sidebar.checkbox("Vaccination", help=radio_markdown_vaccination) 

    return N, N_P, tmax, geste_barriere, confinement, vaccination
    
@st.cache   
def load_data():
    img_SIRCVD_url = "https://drive.google.com/file/d/1myh4YciJwkVKV5PmWdeVY-NuO4wWcK-j/view?usp=sharing"
    img_SIRCVD_path = 'https://drive.google.com/uc?export=download&id='+img_SIRCVD_url.split('/')[-2]
    
    img_SIRCVD_echange_url = "https://drive.google.com/file/d/1AGsBw92clogylm1JyJrzlOJ0SYhZynMK/view?usp=sharing"
    img_SIRCVD_echange_path = 'https://drive.google.com/uc?export=download&id='+img_SIRCVD_echange_url.split('/')[-2]

    return img_SIRCVD_path, img_SIRCVD_echange_path

def IHM():
    st.set_page_config(layout="wide")  

    st.sidebar.title("Menu")

    img_SIRCVD, img_SIRCVD_echange_path = load_data()

    tab1, tab2, tab3, tab4 = st.tabs(["Présentation", "Modèle SIR classique", "Modèle SIRCVD", "Modèle SIRCVD avec échange"])

    N = 500
    N, N_P, tmax, geste_barriere, confinement, vaccination = sidebar()

    with tab1 :
        load_page_accueil()
    with tab2 :         
        load_page_sir(N, tmax, geste_barriere, confinement, vaccination)
    with tab3 :
        load_page_sircvd(N, tmax, geste_barriere, confinement, vaccination, img_SIRCVD)
    with tab4 :
        load_page_sircvd_echange(N, N_P, tmax, geste_barriere, confinement, vaccination, img_SIRCVD_echange_path)

def load_page_accueil():
    st.header("Projet Simulation IMA 5")
    st.markdown("> Clovis Delêtre & Charles Vitry")
    
    st.write("L'ambition de ce projet est de simuler par une approche Systèmes Dynamiques l'évolution d'une épidémie. "
             "Cette méthode se base sur les équations différentielles déterminites établies sur un modèle SIR donnée")
    
    st.write("Dans le cadre de ce projet nous avons utilisés trois modélès : ")
    st.write("> - Un SIR classique,")
    st.write("> - Un SIR que nous avons modifiés (SIRCVD),")
    st.write("> - Un SIRCVD avec deux populations.")

    st.write("De plus, nous avons ajouté des paramètres supplémentaires pour rendre le modèle plus réaliste : ")
    st.write("> - Des geste barrière, qui consiste à diminuer le taux de transmission de l'épidémie à chaque instant t par 0.4")
    st.write("> - Un confinement, qui consiste à diminuer fortement le taux de transmission par 0.05 de l'épidémie entre deux instants t non consécutifs")
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

    st.write("#### Point notable : ")
    st.write("Il est important de garder à l'esprit les limitations de nos modèles. Ces modèles dit *'à compartiments'* sont des modèles simplifiés qui ne prennent pas en compte de nombreux facteurs.")
    st.write("> - Dans notre cas on suppose une population homogène et que les individus se déplacent à taux fixe entre les différents compartiments, ce qui n'est pas le cas dans la réalité.",
            "En effet, les populations sont hétérogènes et la propagation d'une maladie peut être influencée par divers facteurs comme l'âge, le sexe, le niveau de vie, etc.")
    st.write("> - Dans notre cas nous utilisons des paramètres à valeurs fictives, dans un cas concret il faudrait les déterminer à partir de données réelles.")

    st.write("En conclusion, ces modèles sont utiles pour comprendre les mécanismes de propagation d'une épidémie et pour simuler l'évolution d'une épidémie dans le temps, mais ils ne sont pas en mesure de saisir toutes les complexités des épidémies du monde réel.")
def load_page_sir(N, tmax, geste_barriere, confinement, vaccination):
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
    col1, col2 = st.columns([2, 2])
    col1.latex(r'''
    \begin{align*}
\frac{dS}{dt} &= -\beta S I \\
\frac{dI}{dt} &= \beta S I - \lambda I \\
\frac{dR}{dt} &= \lambda I
\end{align*}
''')
    col2.write("Avec : ")
    col2.write("> - β : taux d'infection")
    col2.write("> - λ : taux de guérison")

    espace_entre_parties()

    st.write("A l'aide de la méthode Runge-Kutta de degré 4 (RK4) on observer la simulation suivante :")
    col1 , col2, col3 = st.columns([2, 2, 15])
    beta = col1.number_input("β = ", 0.4)
    lamb = col2.number_input("λ = ", 10)

    I0 = col1.number_input("I0 = ", 1)
    R0 = col2.number_input("R0 = ", 0)

    if geste_barriere :
        beta = beta *(1 - 0.4)   
    facteur = [beta, lamb] 
    X0 = [N-I0-R0, I0, R0]
    t = np.linspace(0, tmax, tmax+1)

    model = model_rk4_sir(N, facteur, X0, t, 2, geste_barriere, confinement, vaccination)
    df, fig = transform_model_sir(t, model)
    col3.plotly_chart(fig, use_container_width=True)


def load_page_sircvd(N, tmax, geste_barriere, confinement, vaccination, img_SIRCVD):
    st.header("Simulation du modèle SIR adapté : SIRCVD")

    try :
        _, col, _ = st.columns([1, 3, 1])
        col.image(img_SIRCVD, caption='img_SIRCVD', use_column_width=True)
    except : 
        pass
        


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
    col2.write("> - η le taux de vaccinés qui ne sont plus immunisés")
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

    st.write("#### A l'aide de la méthode Runge-Kutta de degré 4 (RK4) on observe la simulation suivante :")

    col1, col2, col3, col4, col5, col6, col7, col8, col9= st.columns([3, 3, 3, 3, 1, 3, 3, 3, 3])
    betaSn = col1.number_input("β(Sn) = ", value = 0.7)
    betaSr = col2.number_input("β(Sr) = ", value =  0.8)
    betaVn = col3.number_input("β(Vn) = ", value = 0.4)
    betaVr = col4.number_input("β(Vr) = ", value = 0.5)
    muSn = col1.number_input("µ(Sn) = ", value = 0.002, step=0.001, format='%.3f')
    muSr = col2.number_input("µ(Sr) = ", value = 0.003, step=0.001, format='%.3f')
    muVn = col3.number_input("µ(Vn) = ", value = 0.001, step=0.001, format='%.3f')
    muVr = col4.number_input("µ(Vr) = ", value = 0.001, step=0.001, format='%.3f')
    lambdaSn = col1.number_input("λ(Sn) = ", value = 15)
    lambdaSr = col2.number_input("λ(Sr) = ", value = 20)
    lambdaVn = col3.number_input("λ(Vn) = ", value = 10)
    lambdaVr = col4.number_input("λ(Vr) = ", value = 15)    
    nu = col1.number_input("ν = ", value = 10)
    eta = col2.number_input("η = ", value = 0.1)
    alpha = col3.number_input("α = ", value = 200)
    tau = col4.number_input("τ = ", value = 50)

    I_Sn_0 = col6.number_input("I(Sn)0 = ", value = 1)
    I_Sr_0 = col7.number_input("I(Sr)0 = ", value = 1)
    I_Vn_0 = col8.number_input("I(Vn)0 = ", value = 1)
    I_Vr_0 = col9.number_input("I(Vr)0 = ", value = 1)
    R_n_0 = col6.number_input("R(n)0 = ", value = 0)
    R_r_0 = col7.number_input("R(r)0 = ", value = 0)
    D_n_0 = col8.number_input("D(n)0 = ", value = 0)
    D_r_0 = col9.number_input("D(r)0 = ", value = 0)
    C_Sn_0 = col6.number_input("C(Sn)0 = ", value = 0)
    C_Sr_0 = col7.number_input("C(Sr)0 = ", value = 0)
    C_Vn_0 = col8.number_input("C(Vn)0 = ", value = 0)
    C_Vr_0 = col9.number_input("C(Vr)0 = ", value = 0)
    V_n_0 = col6.number_input("V(n)0 = ", value = 0)
    V_r_0 = col7.number_input("V(r)0 = ", value = 0)
    
    S_n_0, S_r_0 = (N - I_Sn_0 - I_Sr_0 - I_Vn_0 - I_Vr_0 - R_n_0 - R_n_0)/2, (N - I_Sn_0 - I_Sr_0 - I_Vn_0 - I_Vr_0 - R_n_0 - R_n_0)/2
    
    X0 = [S_n_0, S_r_0, V_n_0, V_r_0, C_Sn_0, C_Sr_0, C_Vn_0, C_Vr_0, I_Sn_0, I_Sr_0, I_Vn_0, I_Vr_0, R_n_0, R_r_0, D_n_0, D_r_0]
    
    B = [betaSn, betaSr, betaVn, betaVr]
    if geste_barriere:                                                          
        B = [B[i]*(1 - 0.4) for i in range(len(B))]

    facteur = [B,  [muSn, muSr, muVn, muVr], [lambdaSn, lambdaSr, lambdaVn, lambdaVr], nu, eta, alpha, tau]
    t = np.linspace(0, tmax, tmax+1)
    model = model_rk4_sircvd(N, facteur, X0, t, 2, geste_barriere, confinement, vaccination)
    df_global, df_global_merge, fig_global_merge, df_global_sans_risque, fig_global_sans_risque, df_global_avec_risque, fig_global_avec_risque = transform_model_sircvd_global(t, model)
    st.plotly_chart(fig_global_merge, use_container_width=True)

    #st.write(facteur)
    #st.write(X0)

    col1, col2 = st.columns([1, 1])
    col1.plotly_chart(fig_global_avec_risque, use_container_width=True)
    col2.plotly_chart(fig_global_sans_risque, use_container_width=True)

    st.button(label = "Voir les données", key = "button_df")
    if st.session_state.get("button_df"):
        st.dataframe(df_global)
    #st.dataframe(df_global)
    #st.write(facteur)

    st.write("#### Analyse du modèle :")
    st.write("Globalement, le modèle défini permet bien de représenter plusieurs populations lors d'une épidémie.",
            "On remarque une forte vague de contamination en début de période qui arrive à se stabiliser au cours du temps.")
    st.write("En comparant les deux populations, on remarque que la population à risque est plus touchée que la population sans risque.",
            "Cela est dû au fait que la population à risque est plus vulnérable au virus et donc plus touchée par celui-ci.",
            "De ce fait, on remarque des vagues plus importantes et plus longue dans le temps. On remarque également une mortalité supérieur.")

    st.write("Le scénario avec gestes barrières montre bien une diminuation de la vague de contamination (~ d'infectation) sur le temps et en intensité.")

    st.write("Le scénario du confinement est quant à lui très efficace, en effet de la période 15 à 30 de confinement la contamination est quasi nulle.",
            "Cependant, un confinement seul n'a pour effet que de décaler la vague de contamination et de limiter légèrement la suivante.")

    st.write("Le scénario avec la vaccination sinusoïdale montre bien de forte montée de la population de vaccinées au cours du temps.", 
            "On remarque que la première vague de vaccination est plus importante que les autres." ,
            "Ce phénomène peut s'expliquer par le fait que la population de *sains* est plus importantes en début de période.",
            "De plus, on remarque de légère diminution de la contamination (et donc d'infection) au cours du temps."
            ) 
    
    st.write("Globalement, on remarque que les gestes barrières et le confinements sont les deux scénarios les plus efficaces pour limiter la propagation du virus.", 
            "La combinaison de plusieurs scénarios restent privilégiable pour limiter la propagation du virus.")



def load_page_sircvd_echange(N, N_P, tmax, geste_barriere, confinement, vaccination, img_SIRCVD_echange_path):
    st.header("Simulation du modèle SIR adapté : SIRCVD avec échange")

    try :
        _, col, _ = st.columns([1, 3, 1])
        col.image(img_SIRCVD_echange_path, caption='img_SIRCVD', use_column_width=True)
    except : 
        pass

    espace_entre_parties()

    st.write("#### Présentation du modèle :")
    st.write("L'ambition de ce modèle est de comparer deux pays / groupe d'individus qui correspondent tous les deux à des modèles SIRCVD vus précédemment.")
    st.write("Dans notre cas, on va estimer que deux pays peuvent s'échanger des individus contaminés entre eux selon une probabilité φ mais pas des autres 'cases'. ",
            "")
    

    espace_entre_parties()
    st.write("#### A l'aide de la méthode Runge-Kutta de degré 4 (RK4) on observe la simulation suivante :")
    st.write("*Pour faciliter l'utilisation, les données et paramètres sont pré-remplis mais peuvent être modifiés en cliquant sur le bouton **Choix des paramètres**.*")
    

    #st.button(label = "Choix des paramètres", key = "button_choix_param")

    #if st.button(label = "Choix des paramètres", key = "button_choix_param"):
    if True : 
        col1, col2, col3, col4, _, col6, col7, col8, col9= st.columns([3, 3, 3, 3, 1, 3, 3, 3, 3])
        
        betaSn = col1.number_input("β(Sn) = ", value = 0.7, key="unique_key1")
        betaSr = col2.number_input("β(Sr) = ", value =  0.8, key="unique_key2")
        betaVn = col3.number_input("β(Vn) = ", value = 0.4, key="unique_key3")
        betaVr = col4.number_input("β(Vr) = ", value = 0.5, key="unique_key4")
        muSn = col1.number_input("µ(Sn) = ", value = 0.002, step=0.001, format='%.3f', key="unique_key5")
        muSr = col2.number_input("µ(Sr) = ", value = 0.003, step=0.001, format='%.3f', key="unique_key6")
        muVn = col3.number_input("µ(Vn) = ", value = 0.001, step=0.001, format='%.3f', key="unique_key7")
        muVr = col4.number_input("µ(Vr) = ", value = 0.001, step=0.001, format='%.3f', key="unique_key8")
        lambdaSn = col1.number_input("λ(Sn) = ", value = 15, key="unique_key9")
        lambdaSr = col2.number_input("λ(Sr) = ", value = 20, key="unique_key10")
        lambdaVn = col3.number_input("λ(Vn) = ", value = 10, key="unique_key11")
        lambdaVr = col4.number_input("λ(Vr) = ", value = 15, key="unique_key12")    
        nu = col1.number_input("ν = ", value = 10, key="unique_key13")
        eta = col2.number_input("η = ", value = 0.1, key="unique_key14")
        alpha = col3.number_input("α = ", value = 200, key="unique_key15")
        tau = col4.number_input("τ = ", value = 50, key="unique_key16")
        phi = col1.number_input("φ = ", value = 0.001, key="unique_key61", format='%.3f')

        betaSn_P = col6.number_input("β_P(Sn) = ", value = 0.7, key="unique_key17")
        betaSr_P = col7.number_input("β_P(Sr) = ", value =  0.8, key="unique_key18")
        betaVn_P = col8.number_input("β_P(Vn) = ", value = 0.4, key="unique_key19")
        betaVr_P = col9.number_input("β_P(Vr) = ", value = 0.5, key="unique_key20")
        muSn_P = col6.number_input("µ_P(Sn) = ", value = 0.002, step=0.001, format='%.3f', key="unique_key21")
        muSr_P = col7.number_input("µ_P(Sr) = ", value = 0.003, step=0.001, format='%.3f', key="unique_key22")
        muVn_P = col8.number_input("µ_P(Vn) = ", value = 0.001, step=0.001, format='%.3f', key="unique_key23")
        muVr_P = col9.number_input("µ_P(Vr) = ", value = 0.001, step=0.001, format='%.3f', key="unique_key24")
        lambdaSn_P = col6.number_input("λ_P(Sn) = ", value = 15, key="unique_key25")
        lambdaSr_P = col7.number_input("λ_P(Sr) = ", value = 20, key="unique_key26")
        lambdaVn_P = col8.number_input("λ_P(Vn) = ", value = 10, key="unique_key27")
        lambdaVr_P = col9.number_input("λ_P(Vr) = ", value = 15, key="unique_key28")
        nu_P = col6.number_input("ν_P = ", value = 10, key="unique_key29")
        eta_P = col7.number_input("η_P = ", value = 0.1, key="unique_key30")
        alpha_P = col8.number_input("α_P = ", value = 200, key="unique_key31")
        tau_P = col9.number_input("τ_P = ", value = 50, key="unique_key32")
        phi_P = phi

        col1, col2, col3, col4, _, col6, col7, col8, col9= st.columns([3, 3, 3, 3, 1, 3, 3, 3, 3])

        I_Sn_0 = col1.number_input("I(Sn)0 = ", value = 0, key="unique_key33")
        I_Sr_0 = col2.number_input("I(Sr)0 = ", value = 0, key="unique_key34")
        I_Vn_0 = col3.number_input("I(Vn)0 = ", value = 0, key="unique_key35")
        I_Vr_0 = col4.number_input("I(Vr)0 = ", value = 0, key="unique_key36")
        R_n_0 = col1.number_input("R(n)0 = ", value = 0, key="unique_key37")
        R_r_0 = col2.number_input("R(r)0 = ", value = 0, key="unique_key38")
        D_n_0 = col3.number_input("D(n)0 = ", value = 0, key="unique_key39")
        D_r_0 = col4.number_input("D(r)0 = ", value = 0, key="unique_key40")
        C_Sn_0 = col1.number_input("C(Sn)0 = ", value = 0, key="unique_key41")
        C_Sr_0 = col2.number_input("C(Sr)0 = ", value = 0, key="unique_key42")
        C_Vn_0 = col3.number_input("C(Vn)0 = ", value = 0, key="unique_key43")
        C_Vr_0 = col4.number_input("C(Vr)0 = ", value = 0, key="unique_key44")
        V_n_0 = col1.number_input("V(n)0 = ", value = 0, key="unique_key45")
        V_r_0 = col2.number_input("V(r)0 = ", value = 0, key="unique_key46")

        I_Sn_0_P = col6.number_input("I_P(Sn)0 = ", value = 1, key="unique_key47")
        I_Sr_0_P = col7.number_input("I_P(Sr)0 = ", value = 1, key="unique_key48")
        I_Vn_0_P = col8.number_input("I_P(Vn)0 = ", value = 1, key="unique_key49")
        I_Vr_0_P = col9.number_input("I_P(Vr)0 = ", value = 1, key="unique_key50")
        R_n_0_P = col6.number_input("R_P(n)0 = ", value = 0, key="unique_key51")
        R_r_0_P = col7.number_input("R_P(r)0 = ", value = 0, key="unique_key52")
        D_n_0_P = col8.number_input("D_P(n)0 = ", value = 0, key="unique_key53")
        D_r_0_P = col9.number_input("D_P(r)0 = ", value = 0, key="unique_key54")
        C_Sn_0_P = col6.number_input("C_P(Sn)0 = ", value = 0, key="unique_key55")
        C_Sr_0_P = col7.number_input("C_P(Sr)0 = ", value = 0, key="unique_key56")
        C_Vn_0_P = col8.number_input("C_P(Vn)0 = ", value = 0, key="unique_key57")
        C_Vr_0_P = col9.number_input("C_P(Vr)0 = ", value = 0, key="unique_key58")
        V_n_0_P = col6.number_input("V_P(n)0 = ", value = 0, key="unique_key59")
        V_r_0_P = col7.number_input("V_P(r)0 = ", value = 0, key="unique_key60")

        S_n_0, S_r_0 = (N - I_Sn_0 - I_Sr_0 - I_Vn_0 - I_Vr_0 - R_n_0 - R_n_0)/2, (N - I_Sn_0 - I_Sr_0 - I_Vn_0 - I_Vr_0 - R_n_0 - R_n_0)/2
        S_n_P_0, S_r_P_0 = (N - I_Sn_0_P - I_Sr_0_P - I_Vn_0_P - I_Vr_0_P - R_n_0_P - R_n_0_P)/2, (N - I_Sn_0_P - I_Sr_0_P - I_Vn_0_P - I_Vr_0_P - R_n_0_P - R_n_0_P)/2
        
        X0 = [S_n_0, S_r_0, I_Sn_0, I_Sr_0, I_Vn_0, I_Vr_0, R_n_0, R_r_0, D_n_0, D_r_0, C_Sn_0, C_Sr_0, C_Vn_0, C_Vr_0, V_n_0, V_r_0, S_n_P_0, S_r_P_0, I_Sn_0_P, I_Sr_0_P, I_Vn_0_P, I_Vr_0_P, R_n_0_P, R_r_0_P, D_n_0_P, D_r_0_P, C_Sn_0_P, C_Sr_0_P, C_Vn_0_P, C_Vr_0_P, V_n_0_P, V_r_0_P]
        #st.write(X0)
        B = [betaSn, betaSr, betaVn, betaVr]
        B_P = [betaSn_P, betaSr_P, betaVn_P, betaVr_P]
        if geste_barriere:                                                          
            B = [B[i]*(1 - 0.4) for i in range(len(B))]
            B_P = [B_P[i]*(1 - 0.4) for i in range(len(B_P))]
        
        facteur = [B,  [muSn, muSr, muVn, muVr], [lambdaSn, lambdaSr, lambdaVn, lambdaVr], nu, eta, alpha, tau, phi, B_P, [muSn_P, muSr_P, muVn_P, muVr_P], [lambdaSn_P, lambdaSr_P, lambdaVn_P, lambdaVr_P], nu_P, eta_P, alpha_P, tau_P, phi_P]
        print(facteur)
        #st.write(facteur)
    else : 
        facteur = [[0.7, 0.8, 0.4, 0.5], [0.002, 0.003, 0.001, 0.001], [15, 20, 10, 15], 10, 0.1, 200, 50, 0.001, [0.7, 0.8, 0.4, 0.5], [0.002, 0.003, 0.001, 0.001], [15, 20, 10, 15], 10, 0.1, 200, 50, 0.1]
        X0 = [N/2, N/2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (N-4)/2, (N-4)/2, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    
    t = np.linspace(0, tmax, tmax+1)
    N = [N, N_P]
    model = model_rk4_sircvd_echange(N, facteur, X0, t, 2, geste_barriere, confinement, vaccination)

    df_global, df_merge, fig_merge, df_classique, fig_classique, df_P, fig_P = transform_model_sircvd_echange(t, model)
    st.plotly_chart(fig_merge, use_container_width=True)

    col1, col2 = st.columns(2)
    col1.plotly_chart(fig_classique, use_container_width=True)
    col2.plotly_chart(fig_P, use_container_width=True)

    espace_entre_parties()

    st.write("#### Conclusion :")

    st.write("Dans notre modèle SIRCDV, notre population nationale ne comportent aucun infecté à t=0, on remarque que c’est lors du pic d’infection dans la population internationale que l’épidémie débute au sein de la population nationale. De plus, les échanges avec la population internationale limitent la rechute du nombre d’infectés au sein de la population nationale.")
    st.write("Les scénarios s’appliquant de la même façon que précédemment on peut en déduire les mêmes retombés sur l’épidémie. ")
    st.write("Ce qui est intéressant est de faire varier les paramètres de l’épidémie par pays. ")