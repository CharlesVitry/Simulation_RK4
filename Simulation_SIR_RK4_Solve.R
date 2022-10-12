library(tidyverse)
source("SIR_implementation.R")
source("IHM.R")

# Variables Globales
tailleDeLaPopulationTotale <- 65 * 10^6 #Taille de la Population 
sous_population <- list(
  sains = tailleDeLaPopulationTotale - 15000, #Taille des sous-populations
  infectes = 15000,
  retablis = 0)  
LongueurPeriodeEtude <- 200 #Quelle est la période d'étude des sous population ? : 6 mois
pasDeTemps <- 0.1 #actualisation des équations différentielles
periodeEtude = rep(seq(from = 0, to = LongueurPeriodeEtude, by = 0.1), 4) # 1 fois par scénario.
facteurAleatoire = 280

# Coefficients dans les équations différentiels
TauxRetablissement <- 0.03
TauxInfection <- 0.000000004
TauxEfficaciteVaccin <- 0.007
#paste("avec confinement le ",DebutConfinement1, "ème jour" )
# Pondération des effets de scénarios
reactiviteDuGouvernementVaccin <- 1 #début de la campagne de vaccination
reactiviteDuGouvernementGestesBarrieres <- 0

ReductionInfectionParGestesBarrieresEnPourcent <- 0.2
ReductionInfectionParConfinementEnPourcent <- 0.95

DebutConfinement1 <- 30 # début le 30ème jour du confinement
DureeConfinement1 <- 40 # 40 jours de confinement

scenarios = expand.grid(0:1, 0:1, 0:1) #Booléens de scénarios
scenarios = scenarios[rep(seq_len(nrow(scenarios)),each = LongueurPeriodeEtude*10+1), ]

# Création de la table
SIR_donnees <- data.frame(periodeEtude,scenarios, TauxRetablissement, TauxInfection,TauxEfficaciteVaccin, sous_population)

#Ajout des données de scénario
SIR_donnees <- SIR_donnees %>%
  mutate(
    Var1 = as.logical(Var1),  #de 0,1 à TRUE/FALSE : passage en booléen   
    Var2 = as.logical(Var2),
    Var3 = as.logical(Var3)
  ) %>% 
  rename(
    Vaccin = Var1,     #On renomme les booléens 
    Gestes_Barrieres = Var2,
    Confinement = Var3 
  ) %>% 
  mutate(
    #Effet Vaccin
    TauxEfficaciteVaccin = ifelse( Vaccin & periodeEtude > reactiviteDuGouvernementVaccin , -(TauxEfficaciteVaccin/periodeEtude) +TauxEfficaciteVaccin , 0 ),

    #Effet Gestes barrieres
    TauxInfection = ifelse(Gestes_Barrieres & periodeEtude > reactiviteDuGouvernementGestesBarrieres , TauxInfection * (1 - ReductionInfectionParGestesBarrieresEnPourcent)  ,TauxInfection),
    
    #Effet Confinement
    TauxInfection = ifelse(Confinement & periodeEtude > DebutConfinement1 & periodeEtude < DebutConfinement1 + DureeConfinement1 , TauxInfection * (1 - ReductionInfectionParConfinementEnPourcent)  ,TauxInfection),
    
    #Variable de Scénario issus de la combinaisons des booléens
    scenario = paste (as.character(factor(Vaccin, labels = c("Pas de vaccin", "Avec Vaccin"))),",",
                      as.character( factor(Gestes_Barrieres, labels = c("sans gestes barrièree", "avec gestes barrière"))),"et",
                      as.character( factor(Confinement, labels = c("sans confinement", "avec confinement")) ))
    # ,
    # #Ajout d'aléatoire
    # TauxInfection = abs(jitter(TauxInfection,factor = facteurAleatoire)) ,
    # TauxRetablissement =abs( jitter(TauxRetablissement,factor = facteurAleatoire))
  )

#Application du modèle SIR
for (i in 1:nrow(SIR_donnees)) { #transformation réalisable avec mutate+lag() mais moins compréhensible pour mon binome
  if (SIR_donnees$periodeEtude[i] == 0) next
  
  Actualisation <- EDO_RK4_Vaccin(
    SIR_donnees$sains[i - 1],
    SIR_donnees$infectes[i - 1],
    SIR_donnees$retablis[i - 1],
    SIR_donnees$TauxInfection[i-1],
    SIR_donnees$TauxRetablissement[i-1],
    SIR_donnees$TauxEfficaciteVaccin[i-1])
  
  SIR_donnees$sains[i] <- Actualisation$NbreSains_KPlusUn
  SIR_donnees$infectes[i] <- Actualisation$NbreInfecte_KPlusUn
  SIR_donnees$retablis[i] <- Actualisation$NbreRetablis_KPlusUn
}
#Transformation données
SIR_donnees <- SIR_donnees %>% 
  pivot_longer(
    cols = c(sains,infectes,retablis),
    names_to = "sous_population",
    values_to = "Nombre")

#Affichage
plot_SIR(SIR_donnees)
ggsave(plot_SIR(SIR_donnees),filename =  "plot.pdf",width = 14, height = 8, device = cairo_pdf)