library(tidyverse) #Opérations sur les données
library("rcartocolor") #palette de couleurs spécifique
source("SIR_implementation.R")
source("IHM.R")

# Variables Globales SIR
tailleDeLaPopulationTotale <- 65 * 10^6 #Taille de la Population 
tailleDeLaPopulationEtrangere <- 130 * 10^6
LongueurPeriodeEtude <- 200 #Quelle est la période d'étude des sous population ? : 6 mois
pasDeTemps <- 0.1 #actualisation des équations différentielles
periodeEtude = rep(seq(from = 0, to = LongueurPeriodeEtude, by = 0.1), 4) # 1 fois par scénario.
facteurAleatoire = 280

# Coefficients dans les équations différentiels SIR
TauxRetablissement <- 0.03
TauxInfectionSIR <- 0.000000004
TauxEfficaciteVaccin <- 0.007

# Coefficients dans les équations différentiels SIRCDV
TauxContamination    = 0.00000001
TauxInfection        = 15
TauxRetablissementSN = 5
TauxRetablissementSR = 70
TauxVaccination      = 3000
TauxMortaliteSN      = 0.000001
TauxMortaliteVN      = 0.00005
TauxMortaliteSR      = 0.0005
TauxMortaliteVR      = 0.0001
TauxImmunite         = 5
TauxEchangeDePopulation = 0.003

# Pondération des effets de scénarios
reactiviteDuGouvernementVaccin <- 1 #début de la campagne de vaccination
reactiviteDuGouvernementGestesBarrieres <- 0

ReductionInfectionParGestesBarrieresEnPourcent <- 0.2
ReductionInfectionParConfinementEnPourcent <- 0.95

DebutConfinement1 <- 30 # début le 30ème jour du confinement
DureeConfinement1 <- 40 # 40 jours de confinement

scenarios = expand.grid(0:1, 0:1, 0:1) #Booléens de scénarios
scenarios = scenarios[rep(seq_len(nrow(scenarios)),each = LongueurPeriodeEtude*10+1), ]

sous_population_SIR <- list(
  sains = tailleDeLaPopulationTotale - 15000, #Taille des sous-populations
  infectes = 15000,
  retablis = 0) 
sous_population_SIRCDV <- list(SN  = tailleDeLaPopulationTotale/2,SR  = tailleDeLaPopulationTotale/2,ISN = 0,IVN = 0,IVR = 0,ISR = 0,RN  = 0,RR  = 0,CSN = 0,CVN = 0,CVR = 0,CSR = 0,VN  = 0,VR  = 0,DN  = 0,DR  = 0,SN_P  = tailleDeLaPopulationEtrangere/2 -150000,SR_P  = tailleDeLaPopulationEtrangere/2 -150000,ISN_P = 150000,IVN_P = 0,IVR_P = 0,ISR_P = 150000,RN_P  = 0,RR_P  = 0,CSN_P = 0,CVN_P = 0,CVR_P = 0,CSR_P = 0,VN_P  = 0,VR_P  = 0,DN_P  = 0,DR_P  = 0)

# Création de la table
SIR_donnees <- data.frame(periodeEtude,scenarios, TauxRetablissement, TauxInfectionSIR,TauxEfficaciteVaccin, sous_population_SIR)
SIRCDV_donnees <- data.frame(periodeEtude,scenarios, TauxContamination, TauxInfection, TauxRetablissementSN, TauxRetablissementSR, TauxVaccination, TauxMortaliteSN, TauxMortaliteVN, TauxMortaliteSR,  TauxMortaliteVR, TauxImmunite,TauxEchangeDePopulation, sous_population_SIRCDV)

#Application des scénarios
SIR_donnees = scenario(SIR_donnees)
SIRCDV_donnees = scenario(SIRCDV_donnees)

#Application du modèle SIR
for (i in 1:nrow(SIR_donnees)) { #transformation réalisable avec mutate+lag() mais moins comprehensible.
  if (SIR_donnees$periodeEtude[i] == 0) next
  
  SIR_donnees[i,c(8:10)] <- EDO_RK4_Vaccin(
    SIR_donnees$sains[i - 1],
    SIR_donnees$infectes[i - 1],
    SIR_donnees$retablis[i - 1],
    SIR_donnees$TauxInfectionSIR[i-1],
    SIR_donnees$TauxRetablissement[i-1],
    SIR_donnees$TauxEfficaciteVaccin[i-1]
    )
}

#Application du modèle SIRCDV
for (i in 1:nrow(SIRCDV_donnees)) {
  if (SIR_donnees$periodeEtude[i] == 0) next
  SIRCDV_donnees[i,c(16:47)] <- SIRCDV_RK4(SIRCDV_donnees[i-1,c(5:47)])}

#Transformation données SIR 
SIR_donnees <- SIR_donnees %>% 
  pivot_longer(
    cols = c(sains,infectes,retablis),
    names_to = "sous_population",
    values_to = "Nombre")

#Transformation données SIRCDV
SIRCDV_donnees <- SIRCDV_donnees %>% 
  mutate(
    Contaminés = CSN + CSR + CVN + CVR,
    Décédés = DN + DR,
    infectes = ISN + ISR + IVN + IVR,
    Rétablis = RN + RR,
    Sains = SN + SR,
    Vaccinnés = VN + VR,
    `Contaminés à l'international` = CSN_P + CSR_P + CVN_P + CVR_P,
    `Infectés à l'international` = ISN_P + ISR_P + IVN_P + IVR_P
  ) %>% 
  pivot_longer(
    cols = "Contaminés":"Infectés à l'international",
    names_to = "sous_population",
    values_to = "Nombre") %>% 
  mutate(
    SecondePopulation = grepl("international" , sous_population )
  )

#Affichage SIR
#plot_SIR(SIR_donnees,FALSE)
#ggsave(plot_SIR(SIR_donnees),filename =  "plot.pdf",width = 14, height = 8, device = cairo_pdf)

#Affichage SIRCDV
plot_SIR(SIRCDV_donnees,TRUE)
ggsave(plot_SIR(SIRCDV_donnees,TRUE),filename =  "plot2.pdf",width = 15, height = 16, device = cairo_pdf)
