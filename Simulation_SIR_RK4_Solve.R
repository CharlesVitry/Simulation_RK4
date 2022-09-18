library(tidyverse)
source("SIR_implementation.R")
source("IHM.R")

#https://blog.tonytsai.name/blog/2014-11-24-rk4-method-for-solving-sir-model/
# Variables Globales
tailleDeLaPopulationTotale <- 65 * 10^6
sous_population <- list(
  sains = tailleDeLaPopulationTotale - 1,
  infectes = 10,
  retablis = 0)
LongueurPeriodeEtude <- 182 # 6 mois
pasDeTemps <- 0.1
periodeEtude = rep(seq(from = 0, to = LongueurPeriodeEtude, by = 0.1), 4) # 1 fois par scénario.
TauxRetablissement <- 0.3
TauxInfection <- 0.6

effetVaccinInfection <- 1.5
effetVaccinRetablissement <- .7
effetPolitiqueSanitaireInfection <- .7 #30% réduction
reactiviteDuGouvernement <- 30 #jours

scenario = c(rep("Scénario_SansVaccinSansPolitiqueSanitaire", LongueurPeriodeEtude*10+1),
             rep("Scénario_AvecVaccinSansPolitiqueSanitaire", LongueurPeriodeEtude*10+1),
             rep("Scénario_SansVaccinAvecPolitiqueSanitaire", LongueurPeriodeEtude*10+1),
             rep("Scénario_AvecVaccinAvecPolitiqueSanitaire", LongueurPeriodeEtude*10+1))

# Création de la table
SIR_donnees <- data.frame(periodeEtude,scenario, TauxRetablissement, TauxInfection, sous_population)

#Ajout des données de scénario
SIR_donnees <- SIR_donnees %>%
  mutate(
    Vaccin = (scenario == "Scénario_AvecVaccinSansPolitiqueSanitaire")| (scenario == "Scénario_AvecVaccinAvecPolitiqueSanitaire") ,
    PolitiqueSanitaire = (scenario == "Scénario_SansVaccinAvecPolitiqueSanitaire")| (scenario == "Scénario_AvecVaccinAvecPolitiqueSanitaire"),
    
    #Effet Politique Sanitaire
    TauxInfection = ifelse(PolitiqueSanitaire & periodeEtude > reactiviteDuGouvernement, TauxInfection * effetPolitiqueSanitaireInfection ,TauxInfection),
    
    #Effet Vaccin
    TauxInfection = ifelse(Vaccin & periodeEtude > reactiviteDuGouvernement, TauxInfection * (1/log10(periodeEtude)*effetVaccinInfection) ,TauxInfection),
    
  )


#Application du modèle SIR
for (i in 1:nrow(SIR_donnees)) { #transformation réalisable avec mutate+lag() mais moins compréhensible pour mon binome
  if (SIR_donnees$periodeEtude[i] == 0) next
  
  Actualisation <- EDO_RK4(
    SIR_donnees$sains[i - 1],
    SIR_donnees$infectes[i - 1],
    SIR_donnees$retablis[i - 1],
    SIR_donnees$TauxInfection[i-1],
    SIR_donnees$TauxRetablissement[i-1])
  SIR_donnees$sains[i] <- Actualisation$NbreSains_PlusUn
  SIR_donnees$infectes[i] <- Actualisation$NbreInfecte_PlusUn
  SIR_donnees$retablis[i] <- Actualisation$NbreRetablis_PlusUn
}
#Transformation données
SIR_donnees <- SIR_donnees %>% 
  pivot_longer(
    cols = c(sains,infectes,retablis),
    names_to = "sous_population",
    values_to = "Nombre"
  )

#Affichage
plot_SIR(SIR_donnees)
ggsave(plot_SIR(SIR_donnees),filename =  "plot.pdf",width = 10, height = 8, device = cairo_pdf)