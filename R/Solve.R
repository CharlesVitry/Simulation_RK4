# Variables Globales
N = 65 * 10^6   #Taille de la Population 
N_E = 5 * 10^6 #Taille de la Population Etrangère
facteurAleatoire = 280 #bruit ajouté dans la variation des taux
tempsEtude = 200 #Quelle est la période d'étude des sous population ? : 6 mois
dt = 0.5 #Pas de temps pour l'actualisation des équations différentielles
t = seq(from = 0, to = tempsEtude, by = dt) ; nb_iter = tempsEtude/dt +1 

#Sous populations
pop_SIR = matrix(data=NA,nrow=nb_iter ,ncol=3)
colnames(pop_SIR) = c("S", "I", "R")
pop_SIR[1,] = c(N - 15000,15000,0)

pop_SIRCDV = matrix(data=NA,nrow=nb_iter ,ncol=32)
colnames(pop_SIRCDV) = c("SN", "SR", "ISN","IVN","IVR","ISR","RN","RR","CSN","CVN","CVR","CSR","VN","VR","DN","DR","SN_P","SR_P","ISN_P","ISR_P","IVN_P","IVR_P","RN_P","RR_P","CSN_P","CVN_P","CVR_P","CSR_P","VN_P","VR_P","DN_P","DR_P")
pop_SIRCDV[1,] =c(N/2,N/2,rep(0,14),N_E/2 -150000,N_E/2 -150000,150000,150000,rep(0,12))

# Méthodes de résolutions RK2 et RK4
rk2 <- function(f, pop, taux){
  k1 = dt * f(pop,taux)
  k2 = dt * f(pop + k1,taux)
  return(pop + (k1 + k2)/2)}

rk4 <- function(f, pop, taux){
  k1 = dt * f(pop,taux)
  k2 = dt * f(pop + k1/2,taux)
  k3 = dt * f(pop + k2/2,taux)
  k4 = dt * f(pop + k3,taux)
  return(pop + (k1 + 2*k2 + 2*k3 + k4)/6)}

#Implémentation des équations différentiels
SIR <- function(p,t){
  d_f = c(
    -t["v"] * p["I"] * p["S"] / N - t["α"] * p["S"] ,
     t["v"] * p["I"] * p["S"] / N - p["I"]/ t["λ"],
     p["I"]/t["λ"]   + t["α"] * p["S"])}

SIRCDV <- function(p, t){
  I   = p["ISN"]   + p["IVN"]   + p["IVR"]   + p["ISR"]
  I_P = p["ISN_P"] + p["IVN_P"] + p["IVR_P"] + p["ISR_P"]
  
  d_f = c(
  #France métropolitaine
  -t["BSN"]* I * p["SN"] / N - p["SN"] / t["α"] + t["η"] * p["VN"]  ,
  -t["BSR"]* I * p["SR"] / N - p["SR"] / t["α"] + t["η"] * p["VR"],
  
  p["CSN"] / t["v"] - p["ISN"] / t["λSN"] - p["ISN"] * t["µSN"],
  p["CVN"] / t["v"] - p["IVN"] / t["λVN"] - p["IVN"] * t["µVN"],
  p["CVR"] / t["v"] - p["IVR"] / t["λVR"] - p["IVR"] * t["µVR"],
  p["CSR"] / t["v"] - p["ISR"] / t["λSR"] - p["ISR"] * t["µSR"],
  
  p["ISN"] / t["λSN"] + p["IVN"] / t["λVN"] - p["RN"] / t["τ"],
  p["ISR"] / t["λSR"] + p["IVR"] / t["λVR"] - p["RR"] / t["τ"],
  
  t["BSN"] * I * p["SN"] / N - p["CSN"] / t["v"] + t["E"] * (p["CSN_P"] - p["CSN"]),
  t["BVN"] * I * p["VN"] / N - p["CVN"] / t["v"] + t["E"] * (p["CVN_P"] - p["CVN"]),
  t["BVR"] * I * p["VR"] / N - p["CVR"] / t["v"] + t["E"] * (p["CVR_P"] - p["CVR"]),
  t["BSR"] * I * p["SR"] / N - p["CSR"] / t["v"] + t["E"] * (p["CSR_P"] - p["CSR"]),
  
  p["RN"] / t["τ"] + p["SN"] / t["α"] - t["BVN"] * I * p["VN"] / N - t["η"] * p["VN"],
  p["RR"] / t["τ"] + p["SR"] / t["α"] - t["BVR"] * I * p["VR"] / N - t["η"] * p["VR"],
  
  p["ISN"] * t["µSN"] + p["IVN"] * t["µVN"],
  p["ISR"] * t["µSR"] + p["IVR"] * t["µVR"],
    
  #Echange internationaux
  - t["BSN"]* I_P * p["SN_P"] / N_E - p["SN_P"] / t["α"] + t["η"] * p["VN_P"],
  - t["BSR"]* I_P * p["SR_P"] / N_E - p["SR_P"] / t["α"] + t["η"] * p["VR_P"],
  
  p["CSN_P"] / t["v"] - p["ISN_P"] / t["λSN"] - p["ISN_P"] * t["µSN"],
  p["CSR_P"] / t["v"] - p["ISR_P"] / t["λSN"] - p["ISR_P"] * t["µVN"],
  p["CVN_P"] / t["v"] - p["IVN_P"] / t["λSR"] - p["IVN_P"] * t["µVR"],
  p["CVR_P"] / t["v"] - p["IVR_P"] / t["λSR"] - p["IVR_P"] * t["µSR"],
  
  p["ISN_P"] / t["λSN"] + p["IVN_P"] / t["λVN"] - p["RN_P"] / t["τ"],
  p["ISR_P"] / t["λSR"] + p["IVR_P"] / t["λVR"] - p["RR_P"] / t["τ"],
  
  t["BSN"] * I_P * p["SN_P"] / N_E - p["CSN_P"] / t["v"] + t["E"] * (p["CSN"] - p["CSN_P"] ),
  t["BVN"] * I_P * p["VN_P"] / N_E - p["CVN_P"] / t["v"] + t["E"] * (p["CVN"] - p["CVN_P"] ),
  t["BVR"] * I_P * p["VR_P"] / N_E - p["CVR_P"] / t["v"] + t["E"] * (p["CVR"] - p["CVR_P"] ),
  t["BSR"] * I_P * p["SR_P"] / N_E - p["CSR_P"] / t["v"] + t["E"] * (p["CSR"] - p["CSR_P"] ),
  
  p["RN_P"] / t["τ"] + p["SN_P"] / t["α"] - t["BVN"] * I_P * p["VN_P"] / N_E - t["η"] * p["VN_P"],
  p["RR_P"] / t["τ"] + p["SR_P"] / t["α"] - t["BVR"] * I_P * p["VR_P"] / N_E - t["η"] * p["VR_P"],
  
  p["ISN_P"] * t["µSN"] + p["IVN_P"] * t["µVN"],
  p["ISR_P"] * t["µSR"] + p["IVR_P"] * t["µVR"])}

#Scénarios
combi = expand.grid( Vaccin = 0:1,Gestes_Barrieres = 0:1,Confinement = 0:1)
combi$label = paste(combi$Vaccin,combi$Gestes_Barrieres,combi$Confinement,sep = "_")
scenarios = split(combi, f= combi$label)

Scenario_Taux <- function(f,scenario,taux,e){
  if (identical(f,SIR)) B = 1 else B = 1:4
  
  if(scenario[["Vaccin"]]){
    taux[which(taux[,"t"] > e["reactiviteDuGouvernementVaccin"]),"α"] = abs(1/2 * sin(taux[which(taux[,"t"] > e["reactiviteDuGouvernementVaccin"]),"t"])) * taux[which(taux[,"t"] > e["reactiviteDuGouvernementVaccin"]),"α"]
  }
  if(scenario[["Gestes_Barrieres"]]){
    taux[which(taux[,"t"] > e["reactiviteDuGouvernementGestesBarrieres"]),B] = taux[which(taux[,"t"] > e["reactiviteDuGouvernementGestesBarrieres"]),B] *  (1 - e["ReductionInfectionParGestesBarrieresEnPourcent"]) 
  }
  if(scenario[["Confinement"]]){
    taux[which(e["DebutConfinement1"] < taux[,"t"] & taux[,"t"] < e["DebutConfinement1"] + e["DureeConfinement1"]) ,B] = taux[which(e["DebutConfinement1"] < taux[,"t"] & taux[,"t"] < e["DebutConfinement1"] + e["DureeConfinement1"]),B] *  (1 - e["ReductionInfectionParConfinementEnPourcent"])
  }
  taux}

exec_scenario <- function(scenario,f,solveur,taux,pop,effets){
  taux = Scenario_Taux(f, scenario, taux,effets)
  
  for (i in 2:nrow(pop)) {
    pop[i,] <- rk4(f, pop[i-1,], taux[i-1,])}
 pop
  }