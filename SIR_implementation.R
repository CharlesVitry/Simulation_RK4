#fonction de scénario

scenario <- function(data){
return( data %>%
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
    TauxInfectionSIR = ifelse(Gestes_Barrieres & periodeEtude > reactiviteDuGouvernementGestesBarrieres , TauxInfectionSIR * (1 - ReductionInfectionParGestesBarrieresEnPourcent)  ,TauxInfectionSIR),
    
    #Effet Confinement
    TauxInfectionSIR = ifelse(Confinement & periodeEtude > DebutConfinement1 & periodeEtude < DebutConfinement1 + DureeConfinement1 , TauxInfectionSIR * (1 - ReductionInfectionParConfinementEnPourcent)  ,TauxInfectionSIR),
    
    #Variable de Scénario issus de la combinaisons des booléens
    scenario = paste (as.character(factor(Vaccin, labels = c("Pas de vaccin", "Avec Vaccin"))),",",
                      as.character( factor(Gestes_Barrieres, labels = c("sans gestes barrièree", "avec gestes barrière"))),"et",
                      as.character( factor(Confinement, labels = c("sans confinement", "avec confinement")) ))
    # ,
    # #Ajout d'aléatoire
    # TauxInfection = abs(jitter(TauxInfection,factor = facteurAleatoire)) ,
    # TauxRetablissement =abs( jitter(TauxRetablissement,factor = facteurAleatoire))
  ))}

################################################
#         SIR modele avec vaccin               #
################################################

f_S_v <- function(S, I,TauxInfectionSIR,TauxEfficaciteVaccin) {
  return(-TauxInfectionSIR * I * S -TauxEfficaciteVaccin * S) 
}

f_I_v <- function(S, I,TauxInfectionSIR,TauxRetablissement) {
  return( (TauxInfectionSIR *I * S) - TauxRetablissement * I)
}

f_R_v <- function(S, I, TauxRetablissement,TauxEfficaciteVaccin) {
  return(TauxRetablissement * I + TauxEfficaciteVaccin * S )
}

EDO_RK4_Vaccin <- function(S, I, R,TauxInfectionSIR,TauxRetablissement,TauxEfficaciteVaccin) {
  
  k1_S <- pasDeTemps * f_S_v(S, I,TauxInfectionSIR,TauxEfficaciteVaccin)
  k1_I <- pasDeTemps * f_I(S, I,TauxInfectionSIR,TauxRetablissement)
  k1_R <- pasDeTemps * f_R_v(S, I, TauxRetablissement,TauxEfficaciteVaccin)
  
  k2_S <- pasDeTemps * f_S_v(S + (k1_S/2), I + (k1_I/2), TauxInfectionSIR, TauxEfficaciteVaccin)
  k2_I <- pasDeTemps * f_I_v(S + (k1_S/2), I + (k1_I/2), TauxInfectionSIR, TauxRetablissement)
  k2_R <- pasDeTemps * f_R_v(S + (k1_S/2), I + (k1_I/2), TauxRetablissement, TauxEfficaciteVaccin)
  
  k3_S <- pasDeTemps * f_S_v(S + (k2_S/2), I + (k2_I/2), TauxInfectionSIR,TauxEfficaciteVaccin)
  k3_I <- pasDeTemps * f_I_v(S + (k2_S/2), I + (k2_I/2), TauxInfectionSIR,TauxRetablissement)
  k3_R <- pasDeTemps * f_R_v(S + (k2_S/2), I + (k2_I/2),  TauxRetablissement,TauxEfficaciteVaccin)
  
  k4_S <- pasDeTemps * f_S_v(S + k3_S, I + k3_I, TauxInfectionSIR,TauxEfficaciteVaccin)
  k4_I <- pasDeTemps * f_I_v(S + k3_S, I + k3_I, TauxInfectionSIR,TauxRetablissement)
  k4_R <- pasDeTemps * f_R_v(S + k3_S, I + k3_I, TauxRetablissement,TauxEfficaciteVaccin)
  
  return(c(
 S    + 1/6 * (k1_S + 2 * k2_S + 2 * k3_S + k4_S) ,
I  + 1/6 * (k1_I + 2 * k2_I + 2 * k3_I + k4_I),
R + 1/6 * (k1_R + 2 * k2_R + 2 * k3_R + k4_R)
  ))
}

################################################
#         SIRCDV*(SIRCDV)' modele              #
################################################

Derivation <- function(SN, SR, ISN, IVN, IVR, ISR, RN, RR,  CSN, CVN, CVR, CSR, VN, VR, DN, DR,SN_P, SR_P, ISN_P, IVN_P, IVR_P, ISR_P, RN_P, RR_P,  CSN_P, CVN_P, CVR_P, CSR_P, VN_P, VR_P, DN_P, DR_P,
                       TauxContamination, TauxInfection, TauxRetablissementSN, TauxRetablissementSR,  TauxVaccination, TauxMortaliteSN, TauxMortaliteVN, TauxMortaliteSR, TauxMortaliteVR, TauxImmunite, TauxEchangeDePopulation
                       ){
  I = ISN + IVN + IVR + ISR
  I_P = ISN_P + IVN_P + IVR_P + ISR_P
  
  #France métropolitaine
  dSN    = - TauxContamination * I * SN - SN / TauxVaccination
  dSR    = - TauxContamination * I * SR - SR/TauxVaccination
  dISN   = CSN / TauxInfection - ISN / TauxRetablissementSN - ISN * TauxMortaliteSN
  dIVN   = CVN / TauxInfection - IVN / TauxRetablissementSN - IVN * TauxMortaliteVN
  dIVR   = CVR / TauxInfection - IVR / TauxRetablissementSR - IVR * TauxMortaliteVR
  dISR   = CSR / TauxInfection - ISR / TauxRetablissementSR - ISR * TauxMortaliteSR
  dRN    = ISN / TauxRetablissementSN + IVN / TauxRetablissementSN - RN / TauxImmunite
  dRR    = ISR / TauxRetablissementSR + IVR / TauxRetablissementSR - RR / TauxImmunite
  dCSN   = TauxContamination * I * SN - CSN / TauxInfection + TauxEchangeDePopulation *( CSN_P - CSN )
  dCVN   = TauxContamination * I * VN - CVN / TauxInfection + TauxEchangeDePopulation *( CVN_P - CVN )
  dCVR   = TauxContamination * I * VR - CVR / TauxInfection + TauxEchangeDePopulation *( CVR_P - CVR )
  dCSR   = TauxContamination * I * SR - CSR / TauxInfection + TauxEchangeDePopulation *( CSR_P - CSR )
  dVN    = RN / TauxImmunite + SN / TauxVaccination - TauxContamination * I * VN
  dVR    = RR / TauxImmunite + SR / TauxVaccination - TauxContamination * I * VR
  dDN    = ISN * TauxMortaliteSN + IVN * TauxMortaliteVN
  dDR    = ISR * TauxMortaliteSR + IVR * TauxMortaliteVR
  
  #Echange internationaux
  dSN_P  = - TauxContamination * I_P * SN_P - SN_P / TauxVaccination
  dSR_P  = - TauxContamination * I_P * SR_P - SR_P/TauxVaccination
  dISN_P = CSN_P / TauxInfection - ISN_P / TauxRetablissementSN - ISN_P * TauxMortaliteSN
  dIVN_P = CVN_P / TauxInfection - IVN_P / TauxRetablissementSN - IVN_P * TauxMortaliteVN
  dIVR_P = CVR_P / TauxInfection - IVR_P / TauxRetablissementSR - IVR_P * TauxMortaliteVR
  dISR_P = CSR_P / TauxInfection - ISR_P / TauxRetablissementSR - ISR_P * TauxMortaliteSR
  dRN_P  = ISN_P / TauxRetablissementSN + IVN_P / TauxRetablissementSN - RN_P / TauxImmunite
  dRR_P  = ISR_P / TauxRetablissementSR + IVR_P / TauxRetablissementSR - RR_P / TauxImmunite
  dCSN_P = TauxContamination * I_P * SN_P - CSN_P / TauxInfection + TauxEchangeDePopulation *( CSN - CSN_P )
  dCVN_P = TauxContamination * I_P * VN_P - CVN_P / TauxInfection + TauxEchangeDePopulation *( CVN - CVN_P )
  dCVR_P = TauxContamination * I_P * VR_P - CVR_P / TauxInfection + TauxEchangeDePopulation *( CVR - CVR_P )
  dCSR_P = TauxContamination * I_P * SR_P - CSR_P / TauxInfection + TauxEchangeDePopulation *( CSR - CSR_P )
  dVN_P  = RN_P / TauxImmunite + SN_P / TauxVaccination - TauxContamination * I_P * VN_P
  dVR_P  = RR_P / TauxImmunite + SR_P / TauxVaccination - TauxContamination * I_P * VR_P
  dDN_P  = ISN_P * TauxMortaliteSN + IVN_P * TauxMortaliteVN
  dDR_P  = ISR_P * TauxMortaliteSR + IVR_P * TauxMortaliteVR
  
return(c("dSN"= dSN ,
         "dSR"= dSR ,
         "dISN"= dISN,
         "dIVN"= dIVN,
         "dIVR"= dIVR,
         "dISR"= dISR,
         "dRN"= dRN ,
         "dRR"= dRR ,
         "dCSN"= dCSN,
         "dCVN"= dCVN,
         "dCVR"= dCVR,
         "dCSR"= dCSR,
         "dVN"= dVN ,
         "dVR"= dVR ,
         "dDN"= dDN ,
         "dDR"= dDR,
         "dSN_P" = dSN_P,
         "dSR_P" = dSR_P,
         "dISN_P" = dISN_P,
         "dIVN_P" = dIVN_P,
         "dIVR_P" = dIVR_P,
         "dISR_P" = dISR_P,
         "dRN_P" = dRN_P,
         "dRR_P" = dRR_P,
         "dCSN_P" = dCSN_P,
         "dCVN_P" = dCVN_P,
         "dCVR_P" = dCVR_P,
         "dCSR_P" = dCSR_P,
         "dVN_P" = dVN_P,
         "dVR_P" = dVR_P,
         "dDN_P" = dDN_P,
         "dDR_P" = dDR_P))}


SIRCDV_RK4 <- function(d){
  k1 = pasDeTemps * Derivation(d$SN, d$SR, d$ISN, d$IVN, d$IVR, d$ISR, d$RN, d$RR,  d$CSN, d$CVN, d$CVR, d$CSR, d$VN, d$VR, d$DN, d$DR,d$SN_P, d$SR_P, d$ISN_P, d$IVN_P, d$IVR_P, d$ISR_P, d$RN_P, d$RR_P,  d$CSN_P, d$CVN_P, d$CVR_P, d$CSR_P, d$VN_P, d$VR_P, d$DN_P, d$DR_P, d$TauxContamination, d$TauxInfection, d$TauxRetablissementSN, d$TauxRetablissementSR,  d$TauxVaccination, d$TauxMortaliteSN, d$TauxMortaliteVN, d$TauxMortaliteSR, d$TauxMortaliteVR, d$TauxImmunite,d$TauxEchangeDePopulation)
  k2 = pasDeTemps * Derivation(d$SN + k1[["dSN"]]/2, d$SR + k1[["dSR"]]/2, d$ISN +k1[["dISN"]]/2, d$IVN +k1[["dIVN"]]/2, d$IVR +k1[["dIVR"]]/2, d$ISR +k1[["dISR"]]/2, d$RN +k1[["dRN"]]/2, d$RR +k1[["dRR"]]/2,  d$CSN +k1[["dCSN"]]/2, d$CVN +k1[["dCVN"]]/2, d$CVR +k1[["dCVR"]]/2, d$CSR +k1[["dCSR"]]/2, d$VN +k1[["dVN"]]/2, d$VR +k1[["dVR"]]/2, d$DN +k1[["dDN"]]/2, d$DR +k1[["dDR"]]/2, d$SN_P + k1[["dSN_P"]]/2, d$SR_P + k1[["dSR_P"]]/2, d$ISN_P +k1[["dISN_P"]]/2, d$IVN_P +k1[["dIVN_P"]]/2, d$IVR_P +k1[["dIVR_P"]]/2, d$ISR_P +k1[["dISR_P"]]/2, d$RN_P +k1[["dRN_P"]]/2, d$RR_P +k1[["dRR_P"]]/2,  d$CSN_P +k1[["dCSN_P"]]/2, d$CVN_P +k1[["dCVN_P"]]/2, d$CVR_P +k1[["dCVR_P"]]/2, d$CSR_P +k1[["dCSR_P"]]/2, d$VN_P +k1[["dVN_P"]]/2, d$VR_P +k1[["dVR_P"]]/2, d$DN_P +k1[["dDN_P"]]/2, d$DR_P +k1[["dDR_P"]]/2, d$TauxContamination, d$TauxInfection, d$TauxRetablissementSN, d$TauxRetablissementSR,  d$TauxVaccination, d$TauxMortaliteSN, d$TauxMortaliteVN, d$TauxMortaliteSR, d$TauxMortaliteVR, d$TauxImmunite,d$TauxEchangeDePopulation)
  k3 = pasDeTemps * Derivation(d$SN + k2[["dSN"]]/2, d$SR + k2[["dSR"]]/2, d$ISN +k2[["dISN"]]/2, d$IVN +k2[["dIVN"]]/2, d$IVR +k2[["dIVR"]]/2, d$ISR +k2[["dISR"]]/2, d$RN +k2[["dRN"]]/2, d$RR +k2[["dRR"]]/2,  d$CSN +k2[["dCSN"]]/2, d$CVN +k2[["dCVN"]]/2, d$CVR +k2[["dCVR"]]/2, d$CSR +k2[["dCSR"]]/2, d$VN +k2[["dVN"]]/2, d$VR +k2[["dVR"]]/2, d$DN +k2[["dDN"]]/2, d$DR +k2[["dDR"]]/2,d$SN_P + k2[["dSN_P"]]/2, d$SR_P + k2[["dSR_P"]]/2, d$ISN_P +k2[["dISN_P"]]/2, d$IVN_P +k2[["dIVN_P"]]/2, d$IVR_P +k2[["dIVR_P"]]/2, d$ISR_P +k2[["dISR_P"]]/2, d$RN_P +k2[["dRN_P"]]/2, d$RR_P +k2[["dRR_P"]]/2, d$CSN_P +k2[["dCSN_P"]]/2, d$CVN_P +k2[["dCVN_P"]]/2, d$CVR_P +k2[["dCVR_P"]]/2, d$CSR_P +k2[["dCSR_P"]]/2, d$VN_P +k2[["dVN_P"]]/2, d$VR_P +k2[["dVR_P"]]/2, d$DN_P +k2[["dDN_P"]]/2, d$DR_P +k2[["dDR_P"]]/2, d$TauxContamination, d$TauxInfection, d$TauxRetablissementSN, d$TauxRetablissementSR,  d$TauxVaccination, d$TauxMortaliteSN, d$TauxMortaliteVN, d$TauxMortaliteSR, d$TauxMortaliteVR, d$TauxImmunite,d$TauxEchangeDePopulation)
  k4 = pasDeTemps * Derivation(d$SN + k3[["dSN"]], d$SR + k3[["dSR"]], d$ISN +k3[["dISN"]], d$IVN +k3[["dIVN"]], d$IVR +k3[["dIVR"]], d$ISR +k3[["dISR"]], d$RN +k3[["dRN"]], d$RR +k3[["dRR"]],  d$CSN +k3[["dCSN"]], d$CVN +k3[["dCVN"]], d$CVR +k3[["dCVR"]], d$CSR +k3[["dCSR"]], d$VN +k3[["dVN"]], d$VR +k3[["dVR"]], d$DN +k3[["dDN"]], d$DR +k3[["dDR"]],d$SN_P + k3[["dSN_P"]], d$SR_P + k3[["dSR_P"]], d$ISN_P +k3[["dISN_P"]], d$IVN_P +k3[["dIVN_P"]], d$IVR_P +k3[["dIVR_P"]], d$ISR_P +k3[["dISR_P"]], d$RN_P +k3[["dRN_P"]], d$RR_P +k3[["dRR_P"]],  d$CSN_P +k3[["dCSN_P"]], d$CVN_P +k3[["dCVN_P"]], d$CVR_P +k3[["dCVR_P"]], d$CSR_P +k3[["dCSR_P"]], d$VN_P +k3[["dVN_P"]], d$VR_P +k3[["dVR_P"]], d$DN_P +k3[["dDN_P"]], d$DR_P +k3[["dDR_P"]], d$TauxContamination, d$TauxInfection, d$TauxRetablissementSN, d$TauxRetablissementSR,  d$TauxVaccination, d$TauxMortaliteSN, d$TauxMortaliteVN, d$TauxMortaliteSR, d$TauxMortaliteVR, d$TauxImmunite,d$TauxEchangeDePopulation) 

  resultat = c()
  
  for (i in colnames(d[12:43])){
    a = paste("d",i,sep="")
    resultat = append(resultat, d[[i]]  + 1/6 * (k1[[a]] + 2 * k2[[a]] + 2 * k3[[a]] + k4[[a]]))
  }
  return(resultat)
}