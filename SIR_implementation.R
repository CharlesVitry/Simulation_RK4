f_S_v <- function(NbreSains, NbreInfecte,TauxInfection,TauxEfficaciteVaccin) {
  return(-TauxInfection * NbreInfecte * NbreSains -TauxEfficaciteVaccin * NbreSains) 
}

f_I_v <- function(NbreSains, NbreInfecte,TauxInfection,TauxRetablissement) {
  return( (TauxInfection *NbreInfecte * NbreSains) - TauxRetablissement * NbreInfecte)
}

f_R_v <- function(NbreSains, NbreInfecte, TauxRetablissement,TauxEfficaciteVaccin) {
  return(TauxRetablissement * NbreInfecte + TauxEfficaciteVaccin * NbreSains )
}

EDO_RK4_Vaccin <- function(NbreSains, NbreInfecte, NbreRetablis,TauxInfection,TauxRetablissement,TauxEfficaciteVaccin) {
  
  k1_S <- pasDeTemps * f_S_v(NbreSains, NbreInfecte,TauxInfection,TauxEfficaciteVaccin)
  k1_I <- pasDeTemps * f_I(NbreSains, NbreInfecte,TauxInfection,TauxRetablissement)
  k1_R <- pasDeTemps * f_R_v(NbreSains, NbreInfecte, TauxRetablissement,TauxEfficaciteVaccin)
  
  k2_S <- pasDeTemps * f_S_v(NbreSains + (k1_S/2), NbreInfecte + (k1_I/2), TauxInfection, TauxEfficaciteVaccin)
  k2_I <- pasDeTemps * f_I_v(NbreSains + (k1_S/2), NbreInfecte + (k1_I/2), TauxInfection, TauxRetablissement)
  k2_R <- pasDeTemps * f_R_v(NbreSains + (k1_S/2), NbreInfecte + (k1_I/2), TauxRetablissement, TauxEfficaciteVaccin)
  
  k3_S <- pasDeTemps * f_S_v(NbreSains + (k2_S/2), NbreInfecte + (k2_I/2), TauxInfection,TauxEfficaciteVaccin)
  k3_I <- pasDeTemps * f_I_v(NbreSains + (k2_S/2), NbreInfecte + (k2_I/2), TauxInfection,TauxRetablissement)
  k3_R <- pasDeTemps * f_R_v(NbreSains + (k2_S/2), NbreInfecte + (k2_I/2),  TauxRetablissement,TauxEfficaciteVaccin)
  
  k4_S <- pasDeTemps * f_S_v(NbreSains + k3_S, NbreInfecte + k3_I, TauxInfection,TauxEfficaciteVaccin)
  k4_I <- pasDeTemps * f_I_v(NbreSains + k3_S, NbreInfecte + k3_I, TauxInfection,TauxRetablissement)
  k4_R <- pasDeTemps * f_R_v(NbreSains + k3_S, NbreInfecte + k3_I, TauxRetablissement,TauxEfficaciteVaccin)
  
  return(list(
    NbreSains_KPlusUn    = NbreSains    + 1/6 * (k1_S + 2 * k2_S + 2 * k3_S + k4_S) ,
    NbreInfecte_KPlusUn  = NbreInfecte  + 1/6 * (k1_I + 2 * k2_I + 2 * k3_I + k4_I),
    NbreRetablis_KPlusUn = NbreRetablis + 1/6 * (k1_R + 2 * k2_R + 2 * k3_R + k4_R)
  ))
}


