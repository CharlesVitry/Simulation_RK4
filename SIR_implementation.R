f_S_v <- function(S, I,TauxInfection,TauxEfficaciteVaccin) {
  return(-TauxInfection * I * S -TauxEfficaciteVaccin * S) 
}

f_I_v <- function(S, I,TauxInfection,TauxRetablissement) {
  return( (TauxInfection *I * S) - TauxRetablissement * I)
}

f_R_v <- function(S, I, TauxRetablissement,TauxEfficaciteVaccin) {
  return(TauxRetablissement * I + TauxEfficaciteVaccin * S )
}

EDO_RK4_Vaccin <- function(S, I, R,TauxInfection,TauxRetablissement,TauxEfficaciteVaccin) {
  
  k1_S <- pasDeTemps * f_S_v(S, I,TauxInfection,TauxEfficaciteVaccin)
  k1_I <- pasDeTemps * f_I(S, I,TauxInfection,TauxRetablissement)
  k1_R <- pasDeTemps * f_R_v(S, I, TauxRetablissement,TauxEfficaciteVaccin)
  
  k2_S <- pasDeTemps * f_S_v(S + (k1_S/2), I + (k1_I/2), TauxInfection, TauxEfficaciteVaccin)
  k2_I <- pasDeTemps * f_I_v(S + (k1_S/2), I + (k1_I/2), TauxInfection, TauxRetablissement)
  k2_R <- pasDeTemps * f_R_v(S + (k1_S/2), I + (k1_I/2), TauxRetablissement, TauxEfficaciteVaccin)
  
  k3_S <- pasDeTemps * f_S_v(S + (k2_S/2), I + (k2_I/2), TauxInfection,TauxEfficaciteVaccin)
  k3_I <- pasDeTemps * f_I_v(S + (k2_S/2), I + (k2_I/2), TauxInfection,TauxRetablissement)
  k3_R <- pasDeTemps * f_R_v(S + (k2_S/2), I + (k2_I/2),  TauxRetablissement,TauxEfficaciteVaccin)
  
  k4_S <- pasDeTemps * f_S_v(S + k3_S, I + k3_I, TauxInfection,TauxEfficaciteVaccin)
  k4_I <- pasDeTemps * f_I_v(S + k3_S, I + k3_I, TauxInfection,TauxRetablissement)
  k4_R <- pasDeTemps * f_R_v(S + k3_S, I + k3_I, TauxRetablissement,TauxEfficaciteVaccin)
  
  return(list(
    NbreSains_KPlusUn    = S    + 1/6 * (k1_S + 2 * k2_S + 2 * k3_S + k4_S) ,
    NbreInfecte_KPlusUn  = I  + 1/6 * (k1_I + 2 * k2_I + 2 * k3_I + k4_I),
    NbreRetablis_KPlusUn = R + 1/6 * (k1_R + 2 * k2_R + 2 * k3_R + k4_R)
  ))
}


