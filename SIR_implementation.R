f_S <- function(NbreSains, NbreInfecte,TauxInfection) {
  return(-TauxInfection * NbreInfecte * NbreSains / tailleDeLaPopulationTotale)
}

f_I <- function(NbreSains, NbreInfecte,TauxInfection,TauxRetablissement) {
  return( (TauxInfection *NbreInfecte * NbreSains / tailleDeLaPopulationTotale) - TauxRetablissement * NbreInfecte)
}

f_R <- function(NbreInfecte,TauxRetablissement) {
  return(TauxRetablissement * NbreInfecte)
}

EDO_RK4 <- function(NbreSains, NbreInfecte, NbreRetablis,TauxInfection,TauxRetablissement) {
  
  k1_S <- f_S(NbreSains, NbreInfecte,TauxInfection)
  k1_I <- f_I(NbreSains, NbreInfecte,TauxInfection,TauxRetablissement)
  k1_R <- f_R(NbreInfecte,TauxRetablissement)
  
  k2_S <- f_S( NbreSains + pasDeTemps / 2 * k1_S, NbreInfecte + pasDeTemps / 2 * k1_I,TauxInfection)
  k2_I <- f_I( NbreSains + pasDeTemps / 2 * k1_S, NbreInfecte + pasDeTemps / 2 * k1_I,TauxInfection,TauxRetablissement)
  k2_R <- f_R( NbreInfecte + pasDeTemps / 2 * k1_I,TauxRetablissement)
  
  k3_S <- f_S(NbreSains + pasDeTemps / 2 * k2_S, NbreInfecte + pasDeTemps / 2 * k2_I,TauxInfection)
  k3_I <- f_I(NbreSains + pasDeTemps / 2 * k2_S, NbreInfecte + pasDeTemps / 2 * k2_I,TauxInfection,TauxRetablissement)
  k3_R <- f_R(NbreInfecte + pasDeTemps / 2 * k2_I,TauxRetablissement)
  
  k4_S <- f_S(NbreSains + pasDeTemps * k3_S, NbreInfecte + pasDeTemps * k3_I,TauxInfection)
  k4_I <- f_I(NbreSains + pasDeTemps * k3_S, NbreInfecte + pasDeTemps * k3_I,TauxInfection,TauxRetablissement)
  k4_R <- f_R(NbreInfecte + pasDeTemps * k3_I,TauxRetablissement)
  
  return(list(
    NbreSains_PlusUn = NbreSains + pasDeTemps / 6 * (k1_S + 2 * k2_S + 2 * k3_S + k4_S),
    NbreInfecte_PlusUn = NbreInfecte + pasDeTemps / 6 * (k1_I + 2 * k2_I + 2 * k3_I + k4_I),
    NbreRetablis_PlusUn = NbreRetablis + pasDeTemps / 6 * (k1_R + 2 * k2_R + 2 * k3_R + k4_R)
  ))
}
