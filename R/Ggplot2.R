library(tidyverse) #Préparation plot
library("rcartocolor") #palette de couleurs spécifique

#Font
windowsFonts(Cabinet = windowsFont("CabinetGrotesk-Extrabold"))

#Labels
VaccinLabs <- c("1" = "Avec Vaccin", "0" = "Sans Vaccin")
GestesLabs <- c("1" = "Avec Gestes Barrières", "0" ="Sans Gestes Barrières")
ConfiLabs <- c("1" = "Avec Confinement", "0" ="Sans Confinement")
PopuLabs <- c("TRUE" = "Population Internationale", "FALSE" ="Population Nationale")

plot <- function(data, SIRCDV){

if(SIRCDV){
#palette de couleurs SIRCDV
palette <- carto_pal(name = "Prism", n = 10)[c(6,6,1,8,8,3,5,4)]
hauteur = 1200
#Transo données SIRCDV 
data <- data %>%
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
    SecondePopulation = grepl("international" , sous_population ))}
else{
    hauteur = 900
    #Transformation données SIR
    data <- data %>%
      pivot_longer(
        cols = c(S,I,R),
        names_to = "sous_population",
        values_to = "Nombre")}
  
  #ggplot
  p <- ggplot(data, aes(x = t,y= Nombre, group = sous_population, color = str_to_title( sous_population)))+
    geom_line(alpha = 0.9, linewidth = 1.1)+
    geom_hline(yintercept=3*10^6, linetype="dashed", color = "red")+
    guides(color = FALSE)+
    theme_minimal(base_family = "Cabinet", base_size = 18)+
    theme(strip.placement = "outside", plot.margin = margin(t = 120))+
    scale_y_continuous(
      name = "Evolution des sous-Populations",
      breaks = c(3*10^6,1:10*10^7),
      labels = function(y) paste0( y/10^6,"M"))+
    labs(
      x = "Jours",
      y = "Populations",
      color = NULL)+
    {if(SIRCDV)scale_color_manual(values = palette)}+
    {if(SIRCDV)facet_grid(Vaccin + SecondePopulation ~ Gestes_Barrieres + Confinement,labeller = labeller(Vaccin = VaccinLabs,Gestes_Barrieres = GestesLabs,Confinement = ConfiLabs, SecondePopulation = PopuLabs))}+
    {if(!SIRCDV)facet_grid(Vaccin ~ Gestes_Barrieres + Confinement, labeller = labeller(Vaccin = VaccinLabs, Gestes_Barrieres = GestesLabs, Confinement = ConfiLabs))}

  ggplotly(p, tooltip = c("Nombre","sous_population"), height = hauteur)}
