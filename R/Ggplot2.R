library(tidyverse) #Préparation plot
library("rcartocolor") #palette de couleurs spécifique

#Font
windowsFonts(Cabinet = windowsFont("CabinetGrotesk-Extrabold"))

#Labels
VaccinLabs <- c("1" = "Avec Vaccin", "0" = "Sans Vaccin")
GestesLabs <- c("1" = "Avec Gestes Barrières", "0" ="Sans Gestes Barrières")
ConfiLabs <- c("1" = "Avec Confinement", "0" ="Sans Confinement")
PopuLabs <- c("TRUE" = "Population Internationale", "FALSE" ="Population Nationale")

plot_SIR <- function(data, SIRCDV){
  
  #palette de couleurs SIRCDV
  if(SIRCDV){
    palette <- carto_pal(name = "Prism", n = 10)[c(6,6,1,8,8,3,5,4)]
    
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
        SecondePopulation = grepl("international" , sous_population )
      )}else{
        #Transformation données SIR
        data <- data %>%
          pivot_longer(
            cols = c(S,I,R),
            names_to = "sous_population",
            values_to = "Nombre")
      }
  
  #Pic de l'épidemie pour chaque scénario
  data_pics <- data %>%
    filter(sous_population == "infectes") %>% 
    group_by(label) %>%
    mutate(PicInfection = max(Nombre)) %>% 
    ungroup() %>% 
    filter(Nombre == PicInfection) %>% 
    mutate( PicInfection = paste(format(round(PicInfection/1000),big.mark = ",", trim = TRUE),"K"))

  #ggplot
  p <-   ggplot(data ,
                aes(x = t,y= Nombre, group = sous_population, color = str_to_title( sous_population)))+
    geom_line(alpha = 0.9, linewidth = 1.1)+
    geom_hline(yintercept=3*10^6, linetype="dashed", color = "red")+
    geom_point(data = data_pics, size = 2.5, show.legend = FALSE)+
    theme_bw(base_family = "Cabinet",base_size = 15)+
    theme(panel.grid.minor = element_blank(),legend.position = "top",axis.text = element_text(),
          plot.title = ggtext::element_textbox_simple(
            margin = margin(t = 12, b = 12),
            padding = margin(rep(12, 4)),
            fill = "grey90",
            box.color = "grey40",
            r = unit(9, "pt"),
            halign = .5,
            face = "bold",
            lineheight = .9),)+
    scale_y_continuous(
      name = "Evolution des sous-Populations",
      breaks = c(3*10^6,1:10*10^7),
      labels = function(y) paste0( y/10^6,"M"))+
    labs(
      x = "Jours",
      y = "Populations",
      title = "Transmission du COVID dans la population \n Française simulé par le modèle SIR",
      subtitle = "Simulation par résolution RK4",
      color = NULL,
      caption = "Recherche en épidémiologie - IMA")+
    {if(SIRCDV)scale_color_manual(values = palette)}+
    ggrepel::geom_text_repel(
      data = data_pics,
      aes(label = PicInfection),
      xlim = c(NA,25),
      ylim = c(4*10^7,NA),
      segment.curvature = .01,
      arrow = arrow(length = unit(.02, "npc"), type = "closed"),
      colour = "black",
      show.legend = FALSE
    )+{if(SIRCDV)facet_grid(Vaccin + SecondePopulation ~ Gestes_Barrieres + Confinement,labeller = labeller(Vaccin = VaccinLabs,Gestes_Barrieres = GestesLabs,Confinement = ConfiLabs, SecondePopulation = PopuLabs), scales="free_y")}+{if(!SIRCDV)facet_grid(Vaccin ~ Gestes_Barrieres + Confinement, labeller = labeller(Vaccin = VaccinLabs, Gestes_Barrieres = GestesLabs, Confinement = ConfiLabs))}
  return(p)
}
