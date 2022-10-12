plot_SIR <- function(data){
  # ggplot(SIR_donnees) +
  #   aes(x = periodeEtude, y = Nombre, colour = sous_population) +
  #   geom_point(size = 1.2) +
  #   scale_color_hue(direction = 1) +
  #   facet_grid(Vaccin ~ Gestes_Barrieres + Confinement) +
  #   theme_minimal()
  
  #Font
  windowsFonts(Cabinet = windowsFont("CabinetGrotesk-Extrabold"))
  
  #Labels
  VaccinLabs <- c('TRUE' = "Avec Vaccin",'FALSE' = "Sans Vaccin")
  GestesLabs <- c('TRUE' = "Avec Gestes Barrières", 'FALSE' ="Sans Gestes Barrières")
  ConfiLabs <- c('TRUE' = "Avec Confinement", 'FALSE' ="Sans Confinement")
  
  #Pic de l'épidemie pour chaque scénario
  data_pics <- data %>%
    filter(sous_population == "infectes") %>% 
    group_by(scenario) %>%
    mutate(PicInfection = max(Nombre)) %>% 
    ungroup() %>% 
    filter(Nombre == PicInfection) %>% 
    mutate( PicInfection = paste(format(round(PicInfection/1000),big.mark = ",", trim = TRUE),"K"))

  #ggplot
  p <-   ggplot(data ,
                aes(x = periodeEtude,y= Nombre, group = sous_population, color =str_to_title( sous_population)))+
    geom_line(alpha = 0.9,size = 1.1)+
    geom_hline(yintercept=3*10^6, linetype="dashed", color = "red")+
    facet_grid(Vaccin ~ Gestes_Barrieres + Confinement, labeller = labeller(Vaccin = VaccinLabs, Gestes_Barrieres = GestesLabs, Confinement = ConfiLabs))+
    geom_point(data = data_pics,size = 2.5,show.legend = FALSE)+
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
    ggrepel::geom_text_repel(
      data = data_pics,
      aes(label = PicInfection),
      xlim = c(NA,25),
      ylim = c(4*10^7,NA),
      segment.curvature = .01,
      arrow = arrow(length = unit(.02, "npc"), type = "closed"),
      colour = "black",
      show.legend = FALSE
    )
  return(p)}