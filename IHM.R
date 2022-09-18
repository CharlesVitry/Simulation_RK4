plot_SIR <- function(data){
  #Font
  windowsFonts(Cabinet = windowsFont("CabinetGrotesk-Extrabold"))
  
  #Labels
  VaccinLabs <- c('TRUE' = "Avec Vaccin",'FALSE' = "Sans Vaccin")
  PolitiLabs <- c('TRUE' = "Avec Politique Sanitaire", 'FALSE' ="Sans Politique Sanitaire")
  
  #Pic de l'épidemie pour chaque scénario
  data_pics <- data %>%
    filter(sous_population == "infectes") %>% 
    group_by(scenario) %>%
    mutate(PicInfection = max(Nombre)) %>% 
    ungroup() %>% 
    filter(Nombre == PicInfection)
  
  
  #ggplot
  p <-   ggplot(data ,
                aes(x = periodeEtude,y= Nombre, group = sous_population, color =str_to_title( sous_population)))+
    geom_line(alpha = 0.9,
              size = 1.1
    )+
    geom_hline(yintercept=3*10^6, linetype="dashed", color = "red")+
    geom_point(
      data = data_pics,
      size = 2.5
    )+
    ggrepel::geom_text_repel(
      data = data_pics,
      aes(label = round(PicInfection)),
      xlim = c(NA,25),
      ylim = c(4*10^7,NA),
      segment.curvature = .01,
      arrow = arrow(length = unit(.02, "npc"), type = "closed")
    )+
    facet_grid(Vaccin ~ PolitiqueSanitaire , labeller = labeller(Vaccin = VaccinLabs, PolitiqueSanitaire = PolitiLabs))+
    theme_bw(
      base_family = "Cabinet",
      base_size = 15
    )+
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "top",
      axis.text = element_text(
        # angle = 20
      ),
      plot.title = ggtext::element_textbox_simple(
        margin = margin(t = 12, b = 12),
        padding = margin(rep(12, 4)),
        fill = "grey90",
        box.color = "grey40",
        r = unit(9, "pt"),
        halign = .5,
        face = "bold",
        lineheight = .9
      ),
    )+
    scale_y_continuous(
      name = "Evolution des sous-Populations",
      breaks = 0:10*10^7,
      labels = function(y) paste0( y/10^6,"M")
    )+
    labs(
      x = "Jours",
      y = "Populations",
      title = "Transmission du COVID dans la population \n Française simulé par le modèle SIR",
      subtitle = "Simulation par résolution RK4",
      color = NULL,
      caption = "Recherche en épidémiologie - IMA"
    )
  return(p)
}