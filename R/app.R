source("Solve.R") # Résolution RK4 des équas diff
source("Ggplot2.R") # plot des calculs
library(plotly) # passage du plot en interactif
library(shiny) # serveur

taux_SIR = c(v = 0.5, λ = 3,  α = 0.0007)
taux_SIR = cbind( t(replicate(nb_iter, taux_SIR)), t)

taux_SIRCDV = c(BSN = 0.4, BSR = 0.6,  BVN = 0.6, BVR = 0.6,
                v = 5,
                λSN = 5, λSR = 5, λVN = 5, λVR = 5,
                α = 500,
                µSN = 0.000001, µVN = 0.000001, µVR = 0.000001, µSR = 0.000001,
                τ = 10,
                E = 0.9,
                η = 0.2)
taux_SIRCDV = cbind( t(replicate(nb_iter, taux_SIRCDV)), t)

# côté client
ui <- fluidPage(
  titlePanel("Simulation par résolution RK4 de la transmission du COVID\n dans la population national et international par scénarios"),
  sidebarLayout(sidebarPanel(
    #Inputs
    selectInput("f","Modèle : ",c("SIR","SIRCDV")),br(),
    sliderInput("reactiviteDuGouvernementVaccin","Réactivité du gouvernement sur la campagne vaccinale",value= 20, min=0 ,max=tempsEtude ),
    sliderInput("reactiviteDuGouvernementGestesBarrieres","Réactivité du gouvernement sur les gestes barrières",value= 10, min=0 ,max=tempsEtude),
    sliderInput("ReductionInfectionParGestesBarrieresEnPourcent","Réduction de l'infection par les gestes barrières",value= 0.2, min= 0,max=0.4),
    sliderInput("ReductionInfectionParConfinementEnPourcent","Réduction de l'infection par le confinement",value= 0.95, min=0.5 ,max=0.99),
    sliderInput("DebutConfinement1","Jour du début du confinement",value=30 , min=0 ,max=tempsEtude),
    sliderInput("DureeConfinement1","Durée en jours du confinement",value=40 ,min=10 ,max=tempsEtude),br()),
    
    #Main
    mainPanel(tabsetPanel(type = "tabs",
                          tabPanel("Plot", plotlyOutput('plot')),
                          tabPanel("Stats", verbatimTextOutput("stats")),
                          tabPanel("Info", tableOutput("info"))))))

# côté serveur
server <- function(input, output) {
  
  #données
  d <- reactive({
    effets = c(
      "reactiviteDuGouvernementVaccin" = input$reactiviteDuGouvernementVaccin,
      "reactiviteDuGouvernementGestesBarrieres" = input$reactiviteDuGouvernementGestesBarrieres,
      "ReductionInfectionParGestesBarrieresEnPourcent" = input$ReductionInfectionParGestesBarrieresEnPourcent,
      "ReductionInfectionParConfinementEnPourcent" = input$ReductionInfectionParConfinementEnPourcent,
      "DebutConfinement1" = input$DebutConfinement1, 
      "DureeConfinement1" = input$DureeConfinement1)
    
    if(input$f == "SIR"){f = SIR ; taux = taux_SIR; pop = pop_SIR}
    else{f = SIRCDV ; taux = taux_SIRCDV ; pop = pop_SIRCDV}
    
    pops = lapply(scenarios, exec_scenario, f= f, taux=taux, pop=pop, effets = effets )
    d =  data.frame(do.call("rbind",pops),t,combi[rep(seq_len(nrow(combi)),each = nb_iter),])})
  
  # plot
  output$plot <- renderPlotly({ 
    ggplotly(plot(d(),input$f == 'SIRCDV'))})
  
  # stats
  output$stats <- renderPrint({
    summary(d())
  })
  
  # Infos
  output$info <- renderTable({
    d()
  })}

#Lancement serveur
shinyApp(ui, server)