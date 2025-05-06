library(shiny)
library(survival)
library(data.table)
library(readxl)
library(survminer)
library(dplyr)
library(plotly)
library(shinyWidgets)

datos<-readRDS("datos.rds")


cox_model <- coxph(Surv(survival_days,evento_muerte)~PSAP+EAP+Filtracion_glomerular+Hemoglobina,data=datos)
HR <- exp(coef(cox_model))

# Intervalos de confianza
IC <- exp(confint(cox_model))

# Combinar en una tabla
tabla_HR <- data.frame(
  Variable = names(HR),
  HR = HR,
  "IC inf" = IC[, 1],
  "IC sup" = IC[, 2],
  "pvalor log rank" = summary(cox_model)$coefficients[,5]
)
tabla_HR$Variable<-c("PSAP <30", "EAP", "Filracion glomerular <30", "Hemoglobina <11.8")

ui <- navbarPage("Supervivencia TAVI",

  tabPanel("Predicciones",
     tags$head(
       tags$style(HTML("
        .input-section {
          height: 10vh;
          display: flex;
          align-items: center;           /* centra horizontalmente */
          justify-content: center;
          margin: 0.5vh;
          margin-left: 30vh;
          margin-right: 30vh;
          padding: 0.5vh;
        }
        .plot-section {
          display: flex;
          flex-direction: column;
          align-items: center;           /* centra horizontalmente */
          justify-content: center;
          height: 60vh;
          overflow: hidden;
          margin: 0.5vh;
          padding: 0.5vh;
          
        }
        .table-section {
          display: flex;
          flex-direction: column;
          align-items: center;           /* centra horizontalmente */
          justify-content: center;
          height: 20vh;
          overflow: auto;
          margin: 0.5vh;
          padding: 0.5vh;
        }
        select {
          width: 100%;
        }
      "))
     ),
    fluidPage(
      div(class = "input-section",
          materialSwitch(inputId = "PSAP", status = "success", label = 'PSAP <30'),
          materialSwitch(inputId = "EAP", status = "danger", label = 'EAP'),
          materialSwitch(inputId = "Filtracion_glomerular", status = "danger", label = 'Filtracion glomerular <30'),
          materialSwitch(inputId = "Hemoglobina", status = "danger", label = 'Hemoglobina <11.8')
      ),
      div(class = "plot-section",
          h3("Funciones supervivencia"),
          plotlyOutput("km_plot", height = "80%", width = "80%")
      ),
      div(class = "table-section",
          h3("Tablas riesgo de fallecimiento"),
          tableOutput("prediccion")
      )
    )
  ),

  tabPanel('Tablas',
    fluidPage(
      fluidRow(
        column(width = 4,
               h4("Tablas coeficientes de riesgo")
        ),
        column(width = 8,
               tableOutput("tabla_riesgos")
        )
      )
    )
  )
)


# Define server logic required to draw a histogram
server <- function(input, output) {
  dato_prediccion <- reactive({
    if(input$PSAP) PSAP = "Si"
    else PSAP = "No"
    if(input$EAP) EAP = "Si"
    else EAP = "No"
    if(input$Filtracion_glomerular) Filtracion_glomerular = "Si"
    else Filtracion_glomerular = "No"
    if(input$Hemoglobina) Hemoglobina = "Si"
    else Hemoglobina = "No"
    data.frame(
      PSAP = PSAP,
      EAP = EAP,
      Filtracion_glomerular = Filtracion_glomerular,
      Hemoglobina= Hemoglobina
    )
  })
  output$tabla_riesgos<-renderTable({
    tabla_HR %>%
      mutate(across(where(is.numeric), ~ formatC(., format = "f", digits = 4)))
  })
  output$prediccion<-renderTable({
      
      surv_fit <- survfit(cox_model, newdata = dato_prediccion())
      km_fit <- survfit(Surv(survival_days, evento_muerte) ~ 1, data = datos) 
      
      tiempos_prediccion <- seq(0, 366, by = 30.44)  
      ajuste_mes_pred <- summary(surv_fit, times = tiempos_prediccion)
      ajuste_mes_datos<- summary(km_fit, times = tiempos_prediccion)
      times <- ajuste_mes_pred$time
      surv_prob <- ajuste_mes_pred$surv
      base_prob <- ajuste_mes_datos$surv
      
      pred_df <- data.frame(
         1-surv_prob,
         1-base_prob
      )
      pred_df_t <- as.data.frame(t(pred_df))
      colnames(pred_df_t) <-as.character(round(times/30.44,0))
      pred_df_t$Meses<-c("Riesgo esperado", "Riesgo base")
      pred_df_t <- pred_df_t[, c("Meses", setdiff(names(pred_df_t), "Meses"))]
      pred_df_t %>%
        mutate(across(where(is.numeric), ~ formatC(., format = "f", digits = 3)))
      pred_df_t
  })
  output$km_plot <- renderPlotly({
    cox_pred <- survfit(cox_model, newdata = dato_prediccion())
    km_fit <- survfit(Surv(survival_days, evento_muerte) ~ 1, data = datos) 
    # Kaplan-Meier en meses
    df_km <- data.frame(
      tiempo_meses = km_fit$time / 30,
      supervivencia = km_fit$surv,
      grupo = "Supervivencia base"
    )
    
    # Cox en meses
    df_cox <- data.frame(
      tiempo_meses = cox_pred$time / 30,
      supervivencia = cox_pred$surv,
      grupo = "Supervivencia predicha"
    )
    
    # Unir las dos curvas
    df_plot <- bind_rows(df_km, df_cox)
    plot <- ggplot(df_plot, aes(x = tiempo_meses, y = supervivencia, color = factor(grupo))) +
      geom_step(size = 1) +
      labs(
        x = "Tiempo (meses)",
        y = "Probabilidad de Supervivencia",
        color = "Curvas"
      ) +
      scale_x_continuous(
        breaks = 0:12,  
        limits = c(0, 12)  
      )+
      ylim(0.5,1)+
      theme_minimal()+
      theme(
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 12),
        legend.position = "top"
      )
    ggplotly(plot) %>%
      layout(
        legend = list(
          orientation = "h",   # horizontal
          x = 0.5,             # centrado
          xanchor = "center",
          y = 1.1              # un poco arriba del gráfico
        ),
        yaxis = list(
          tickmode = "linear",
          dtick = 0.05
        )
      )
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
