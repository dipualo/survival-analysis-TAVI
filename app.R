library(shiny)
library(survival)
library(data.table)
library(readxl)
library(survminer)
library(dplyr)
library(plotly)
library(shinyWidgets)
library(shinyjs)

riesgo_percentiles_datos<-readRDS("riesgo_percentiles_datos.rds")
riesgo_percentiles_datos <- riesgo_percentiles_datos[order(riesgo_percentiles_datos$percentil), ]
riesgo_percentiles_datos<- rbind(data.frame(percentil = 0, riesgo = 0), riesgo_percentiles_datos)

datos<-readRDS("datos.rds")

cox_model <- coxph(Surv(survival_days,evento_muerte)~PSAP+EAP+Filtracion_glomerular+Hemoglobina,data=datos, ties="efron")

ui <- fluidPage(
      useShinyjs(),
      tags$script(HTML("
        function sendWidth() {
          var w = window.innerWidth || document.documentElement.clientWidth || document.body.clientWidth;
          Shiny.setInputValue('window_width', w, {priority: 'event'});
        }
    
        // Ejecutar al cargar y al redimensionar
        $(document).on('shiny:connected', function() {
          sendWidth();
        });
    
        $(window).resize(function() {
          sendWidth();
        });  
        $(document).on('shiny:connected', function() {
          // Inicializar con Kaplan-Meier como activo por defecto
          $('#km_button').addClass('active btn-primary').removeClass('btn-default');
          $('#cox_button').removeClass('active btn-primary').addClass('btn-default');
          Shiny.setInputValue('plot_type', 'Funciones supervivencia', {priority: 'event'});
        });
      
        $(document).on('click', '#km_button', function() {
          Shiny.setInputValue('plot_type', 'Funciones supervivencia', {priority: 'event'});
          $('#km_button').addClass('active btn-primary').removeClass('btn-default');
          $('#cox_button').removeClass('active btn-primary').addClass('btn-default');
        });
      
        $(document).on('click', '#cox_button', function() {
          Shiny.setInputValue('plot_type', 'Cox', {priority: 'event'});
          $('#cox_button').addClass('active btn-primary').removeClass('btn-default');
          $('#km_button').removeClass('active btn-primary').addClass('btn-default');
        });
      ")),
     tags$head(
       tags$style(HTML("
        .input-section {
          min-height: 10vh;
          display: flex;
          align-items: center;           /* centra horizontalmente */
          justify-content: center;
          gap: 10px; 
          overflow: auto;
        }
        .plot-section {
          display: flex;
          flex-direction: column;
          align-items: center;           /* centra horizontalmente */
          justify-content: center;
          height: 55vh;
          width: 80%;
          margin-left: auto;
          margin-right: auto;
          overflow: auto;
        }
        .table-section {
          display: flex;
          flex-direction: column;
          align-items: center;           /* centra horizontalmente */
          justify-content: center;
          min-height: 20vh;
        }
        select {
          width: 100%;
        }
        @media (max-width: 768px) {
          .input-section {
              display: flex;
              flex-direction: column; /* Pone los elementos en columna */
              margin:5px;
          }
          .plot-section{
            width:100%
          }
        }
      "))
     ),
     titlePanel(title = div("Variables predictoras de supervivencia TAVI", style = "text-align: center; margin-top: 3vh;margin-bottom 6vh; ")),
     fluidRow(
      div(class = "input-section",
          materialSwitch(inputId = "PSAP", status = "danger", label = 'Presión sistolica arterial pulmonar >=30'),
          materialSwitch(inputId = "EAP", status = "danger", label = 'Enfermedad arterial períferica'),
          materialSwitch(inputId = "Filtracion_glomerular", status = "danger", label = 'Filtracion glomerular <30'),
          materialSwitch(inputId = "Hemoglobina", status = "danger", label = 'Hemoglobina <11.8')
      )
    ),
    fluidRow(
      div(class = "plot-section",
          div(class = "input-section",  
              actionButton("km_button", "Funciones supervivencia", class = "btn btn-primary active"),
              actionButton("cox_button", "Percentil riesgo poblacional", class = "btn btn-default")
          ),
          plotlyOutput("km_plot", height = "100%", width = "95%")
      )
    ),
    fluidRow(
      div(class = "table-section",
          h3("Tabla probabilidad de supervivencia"),
          tableOutput("prediccion")
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
  
  output$prediccion<-renderTable({
    
    surv_fit <- survfit(cox_model, newdata = dato_prediccion())
    km_fit <- survfit(Surv(survival_days, evento_muerte) ~ 1, data = datos) 
    
    tiempos_prediccion <- seq(0, 366, by = 30.44)  
    ajuste_mes_pred <- summary(surv_fit, times = tiempos_prediccion)
    ajuste_mes_datos<- summary(km_fit, times = tiempos_prediccion)
    times <- ajuste_mes_pred$time
    surv_prob <- ajuste_mes_pred$surv
    base_prob <- ajuste_mes_datos$surv
    
    if (input$window_width < 768) {
      pred_df <- data.frame(
        as.character(round(times/30.44,0)),
        surv_prob*100,
        base_prob*100
      )
      colnames(pred_df)<-c("Meses","Supervivencia predicha %", "Supervivencia poblacional %")
      pred_df[] <- lapply(pred_df, function(x) {
        if (is.numeric(x)) formatC(x, format = "f", digits = 1) else x
      })
      pred_df
    } else {
      pred_df <- data.frame(
        surv_prob*100,
        base_prob*100
      )
      
      pred_df_t <- as.data.frame(t(pred_df))
      colnames(pred_df_t) <-as.character(round(times/30.44,0))
      pred_df_t$Meses<-c("Supervivencia predicha %", "Supervivencia poblacional %")
      pred_df_t <- pred_df_t[, c("Meses", setdiff(names(pred_df_t), "Meses"))]
      pred_df_t[] <- lapply(pred_df_t, function(x) {
        if (is.numeric(x)) formatC(x, format = "f", digits = 1) else x
      })
      pred_df_t
    }
  })
  output$km_plot <- renderPlotly({
    
    cox_pred <- survfit(cox_model, newdata = dato_prediccion())
    
    if(input$plot_type == "Funciones supervivencia") {
      km_fit <- survfit(Surv(survival_days, evento_muerte) ~ 1, data = datos) 
      # Kaplan-Meier en meses
      df_km <- data.frame(
        tiempo_meses = km_fit$time / 30,
        supervivencia = km_fit$surv,
        grupo = "Supervivencia poblacional"
      )
      
      # Cox en meses
      df_cox <- data.frame(
        tiempo_meses = cox_pred$time / 30,
        supervivencia = cox_pred$surv,
        grupo = "Supervivencia predicha"
      )
      
      # Unir las dos curvas
      df_plot <- bind_rows(df_km, df_cox)
      plot <- ggplot(df_plot, aes(x = tiempo_meses, y = supervivencia,  color = grupo,
                                  group = grupo,
                                  text = paste0("Tiempo(meses): ", round(tiempo_meses, 1),
                                                "\nProbabilidad supervivencia: ", round(supervivencia, 3)))
                     ) +
      geom_step(size = 1) +
      labs(
        x = "Tiempo (meses)",
        y = "Probabilidad de Supervivencia",
        color = "Curvas K-M"
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
    ggplotly(plot, tooltip="text") %>%
      layout(
        hovermode="x unified",
        legend = list(
          orientation = "h",   # horizontal
          x = 0.5,             # centrado
          xanchor = "center",
          y = 1.2              # un poco arriba del gráfico
        ),
        yaxis = list(
          tickmode = "linear",
          dtick = 0.05
        )
      ) %>%
      config(
        modeBarButtonsToRemove = c(
          "zoom2d", "pan2d", "select2d", "lasso2d",
          "zoomIn2d", "zoomOut2d", "hoverClosestCartesian",
          "hoverCompareCartesian","toImage"
        ),
        displaylogo = FALSE  # Oculta el logo de Plotly
      )
    }
    else{
      
      riesgos <- predict(cox_model, type = "risk")
      riesgo_ind <- predict(cox_model, newdata = dato_prediccion(), type = "risk")
      percentil_ind <- ecdf(riesgos)(riesgo_ind) * 100
      
      # Crear data frame para ggplot
      df_riesgos <- data.frame(riesgo = riesgos) %>%
        mutate(percentil = percent_rank(riesgo) * 100) %>%  # Percentil empírico
        arrange(percentil)
      
      # Plot con ggplot2
      plot<-ggplot(df_riesgos, aes(x = percentil, y = riesgo)) +
        geom_step(direction = "hv", color = "steelblue", size = 1) +
        geom_vline(xintercept = percentil_ind, color = "red", linetype = "dashed", linewidth = 1) +
        labs(
          title = "Riesgo relativo según percentil de riesgo de la población",
          x = "Percentil de riesgo (%)",
          y = "Riesgo relativo (exp(Xβ))"
        ) +
        theme_minimal()+
        theme(
          plot.title = element_text(hjust = 0.5, size = 18),  
        )
       
      
      g<-ggplotly(plot) %>%
        layout(
          hovermode="x unified",
          legend = list(
            orientation = "h",   # horizontal
            x = 10,             # centrado
            xanchor = "center",
            y = 1.2              # un poco arriba del gráfico
          ),
          yaxis = list(
            tickmode = "linear",
            dtick = 1
          ),
          xaxis = list(
            tickmode = "linear",
            dtick = 10
          )
        ) %>%
        config(
          modeBarButtonsToRemove = c(
            "zoom2d", "pan2d", "select2d", "lasso2d",
            "zoomIn2d", "zoomOut2d", "autoScale2d",
            "hoverClosestCartesian",
            "hoverCompareCartesian", "toImage"
          ),
          displaylogo = FALSE  # Oculta el logo de Plotly
        )
        g<-style(g, hoverinfo = "none", traces = 2)
        g
    }
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
