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
          min-height: 5vh;
          display: flex;
          align-items: center;           /* centra horizontalmente */
          justify-content: center;
          overflow: auto;
        }
        .plot-section {
          display: flex;
          flex-direction: column;
          align-items: center;           /* centra horizontalmente */
          justify-content: center;
          height: 60vh;
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
            display: inline;
          }
        }
      "))
     ),
     titlePanel(title = div("Variables predictoras de supervivencia TAVI", style = "text-align: center; margin-top: 5vh;margin-bottom 3vh;")),
     fluidRow(
      div(class = "input-section",
          materialSwitch(inputId = "PSAP", status = "success", label = 'Presión sistolica arterial pulmonar <30'),
          materialSwitch(inputId = "EAP", status = "danger", label = 'Enfermedad arterial períferica'),
          materialSwitch(inputId = "Filtracion_glomerular", status = "danger", label = 'Filtracion glomerular <30'),
          materialSwitch(inputId = "Hemoglobina", status = "danger", label = 'Hemoglobina <11.8')
      )
    ),
    fluidRow(
      div(class = "plot-section",
          div(class = "input-section",  
              actionButton("km_button", "Funciones supervivencia", class = "btn btn-primary active"),
              actionButton("cox_button", "Percentil riesgo en los datos", class = "btn btn-default")
          ),
          plotlyOutput("km_plot", height = "100%", width = "95%")
      )
    ),
    fluidRow(
      div(class = "table-section",
          h3("Tablas riesgo de fallecimiento"),
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

      if (input$window_width < 768) {
        pred_df <- data.frame(
          as.character(round(times/30.44,0)),
          1-surv_prob,
          1-base_prob
        )
        colnames(pred_df)<-c("Meses","Riesgo esperado", "Riesgo base")
        pred_df[] <- lapply(pred_df, function(x) {
          if (is.numeric(x)) formatC(x, format = "f", digits = 4) else x
        })
        pred_df
      } else {
        pred_df <- data.frame(
          1-surv_prob,
          1-base_prob
        )

        pred_df_t <- as.data.frame(t(pred_df))
        colnames(pred_df_t) <-as.character(round(times/30.44,0))
        pred_df_t$Meses<-c("Riesgo esperado", "Riesgo base")
        pred_df_t <- pred_df_t[, c("Meses", setdiff(names(pred_df_t), "Meses"))]
        pred_df_t[] <- lapply(pred_df_t, function(x) {
          if (is.numeric(x)) formatC(x, format = "f", digits = 4) else x
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
      plot <- ggplot(df_plot, aes(x = tiempo_meses, y = supervivencia,  color = grupo,
                                  group = grupo,
                                  text = paste0("Tiempo(meses): ", round(tiempo_meses, 1),
                                                "\nProbabilidad supervivencia: ", round(supervivencia, 3)))
                     ) +
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
    ggplotly(plot, tooltip="text") %>%
      layout(
        hoverifo = "x",
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
          "hoverCompareCartesian"
        ),
        displaylogo = FALSE  # Oculta el logo de Plotly
      )
    }
    else{

      base_surv_en_t <- summary(cox_pred, times = 366)$surv
      valor_riesgo <- 1 - base_surv_en_t
      print(valor_riesgo)
      plot<-ggplot(riesgo_percentiles_datos, aes(x = riesgo, y = percentil)) +
        geom_area(fill = "steelblue", size = 1.2) +
        geom_point(size = 2) +
        geom_vline(xintercept = valor_riesgo, linetype = "dashed", color = "red", size = 1) +
        labs(
          x = "Riesgo estimado",
          y = "Percentil"
        ) +
        scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
        theme_minimal()
      ggplotly(plot) %>%
        config(
          modeBarButtonsToRemove = c(
            "zoom2d", "pan2d", "select2d", "lasso2d",
            "zoomIn2d", "zoomOut2d", "autoScale2d",
            "resetScale2d", "hoverClosestCartesian",
            "hoverCompareCartesian", "toImage"
          ),
          displaylogo = FALSE  # Oculta el logo de Plotly
        )
    }
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
