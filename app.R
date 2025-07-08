# Load required libraries
library(shiny)            # Web app framework
library(survival)         # Survival analysis tools
library(data.table)       # Fast data manipulation
library(readxl)           # Read Excel files
library(survminer)        # Enhanced survival plots
library(dplyr)            # Data manipulation
library(plotly)           # Interactive plots
library(shinyWidgets)     # Additional UI widgets
library(shinyjs)          # JavaScript integration for Shiny

# Load pre-calculated percentiles data for risk
risk_and_percentile_in_data <- readRDS("risk_and_percentile_in_data.rds")
risk_and_percentile_in_data <- risk_and_percentile_in_data[order(risk_and_percentile_in_data$percentil), ]

# Add percentile 0 and last risk value (hardcoded 0.4433 as percentile 1 value)
risk_and_percentile_in_data <- as.data.frame(cbind(
  c(0, round(risk_and_percentile_in_data$percentil, digits=4)),
  c(round(risk_and_percentile_in_data$riesgo, digits=4), 0.4433)
))
colnames(risk_and_percentile_in_data) <- c("percentil", "riesgo")

# Load patient data
data <- readRDS("data.rds")

# Fit Cox proportional hazards model with predictors
cox_model <- coxph(Surv(survival_days, evento_muerte) ~ PSAP + EAP + Filtracion_glomerular + Hemoglobina,
                   data = data, ties = "efron")

# Define UI
ui <- fluidPage(
  useShinyjs(),  # Enable shinyjs functions
  
  # JavaScript: detect window width and handle button toggling for plot type
  tags$script(HTML("
    function sendWidth() {
      var w = window.innerWidth || document.documentElement.clientWidth || document.body.clientWidth;
      Shiny.setInputValue('window_width', w, {priority: 'event'});
    }

    $(document).on('shiny:connected', function() {
      sendWidth();
      // Set default active button to Kaplan-Meier
      $('#km_button').addClass('active btn-primary').removeClass('btn-default');
      $('#cox_button').removeClass('active btn-primary').addClass('btn-default');
      Shiny.setInputValue('plot_type', 'Funciones supervivencia', {priority: 'event'});
    });

    $(window).resize(function() { sendWidth(); });

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
  
  # Page title and custom CSS
  tags$head(
    tags$title("TAVI Survival"),
    tags$style(HTML("
      .input-section { min-height: 8vh; display: flex; align-items: center; justify-content: center; gap: 6px; overflow: auto; }
      .plot-section { display: flex; flex-direction: column; align-items: center; justify-content: center; height: 58vh; width: 80%; margin: auto; overflow: auto; }
      .table-section { display: flex; flex-direction: column; align-items: center; justify-content: center; min-height: 10vh; }
      @media (max-width: 768px) {
        .input-section { flex-direction: column; padding: 5px; gap: 0px; }
        .plot-section { height: initial; min-height: 50vh; width: 100%; }
      }
    "))
  ),
  
  # Main UI layout
  titlePanel(div("Predictive Variables for TAVI Survival", style = "text-align: center; margin-top: 2vh;")),
  
  fluidRow(
    div(class = "input-section",
        materialSwitch(inputId = "PSAP", status = "danger", label = 'Pulmonary systolic pressure >=30'),
        materialSwitch(inputId = "EAP", status = "danger", label = 'Peripheral artery disease'),
        materialSwitch(inputId = "Filtracion_glomerular", status = "danger", label = 'Glomerular filtration <30'),
        materialSwitch(inputId = "Hemoglobina", status = "danger", label = 'Hemoglobin <11.8')
    )
  ),
  
  fluidRow(
    div(class = "plot-section",
        div(class = "input-section",  
            actionButton("km_button", "Survival Functions", class = "btn btn-primary active"),
            actionButton("cox_button", "Sample Risk Percentile", class = "btn btn-default")
        ),
        plotlyOutput("km_plot", height = "100%", width = "95%")
    )
  ),
  
  fluidRow(
    div(class = "table-section",
        h3("Predicted Survival Table"),
        tableOutput("prediccion")
    )
  )
)

# Server logic
server <- function(input, output) {
  
  # Reactive block: create prediction data based on user inputs
  dato_prediccion <- reactive({
    data.frame(
      PSAP = ifelse(input$PSAP, "Si", "No"),
      EAP = ifelse(input$EAP, "Si", "No"),
      Filtracion_glomerular = ifelse(input$Filtracion_glomerular, "Si", "No"),
      Hemoglobina = ifelse(input$Hemoglobina, "Si", "No")
    )
  })
  
  # Render survival probability table
  output$prediccion <- renderTable({
    surv_fit <- survfit(cox_model, newdata = dato_prediccion())
    km_fit <- survfit(Surv(survival_days, evento_muerte) ~ 1, data = data)
    
    tiempos_prediccion <- seq(0, 366, by = 30.44)  # Monthly timepoints
    ajuste_mes_pred <- summary(surv_fit, times = tiempos_prediccion)
    ajuste_mes_data <- summary(km_fit, times = tiempos_prediccion)
    
    times <- ajuste_mes_pred$time
    surv_prob <- ajuste_mes_pred$surv
    base_prob <- ajuste_mes_data$surv
    
    # Adjust table format for mobile vs desktop
    if (input$window_width < 768) {
      pred_df <- data.frame(
        as.character(round(times/30.44, 0)),
        surv_prob * 100,
        base_prob * 100
      )
      colnames(pred_df) <- c("Months", "Predicted Survival %", "Sample Survival %")
      
      pred_df[] <- lapply(pred_df, function(x) {
        if (is.numeric(x)) {
          formatted <- formatC(x, format = "f", digits = 1)
        } else x
      })
      pred_df
    } else {
      pred_df <- data.frame(
        surv_prob * 100,
        base_prob * 100
      )
      pred_df_t <- as.data.frame(t(pred_df))
      colnames(pred_df_t) <- as.character(round(times / 30.44, 0))
      pred_df_t$Months <- c("Predicted Survival %", "Sample Survival %")
      pred_df_t <- pred_df_t[, c("Months", setdiff(names(pred_df_t), "Months"))]
      pred_df_t
    }
  })
  
  # Render survival or risk percentile plot
  output$km_plot <- renderPlotly({
    cox_pred <- survfit(cox_model, newdata = dato_prediccion())
    
    if (input$plot_type == "Funciones supervivencia") {
      # Kaplan-Meier vs. Cox plot
      km_fit <- survfit(Surv(survival_days, evento_muerte) ~ 1, data = data)
      df_km <- data.frame(
        tiempo_meses = km_fit$time / 30,
        supervivencia = km_fit$surv,
        group = "Sample Survival"
      )
      df_cox <- data.frame(
        tiempo_meses = cox_pred$time / 30,
        supervivencia = cox_pred$surv,
        group = "Predicted Survival"
      )
      df_plot <- bind_rows(df_km, df_cox)
      df_plot <- df_plot %>%
        mutate(
          tiempo_meses_text = gsub("\\.", ",", formatC(tiempo_meses, format = "f", digits = 1)),
          supervivencia_text = gsub("\\.", ",", formatC(supervivencia, format = "f", digits = 3)),
          tooltip_text = paste0("Time (months): ", tiempo_meses_text, "\nSurvival Probability: ", supervivencia_text)
        )
      
      plot <- ggplot(df_plot, aes(x = tiempo_meses, y = supervivencia, color = group, group = group, text = tooltip_text)) +
        geom_step(aes(linetype = group, alpha = group, size = group)) +
        scale_linetype_manual(name = "K-M curves",values = c("Sample Survival" = "dashed", "Predicted Survival" = "solid")) +
        scale_alpha_manual(name = "K-M curves", values = c("Sample Survival" = 0.4, "Predicted Survival" = 1)) +
        scale_size_manual(name = "K-M curves",values = c("Sample Survival" = 0.7, "Predicted Survival" = 1.2)) +
        scale_color_manual(name = "K-M curves",values = c("Sample Survival" = "gray40", "Predicted Survival" = "steelblue")) +
        labs(x = "Time (months)", y = "Survival Probability") +
        scale_x_continuous(breaks = 0:12, limits = c(0, 12)) +
        ylim(0.5, 1) +
        theme_minimal() +
        theme(
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 12),
          legend.position = "top"
        )
      ggplotly(plot, tooltip="text") %>%
        layout(
          hovermode = "x unified",
          legend = list(
            orientation = "h",
            x = 0.5,
            xanchor = "center",
            y = 1.2
          )
        ) %>%
        config(
          modeBarButtonsToRemove = c(
            "zoom2d", "pan2d", "select2d", "lasso2d",
            "zoomIn2d", "zoomOut2d", "hoverClosestCartesian",
            "hoverCompareCartesian", "toImage"
          ),
          displaylogo = FALSE
        )
    } else {
      # Risk percentile plot
      base_surv_en_t <- summary(cox_pred, times = 365)$surv
      valor_riesgo <- round(1 - base_surv_en_t, 4)
      percentil_riesgo <- approx(x = risk_and_percentile_in_data$riesgo, y = risk_and_percentile_in_data$percentil, xout = valor_riesgo)$y
      percentil_riesgo <- round(percentil_riesgo, 4)
      
      if (percentil_riesgo < 0.996)
        percentil_riesgo <- risk_and_percentile_in_data$percentil[which(risk_and_percentile_in_data$percentil == percentil_riesgo) + 1]
      else
        percentil_riesgo <- 1
      
      plot <- ggplot(risk_and_percentile_in_data, aes(x = percentil, y = riesgo)) +
        geom_step(direction = "hv", color = "steelblue", size = 1) +
        geom_point(aes(text = paste0("Percentile: ", percentil, "\nRisk: ", riesgo)), size = 1, color = "black") +
        geom_vline(xintercept = percentil_riesgo - 0.0001, linetype = "dashed", color = "red", size = 1) +
        labs(x = "Sample death risk percentile", y = "Predicted death risk",
             title = "Predicted risk vs. population percentile") +
        annotate("text", x = percentil_riesgo - 0.05, y = 0.3,
                 label = paste0("Percentile ", round(percentil_riesgo * 100, 2), "%"), color = "red") +
        scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
        theme_minimal()+
        theme(plot.title = element_text(hjust = 0.5, size=15, face = "bold"))
      
      ggplotly(plot, tooltip = "text") %>%
        layout(
          hovermode = "closes"
        )%>%
        config(
          modeBarButtonsToRemove = c(
            "zoom2d", "pan2d", "select2d", "lasso2d",
            "zoomIn2d", "zoomOut2d", "autoScale2d",
            "hoverClosestCartesian",
            "hoverCompareCartesian", "toImage"
          ),
          displaylogo = FALSE
        )
    }
  })
}

# Run the Shiny app
shinyApp(ui = ui, server = server)
