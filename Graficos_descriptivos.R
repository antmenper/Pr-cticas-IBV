grafico1<-function(mix){
  bin_long<-mix%>%
    pivot_longer(
      cols = -c(id, pheno, mic, R, NWT),
      names_to = "variable",
      values_to = "valor")
  bin_summary <- bin_long %>%
    group_by(variable) %>%
    summarise(
      casos = sum(valor == 1, na.rm = TRUE),      
      total = n(),                               
      proporcion = mean(valor, na.rm = TRUE))
  
  ggplot(bin_summary, aes(x = reorder(variable, proporcion), y = proporcion)) +
    geom_col() +
    geom_text(aes(label = casos), hjust = -0.2, size = 3) +   # etiqueta con # de casos
    coord_flip() +
    ylim(0, max(bin_summary$proporcion) + 0.1) +              # espacio para la etiqueta
    labs(x = "Variable", y = "Proporción de 1s")
  
}

grafico2<-function(mix){
  bin_long <- mix %>% 
    pivot_longer(
      cols =  -c(id, pheno, mic, R, NWT),
      names_to = "variable",
      values_to = "valor"
    )
  bin_wide <- mix %>%
    select( -c(id, pheno, mic, R, NWT))  # excluir id
  bin_wide <- bin_wide %>%
    mutate(combinacion = apply(., 1, function(x) {
      paste(names(x)[x == 1], collapse = " + ")
    }))
  bin_wide$combinacion[bin_wide$combinacion == ""] <- "Ninguna"
  combo_summary <- bin_wide %>%
    group_by(combinacion) %>%
    summarise(casos = n()) %>%
    arrange(desc(casos))
  ggplot(combo_summary, aes(x = reorder(combinacion, casos), y = casos / sum(casos))) +
    geom_col(fill = "steelblue") +
    coord_flip() +
    labs(x = "Combinación de variables", y = "Proporción de casos")+theme_minimal()+geom_text(aes(label = paste0(casos, " (", round(casos/sum(casos)*100,1), "%)")),
                                                                                              hjust = -0.1) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15)))
  
}

grafico4<- function(mix){
  # Asegúrate de que mic es numérico
  mix <- mix %>%
    mutate(mic = as.numeric(as.character(mic)))
  
  # Identificar variables categóricas binarias
  vars_categoricas <- mix %>%
    select(-id, -pheno, -mic, -R, -NWT) %>%  # ajusta según tus columnas
    names()
  
  # Crear combinación única para cada perfil
  datos_preparados <- mix %>%
    # Crear identificador de combinación
    unite("combinacion", all_of(vars_categoricas), sep = "-", remove = FALSE) %>%
    # Crear etiqueta con solo las variables presentes
    rowwise() %>%
    mutate(
      vars_presentes = paste(
        vars_categoricas[c_across(all_of(vars_categoricas)) == 1],
        collapse = " + "
      ),
      vars_presentes = ifelse(vars_presentes == "", "No determinant", vars_presentes)
    ) %>%
    ungroup()
  
  # Calcular conteos por combinación y valor de MIC
  datos_grafico <- datos_preparados %>%
    count(mic, vars_presentes) %>%
    arrange(mic, vars_presentes)
  datos_grafico <- datos_grafico %>%
    mutate(mic_fac = factor(mic, levels = c(0.5,1,2,4,8,16,32,64,128,256)))
  
  ggplot(datos_grafico, aes(x = mic_fac, y = n, fill = vars_presentes)) +
    geom_col(position = "stack") +
    scale_x_discrete(
      labels = c("≤2", "4", "", "8", "", "16", "", "≥32", "", "")
    )
  
  # Crear el gráfico de barras apiladas
  ggplot(datos_grafico, aes(x = factor(mic), y = n, fill = vars_presentes)) +
    geom_col(position = "stack", width = 0.7) +
    scale_x_discrete(breaks = unique(datos_grafico$mic)) +
    labs(
      x = "MIC",
      y = "Isolates (n)",
      fill = "Combinación"
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      legend.position = "right",
      legend.text = element_text(size = 8)
    )
  ggplot(datos_grafico, aes(x = mic, y = n, fill = vars_presentes)) +
    geom_col(position = "stack", width = 0.15) +
    scale_x_log10(
      breaks = c(0.5, 1, 2, 4, 8, 16, 32, 64, 128, 256),
      labels = c("≤2", "4", "", "8", "", "16", "", "≥32", "", "")
    ) +
    labs(
      x = "MIC",
      y = "Isolates (n)",
      fill = ""
    ) +
    theme_bw() +
    theme(
      legend.position = "top",
      legend.text = element_text(size = 6),        # Texto más pequeño
      legend.title = element_text(size = 7),       # Título más pequeño
      legend.key.size = unit(0.4, "cm"),           # Keys más pequeñas
      legend.spacing.x = unit(0.1, "cm"),          # Espaciado horizontal
      legend.box.spacing = unit(0.2, "cm"),        # Espaciado de la caja
      panel.grid.minor = element_blank()
    )
}

save(grafico1,grafico2,grafico4,file="graficos_descriptivos.RData")
