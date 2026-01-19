analizar_antibiotico <- function(antibiotic, 
                                 drug_class,
                                 geno_table = geno2,
                                 pheno_table = y_largo,
                                 geno_sample_col = "Name",
                                 pheno_sample_col = "biosample",
                                 sir_col = "pheno",
                                 marker_col = "marker",
                                 maf = 10,
                                 min_set_size = 5,
                                 mostrar_graficos = TRUE,
                                 verbose = TRUE) {
  
  # Mensaje informativo
  if(verbose) {
    cat("\n========================================\n")
    cat("Análisis para:", antibiotic, "\n")
    cat("Clase de droga:", drug_class, "\n")
    cat("========================================\n\n")
  }
  
  # 1. Crear matriz binaria
  if(verbose) cat("1. Generando matriz binaria...\n")
  
  mix_data <- get_binary_matrix(
    geno_table,
    pheno_table,
    antibiotic = antibiotic,
    drug_class_list = drug_class,
    geno_sample_col = geno_sample_col,
    pheno_sample_col = pheno_sample_col,
    keep_assay_values = TRUE,
    keep_assay_values_from = "mic",
    sir_col = sir_col,
    ecoff = NULL
  )
  
  # 2. Análisis de resistentes
  if(verbose) cat("\n2. Análisis de resistentes...\n")
  
  res_resistentes <- resistentes(mix_data)
  print(res_resistentes)
  
  if(mostrar_graficos) {
    cat("\n")
    g <- grafico_resistentes(mix_data)
    if(!is.null(g)) print(g)
  }
  
  # 3. Análisis NWT
  if(verbose) cat("\n3. Análisis NWT...\n")
  
  res_nwt <- NWT(mix_data)
  print(res_nwt)
  
  if(mostrar_graficos) {
    cat("\n")
    g <- grafico_NWT(mix_data)
    if(!is.null(g)) print(g)
  }
  
  # 4. Histograma MIC
  if(mostrar_graficos) {
    cat("\n4. Histograma MIC...\n")
    hist(log2(mix_data$mic), 
         main = paste("Distribución MIC -", antibiotic),
         xlab = expression(log[2](MIC)))
  }
  
  # 5. Análisis MIC
  if(verbose) cat("\n5. Análisis MIC...\n")
  
  res_mic <- Mic(mix_data)
  print(res_mic)
  
  # 6. Mutaciones
  if(verbose) cat("\n6. Análisis de mutaciones...\n")
  
  mut_r <- mutaciones_R(mix_data)
  cat("\nMutaciones en Resistentes:\n")
  print(mut_r)
  
  mut_nwt <- mutaciones_NWT(mix_data)
  cat("\nMutaciones en NWT:\n")
  print(mut_nwt)
  
  # 7. Gráficos adicionales
  if(mostrar_graficos) {
    cat("\n7. Gráficos adicionales...\n")
    
    cat("\nGráfico 1:\n")
    g1 <- grafico1(mix_data)
    if(!is.null(g1)) print(g1)
    
    cat("\nGráfico 2:\n")
    g2 <- grafico2(mix_data)
    if(!is.null(g2)) print(g2)
    
    cat("\nGráfico 4:\n")
    g4 <- grafico4(mix_data)
    if(!is.null(g4)) print(g4)
  }
  
  # 8. Modelos logísticos (CON MANEJO DE ERRORES)
  if(verbose) cat("\n8. Ajustando modelos logísticos...\n")
  
  models <- tryCatch({
    amr_logistic(
      geno_table = geno_table,
      pheno_table = pheno_table,
      antibiotic = antibiotic,
      drug_class_list = drug_class,
      geno_sample_col = geno_sample_col,
      pheno_sample_col = pheno_sample_col,
      maf = maf,
      sir_col = sir_col,
      marker_col = marker_col
    )
  }, error = function(e) {
    cat("ERROR en modelos logísticos:", e$message, "\n")
    cat("Continuando con el análisis...\n")
    return(NULL)  # Retornar NULL si hay error
  }, warning = function(w) {
    cat("ADVERTENCIA en modelos logísticos:", w$message, "\n")
    # Intentar obtener el resultado a pesar de la advertencia
    suppressWarnings(
      amr_logistic(
        geno_table = geno_table,
        pheno_table = pheno_table,
        antibiotic = antibiotic,
        drug_class_list = drug_class,
        geno_sample_col = geno_sample_col,
        pheno_sample_col = pheno_sample_col,
        maf = maf,
        sir_col = sir_col,
        marker_col = marker_col
      )
    )
  })
  
  if(!is.null(models)) {
    cat("Modelos logísticos ajustados correctamente\n")
  }
  
  # 9. Análisis solo PPV (CON MANEJO DE ERRORES)
  if(verbose) cat("\n9. Análisis PPV...\n")
  
  soloPPV <- tryCatch({
    solo_ppv_analysis(
      geno_table = geno_table,
      pheno_table = pheno_table,
      antibiotic = antibiotic,
      drug_class_list = drug_class,
      geno_sample_col = geno_sample_col,
      pheno_sample_col = pheno_sample_col,
      sir_col = sir_col,
      marker_col = marker_col
    )
  }, error = function(e) {
    cat("ERROR en análisis PPV:", e$message, "\n")
    cat("Continuando con el análisis...\n")
    return(NULL)  # Retornar NULL si hay error
  }, warning = function(w) {
    cat("ADVERTENCIA en análisis PPV:", w$message, "\n")
    # Intentar obtener el resultado a pesar de la advertencia
    suppressWarnings(
      solo_ppv_analysis(
        geno_table = geno_table,
        pheno_table = pheno_table,
        antibiotic = antibiotic,
        drug_class_list = drug_class,
        geno_sample_col = geno_sample_col,
        pheno_sample_col = pheno_sample_col,
        sir_col = sir_col,
        marker_col = marker_col
      )
    )
  })
  
  if(!is.null(soloPPV)) {
    cat("Análisis PPV completado correctamente\n")
  }
  
  # 10. Análisis UpSet (CON MANEJO DE ERRORES)
  if(verbose) cat("\n10. Generando gráfico UpSet...\n")
  
  mic_upset <- tryCatch({
    amr_upset(
      mix_data,
      min_set_size = min_set_size,
      order = "value",
      assay = "mic"
    )
  }, error = function(e) {
    cat("ERROR en gráfico UpSet:", e$message, "\n")
    cat("Continuando...\n")
    return(NULL)
  })
  
  if(verbose) cat("\nAnálisis completado\n\n")
  
  # Retornar todos los resultados en una lista
  resultados <- list(
    antibiotic = antibiotic,
    drug_class = drug_class,
    mix_data = mix_data,
    resistentes = res_resistentes,
    nwt = res_nwt,
    mic = res_mic,
    mutaciones_R = mut_r,
    mutaciones_NWT = mut_nwt,
    models = models,  # Puede ser NULL si hubo error
    ppv_analysis = soloPPV,  # Puede ser NULL si hubo error
    upset_plot = mic_upset  # Puede ser NULL si hubo error
  )
  
  return(invisible(resultados))
}
save(analizar_antibiotico,file= "Analisis")
