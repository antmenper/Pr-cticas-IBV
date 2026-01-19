analizar_lm_mic_automatico <- function(antibiotic, 
                                       drug_class,
                                       geno_table = geno2,
                                       pheno_table = y_largo,
                                       geno_sample_col = "Name",
                                       pheno_sample_col = "biosample",
                                       sir_col = "pheno",
                                       min_count = 5,
                                       verbose = TRUE,
                                       mostrar_summary = TRUE) {
  
  if(verbose) {
    cat("\n╔════════════════════════════════════════════════╗\n")
    cat("║     ANÁLISIS LM MIC AUTOMÁTICO                 ║\n")
    cat("╚════════════════════════════════════════════════╝\n\n")
    cat("Antibiótico:", antibiotic, "\n")
    cat("Clase:", drug_class, "\n")
    cat(rep("-", 50), "\n\n", sep="")
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
  
  if(verbose) cat("   Variables iniciales:", ncol(mix_data), "\n")
  
  # 2. Preparar datos: eliminar columnas no necesarias
  if(verbose) cat("\n2. Preparando datos...\n")
  
  columnas_eliminar <- c("id", "pheno", "NWT", "R")
  data_mic <- mix_data[, setdiff(names(mix_data), c("id","pheno","NWT","R"))]
  
  if(verbose) cat("   Columnas eliminadas (id, pheno, NWT, R)\n")
  
  # 3. Transformar MIC
  if(verbose) cat("\n3. Transformando variable MIC...\n")
  
  data_mic$mic <- as.numeric(data_mic$mic)
  data_mic <- as.data.frame(data_mic)
  data_mic$mic <- log2(data_mic$mic)
  
  if(verbose) cat("   ✓ MIC convertida a log2(MIC)\n")
  
  # 4. Filtrar por conteo mínimo
  if(verbose) cat("\n4. Filtrando por conteo mínimo (>", min_count, ")...\n")
  
  conteos <- colSums(data_mic, na.rm = TRUE)
  data_filtered <- data_mic[, conteos > min_count]
  
  if(verbose) {
    cat("   Variables eliminadas por bajo conteo:", 
        ncol(data_mic) - ncol(data_filtered), "\n")
    cat("   Variables restantes:", ncol(data_filtered), "\n")
  }
  
  # 5. Ajustar modelo lineal
  if(verbose) cat("\n5. Ajustando modelo lineal...\n")
  
  modelo_lm <- tryCatch({
    lm(mic ~ ., data = data_filtered)
  }, error = function(e) {
    cat("   ✗ ERROR al ajustar modelo:", e$message, "\n")
    return(NULL)
  })
  
  if(is.null(modelo_lm)) {
    cat("   ✗ No se pudo ajustar el modelo\n")
    return(NULL)
  }
  
  # 6. Identificar variables con coeficientes NA
  if(verbose) cat("\n6. Identificando variables con coeficientes NA...\n")
  
  # Obtener todos los coeficientes del modelo
  coef_modelo <- coef(modelo_lm)
  
  # Eliminar el intercepto
  coef_modelo <- coef_modelo[names(coef_modelo) != "(Intercept)"]
  
  # Identificar cuáles coeficientes son NA
  vars_con_coef_na <- names(coef_modelo)[is.na(coef_modelo)]
  
  # Limpiar nombres (quitar backticks si existen)
  vars_na <- gsub("`", "", vars_con_coef_na)
  
  if(length(vars_na) > 0) {
    if(verbose) {
      cat("   ⚠️  Variables con coeficientes NA detectadas:", length(vars_na), "\n")
      cat("   Variables a eliminar:\n")
      for(v in vars_na) {
        cat("      •", v, "\n")
      }
    }
    
    # Eliminar variables con coeficientes NA
    data_final <- data_filtered[, !names(data_filtered) %in% vars_na]
    
    if(verbose) {
      cat("\n   ✓ Variables eliminadas:", length(vars_na), "\n")
      cat("   ✓ Variables finales:", ncol(data_final) - 1, "\n")
    }
  } else {
    if(verbose) cat("   ✓ No se detectaron coeficientes NA\n")
    data_final <- data_filtered
  }
  
  # 7. Verificaciones adicionales
  if(verbose) cat("\n7. Verificaciones adicionales...\n")
  
  # Eliminar variables con varianza cero
  vars_no_mic <- setdiff(names(data_final), "mic")
  variances <- sapply(data_final[, vars_no_mic], var, na.rm = TRUE)
  vars_zero_var <- names(variances)[variances == 0 | is.na(variances)]
  
  if(length(vars_zero_var) > 0) {
    if(verbose) {
      cat("   ⚠️  Eliminando variables con varianza cero:", length(vars_zero_var), "\n")
      for(v in vars_zero_var) cat("      •", v, "\n")
    }
    data_final <- data_final[, !names(data_final) %in% vars_zero_var]
  }
  
  if(verbose) {
    cat("   ✓ Variables finales después de limpieza:", ncol(data_final) - 1, "\n")
  }
  
  # 8. Ajustar modelo final con datos limpios
  if(verbose) cat("\n8. Ajustando modelo lineal final...\n")
  
  modelo_final <- tryCatch({
    lm(mic ~ ., data = data_final)
  }, error = function(e) {
    cat("   ✗ ERROR al ajustar modelo final:", e$message, "\n")
    return(NULL)
  })
  
  if(!is.null(modelo_final)) {
    if(verbose) cat("   ✓ Modelo ajustado exitosamente\n")
    
    # Mostrar summary
    if(mostrar_summary) {
      cat("\n", rep("=", 50), "\n", sep="")
      cat("SUMMARY MODELO FINAL\n")
      cat(rep("=", 50), "\n", sep="")
      print(summary(modelo_final))
      cat("\n")
      
      # Información adicional
      cat("R²:", round(summary(modelo_final)$r.squared, 4), "\n")
      cat("R² ajustado:", round(summary(modelo_final)$adj.r.squared, 4), "\n")
      cat("p-valor F:", format.pval(summary(modelo_final)$fstatistic), "\n")
    }
  }
  
  if(verbose) cat("\n¡Análisis completado!\n\n")
  
  # 9. Retornar resultados
  resultados <- list(
    antibiotic = antibiotic,
    drug_class = drug_class,
    mix_data = mix_data,
    data_filtered = data_filtered,
    data_final = data_final,
    variables_eliminadas_conteo = setdiff(names(data_mic), names(data_filtered)),
    variables_eliminadas_na = vars_na,
    variables_eliminadas_var_cero = if(exists("vars_zero_var")) vars_zero_var else character(0),
    modelo_lm = modelo_final,
    n_variables_inicial = ncol(mix_data),
    n_variables_final = ncol(data_final) - 1,
    r_squared = if(!is.null(modelo_final)) summary(modelo_final)$r.squared else NA,
    adj_r_squared = if(!is.null(modelo_final)) summary(modelo_final)$adj.r.squared else NA
  )
  
  return(invisible(resultados))
}
save(analizar_lm_mic_automatico,file = "Calculo mic lm.RData")
