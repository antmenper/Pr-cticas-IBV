analizar_logistf_automatico <- function(antibiotic, 
                                        drug_class,
                                        geno_table = geno2,
                                        pheno_table = y_largo,
                                        geno_sample_col = "Name",
                                        pheno_sample_col = "biosample",
                                        sir_col = "pheno",
                                        min_count = 1,
                                        verbose = TRUE,
                                        mostrar_summary = TRUE) {
  
  if(verbose) {
    cat("    ANÁLISIS LOGISTF AUTOMÁTICO               \n")
    
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
    sir_col = sir_col
  )
  
  if(verbose) cat("   Variables iniciales:", ncol(mix_data), "\n")
  
  # 2. Preparar datos: eliminar columnas no necesarias
  if(verbose) cat("\n2. Preparando datos...\n")
  
  columnas_eliminar <- c("id", "pheno", "mic", "NWT")
  data_clean <- mix_data[, setdiff(names(mix_data), columnas_eliminar)]
  
  if(verbose) cat("Columnas eliminadas (id, pheno, mic, NWT)\n")
  cat("Variables restantes:", ncol(data_clean), "\n")
  
  # 3. Filtrar por conteo mínimo
  if(verbose) cat("\n3. Filtrando por conteo mínimo (>", min_count, ")...\n")
  
  conteos <- colSums(data_clean, na.rm = TRUE)
  data_filtered <- data_clean[, conteos > min_count]
  
  if(verbose) {
    cat("Variables eliminadas por bajo conteo:", 
        ncol(data_clean) - ncol(data_filtered), "\n")
    cat("Variables restantes:", ncol(data_filtered), "\n")
  }
  
  # 4. Ajustar GLM para identificar combinaciones lineales exactas
  if(verbose) cat("\n4. Ajustando GLM para detectar combinaciones lineales...\n")
  
  glm_model <- tryCatch({
    glm(R ~ ., data = data_filtered, family = binomial())
  }, error = function(e) {
    return(NULL)
  }, warning = function(w) {
    # Suprimir warnings porque esperamos algunos
    suppressWarnings(glm(R ~ ., data = data_filtered, family = binomial()))
  })
  
  if(is.null(glm_model)) {
    return(NULL)
  }
  
  if(verbose) cat("\n5. Identificando variables con coeficientes NA (combinaciones lineales)...\n")
  
  # Obtener todos los coeficientes del modelo
  coef_modelo <- coef(glm_model)
  
  # Eliminar el intercepto
  coef_modelo <- coef_modelo[names(coef_modelo) != "(Intercept)"]
  
  # Identificar cuáles coeficientes son NA
  vars_con_coef_na <- names(coef_modelo)[is.na(coef_modelo)]
  
  # Limpiar nombres (quitar backticks si existen)
  vars_na <- gsub("`", "", vars_con_coef_na)
  
  if(length(vars_na) > 0) {
    if(verbose) {
      cat("Variables con coeficientes NA detectadas:", length(vars_na), "\n")
      cat("   Variables a eliminar:\n")
      for(v in vars_na) {
        cat("      •", v, "\n")
      }
    }
    
    # Eliminar variables con coeficientes NA
    data_final <- data_filtered[, !names(data_filtered) %in% vars_na]
    
    if(verbose) {
      cat("\n Variables eliminadas:", length(vars_na), "\n")
      cat("Variables finales para logistf:", ncol(data_final) - 1, "\n")
    }
  } else {
    if(verbose) cat("No se detectaron coeficientes NA\n")
    data_final <- data_filtered
  }
  # 6. Mostrar summary del GLM si se solicita
  if(mostrar_summary) {
    cat("\n", rep("=", 50), "\n", sep="")
    cat("SUMMARY GLM (antes de eliminar NAs)\n")
    cat(rep("=", 50), "\n", sep="")
    print(summary(glm_model))
    cat("\n")
  }
  
  # 7. Verificar y limpiar variables problemáticas adicionales
  if(verbose) cat("\n6. Verificaciones adicionales antes de logistf...\n")
  
  # Eliminar variables con varianza cero
  vars_no_r <- setdiff(names(data_final), "R")
  variances <- sapply(data_final[, vars_no_r], var, na.rm = TRUE)
  vars_zero_var <- names(variances)[variances == 0 | is.na(variances)]
  
  if(length(vars_zero_var) > 0) {
    if(verbose) {
      cat("Eliminando variables con varianza cero:", length(vars_zero_var), "\n")
      for(v in vars_zero_var) cat("      •", v, "\n")
    }
    data_final <- data_final[, !names(data_final) %in% vars_zero_var]
  }
  
  # Verificar multicolinealidad perfecta con correlaciones
  if(ncol(data_final) > 2) {  # Más de R + 1 variable
    vars_no_r <- setdiff(names(data_final), "R")
    if(length(vars_no_r) > 1) {
      cor_matrix <- cor(data_final[, vars_no_r], use = "complete.obs")
      
      # Buscar correlaciones perfectas (1 o -1, excluyendo diagonal)
      diag(cor_matrix) <- 0
      perfect_cor <- which(abs(cor_matrix) > 0.999, arr.ind = TRUE)
      
      if(nrow(perfect_cor) > 0) {
        # Eliminar una de cada par correlacionado
        vars_eliminar_cor <- unique(rownames(perfect_cor))
        if(verbose) {
          cat("Eliminando variables con correlación perfecta:", 
              length(vars_eliminar_cor), "\n")
          for(v in vars_eliminar_cor) cat("      •", v, "\n")
        }
        data_final <- data_final[, !names(data_final) %in% vars_eliminar_cor]
      }
    }
  }
  
  if(verbose) {
    cat("Variables finales después de limpieza:", ncol(data_final) - 1, "\n")
  }
  
  # 8. Ajustar modelo logistf con datos limpios
  if(verbose) cat("\n7. Ajustando modelo logistf...\n")
  
  model_logistf <- tryCatch({
    logistf::logistf(R ~ ., 
                     data = data_final, 
                     pl = FALSE, 
                     firth = TRUE)
  }, error = function(e) {
    cat("ERROR al ajustar logistf:", e$message, "\n")
    cat("Intentando sin penalización de Firth...\n")
    
    tryCatch({
      logistf::logistf(R ~ ., 
                       data = data_final, 
                       pl = FALSE, 
                       firth = FALSE)
    }, error = function(e2) {
      cat("También falló sin Firth. Intentando GLM estándar...\n")
      tryCatch({
        glm(R ~ ., data = data_final, family = binomial())
      }, error = function(e3) {
        cat("Todos los métodos fallaron\n")
        return(NULL)
      })
    })
  })
  
  if(!is.null(model_logistf)) {
    if(verbose) cat("Modelo ajustado exitosamente\n")
    
    # Mostrar summary
    if(mostrar_summary) {
      cat("\n", rep("=", 50), "\n", sep="")
      cat("SUMMARY MODELO FINAL\n")
      cat(rep("=", 50), "\n", sep="")
      print(summary(model_logistf))
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
    variables_eliminadas_conteo = setdiff(names(data_clean), names(data_filtered)),
    variables_eliminadas_na = vars_na,
    variables_eliminadas_var_cero = if(exists("vars_zero_var")) vars_zero_var else character(0),
    variables_eliminadas_correlacion = if(exists("vars_eliminar_cor")) vars_eliminar_cor else character(0),
    glm_model = glm_model,
    logistf_model = model_logistf,
    n_variables_inicial = ncol(mix_data),
    n_variables_final = ncol(data_final) - 1
  )
  
  return(invisible(resultados))
}
save(analizar_logistf_automatico,file= "Análisis numérico automatizado")
