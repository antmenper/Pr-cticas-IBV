Mic<-function(x){
  x%>% # sustituimos x por la combinación resultante de tabla fenotipica y genotipica
    count(mic)%>%
    mutate(freq_relativa = n / sum(n), # Calculo de la frecuencia relativa
           porcentaje = freq_relativa * 100, # Calculo del porcentaje 
           freq_acumulada = cumsum(n), # Frecuencia acumulada nos permite comprobar de forma rápida que se tienen en cuenta todos los datos
           porc_acumulado = cumsum(porcentaje) # Porcentaje acumulado.
    )
}

save(Mic,file="Mic_Descriptivo.RData")
