NWT<-function(x){
  x%>%
    count(NWT)%>%
    mutate(freq_relativa = n / sum(n),
           porcentaje = freq_relativa * 100,
           freq_acumulada = cumsum(n),
           porc_acumulado = cumsum(porcentaje)
    )
}

grafico_NWT<-function(x){
  ggplot(x, aes(x = NWT)) +
    geom_bar(aes(y = after_stat(count / sum(count))))
}

  save(NWT,grafico_NWT,file="NWT_Descriptivos.RData")