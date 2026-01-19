mutaciones_R<-function(x){
  x%>%
    pivot_longer(cols= -c(id, pheno, mic, R, NWT), names_to="variable", values_to= "valor")%>% 
    count(R,variable,valor)%>%
    group_by(R, variable) %>%     
    mutate(
      freq_relativa = n / sum(n),
      porcentaje = freq_relativa * 100,
      freq_acumulada = cumsum(n),
      porc_acumulado = cumsum(porcentaje))
}

mutaciones_NWT<-function(x){
  x%>%
    pivot_longer(cols= -c(id, pheno, mic, R, NWT), names_to="variable", values_to= "valor")%>% 
    count(NWT,variable,valor)%>%
    group_by(NWT, variable) %>%     
    mutate(
      freq_relativa = n / sum(n),
      porcentaje = freq_relativa * 100,
      freq_acumulada = cumsum(n),
      porc_acumulado = cumsum(porcentaje))
}

save(mutaciones_R,mutaciones_NWT,file="Mutaciones descriptivos.RData")
