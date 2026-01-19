# Prácticas-IBV

El archivo Analisis.R es un módulo contenedor que agrupa todas las funciones responsables del análisis descriptivo y la generación de gráficos mediante las herramientas ggplot2 y AMRgen.
El objetivo de este módulo es automatizar completamente el proceso de análisis, dada la gran cantidad de datos que se manejan en el proyecto. Su función principal es analizar_antibiotico(), que orquesta la ejecución del análisis completo.

El archivo Análisis numérico automatizado.R es un módulo que ejecuta el análisis mediante un modelo lineal generalizado (GLM) con distribución binomial, aplicando la corrección de sesgo de Firth para evitar problemas de separación completa que surgen frecuentemente con este tipo de datos. Su función principal es analizar_logistf_automatico(), que ejecuta el análisis completo de forma automatizada.
