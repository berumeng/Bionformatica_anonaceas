#SCRIPT 4
#Dr. Guillermo Berumen Varela
#MAPA DE CALOR Y MAPA DE CORRELACION
####

#Instalar paqueterías
install.packages("pheatmap")  # Solo si no está instalada

#Cargar paqueterías
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
setwd("/home/biocata/Documentos/CONGRESO")

# A partir del archivo genes_homolog_chirimoya_PPO.txt
# Usar los ID y colocarlos en la plataforma https://ihsmsubtropicals.uma.es/easy_gdb/tools/expression/expression_input.php
#Anche102Chr1g0114080.1
#Anche102Chr1g0113950.1
#Anche102Chr7g0098410.1
#Anche102Chr6g0075250.1
#Anche102Chr3g0051290.1
#Anche102Chr6g0016200.1
#Anche102Chr6g0016350.1
#Anche102Chr6g0015750.1
#Anche102Chr6g0015850.1
#Anche102Scf0007g0001170.1
#Anche102Chr6g0016150.1
#Anche102Chr6g0016160.1
#Anche102Chr7g0098120.1

#Descargar el archivo con datos de expresión.

# Leer archivo CSV con filas = genes y columnas = estados
expresion <- read.csv("IHSM Subtropicals_POD_tom.csv", row.names = 1)

# Asegúrate de que las filas sean los genes y columnas los estados
head(expresion)

#Seleccionar solo las columnas que son de interés (solo las primeras 10)
library(dplyr)
expresion <- expresion %>% 
  select(1:10) %>%
filter(rowSums(select(., 1:10) != 0) > 0)

#mapa de calor expresion por gen 
Mapa_calor <- pheatmap(expresion,
         scale = "row",               # escalar por gen (opcional pero común)
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         color = colorRampPalette(c("royalblue4", "white", "firebrick"))(100),
         fontsize_row = 10,
         fontsize_col = 8,
         border_color = "black",
         main = "Expresión de genes que codifican a PPO en chirimoya en distintas etapas de maduración")
print(Mapa_calor)
ggsave("mapa de calor_POD.jpg", Mapa_calor, width = 10, height = 10, units="in", dpi = 300)
dev.off()

#Calcular la matriz de correlación entre genes
#Se usa t() para la correlación entre filas 
matriz_cor_genes <- cor(t(expresion), method = "pearson")

# Generar mapa de calor de correlación entre genes
Correlacion <- pheatmap(matriz_cor_genes,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         color = colorRampPalette(c("royalblue4", "white", "firebrick"))(100),
         main = "Mapa de calor de correlación entre genes",
         fontsize_row = 6,
         fontsize_col = 6)
print(Correlacion)
ggsave("correlacion_POD.jpg", Correlacion, width = 10, height = 10, units="in", dpi = 300)
dev.off()
