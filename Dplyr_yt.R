#Instalar paqueterías (en caso no la tengan instalada)
install.packages("dplyr")

#cargar librería (paquetería)
library(dplyr)
# Crear un dataframe de ejemplo
datos <- tibble(
  gen = c("GeneA", "GeneB", "GeneC", "GeneD", "GenF"),
  expresion = c(10.2, 5.3, 20.1, 0.5, 2.0),
  categoria = c("metabolismo", "señalización", "estructura", "metabolismo", "metabolismo")
)


#filtrar marco de datos por filas con expresión mayor a 5
filter(datos, expresion > 5)

#seleccionar columnas
select(datos, gen, expresion)
select(datos, expresion, categoria)

#crear nuevas columnas
mutate(datos, log2_exp = log2(expresion))

#ordenar datos
arrange(datos, desc(expresion))

#agrupar datos y resumir los datos
#uso del operador pipe %>% (tubería)
datos %>%
  group_by(categoria) %>%
  summarise(promedio = mean(expresion))

#operador PIPE
#permite encadenar funciones y escribir código más limpio y legible
datos %>%
  filter(expresion > 10) %>%
  arrange(desc(expresion))

#creacion de marco de datos
expresion_genes <- tibble(
  gen = c("RGL1", "RGL2", "RGL3", "RGL4"),
  TPM = c(25.4, 12.0, 8.1, 0.2),
  condición = c("control", "tratado", "control", "tratado")
)

expresion_final <- expresion_genes %>%
  filter(TPM > 5) %>%
  mutate(log2TPM = log2(TPM)) %>%
  arrange(desc(log2TPM))
