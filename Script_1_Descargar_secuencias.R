#SCRIPT 1
#Dr. Guillermo Berumen Varela
#Taller Fundamentos de bioinformática
#SCRIPT para descargar secuencias usando rentrez de nucleótidos
#tomate 

#conocer la ruta de trabajo (getwd)
getwd()

# Establecer ruta de trabajo (setwd)
setwd("/home/biocata/Documentos/")

#Crear carpeta y cambiar directorio
dir.create("CONGRESO")

#Establecer la ruta de trabajo para todo el taller
setwd("/home/biocata/Documentos/CONGRESO")

#Verificar la ruta de trabajo
getwd()

#Verificar si se tiene la paquetería, sino esta instalada, se installa
if(!requireNamespace("rentrez", quietly = TRUE)) {
  install.packages("rentrez")
}

#Cargar paquetería (importante, siempre cargar paqueterías), accede a secuencias del NCBI, interfaz para usar entrez
library(rentrez)

# Definir término de búsqueda: PPO en mango (nucleótidos)
# Definir nombre del objeto (puede ser español o ingles) y posteriormente generar la búsqueda
# Nombre del gen Y nomber del organismo (INGLES)
search_term <- "polyphenol oxidase[Title] AND Mangifera indica[Organism]"

# Buscar IDs en la base de datos nucleotide
search_res <- entrez_search(db = "nucleotide", term = search_term, retmax = 10)

# Ver IDs encontrados
print(search_res$ids)

# Descargar secuencias en formato FASTA
fasta_seqs <- entrez_fetch(db = "nucleotide", id = search_res$ids, rettype = "fasta", retmode = "text")

# Guardar en archivo
write(fasta_seqs, file = "PPO_mango_gene.fasta")

# Leer las primeras líneas del archivo
file <- "PPO_mango_gene.fasta"
lineas <- readLines(file)

# Ver las primeras 10 líneas
head(lineas, 10)


