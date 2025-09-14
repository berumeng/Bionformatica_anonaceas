#SCRIPT 2 
#Dr. Guillermo Berumen Varela
#BLAST NUCLEÓTIDOS
####

# Descargar CDS
#system("curl -O https://ihsmsubtropicals.uma.es/downloads/Annona%20cherimola/Sequences/anche102_CDS_annot.fasta.gz")

# Descomprimir
#system("gzip -d anche102_CDS_annot.fasta.gz")

#Archivo descomprimido que contiene los CDS de chirimoya

# Crear base de datos BLAST 
#-in colocar el archivo fasta de donde se realizara la base de datos
#-dbtype tipo de base de datos (nucleótidos)
#-out nombre del archivo de salida (cualquiera)
system("makeblastdb -in anche102_CDS_annot.fasta -dbtype nucl -out chirimoyaCDS_db")

# Verificar si la base de datos de chirimoya existe (deberías tener .nin, .nsq, .nhr)
file.exists("chirimoyaCDS_db.nsq")

#Realizar un BLAST
# tblastx Para encontrar homólogos funcionales que no se detectan bien a nivel de ADN.
# mejor si los genes están conservados funcionalmente pero divergentes a nivel de secuencia de ADN.
#-query es el archivo que contiene las secuencias que queremos identificar homólogos
#-db es la base de datos que creamos
#-out nombre del archivo que queremos colocar
#-evalue valor de e (estadistico)
#query id = identificador de tu secuencia query
#subject id = identificador del hit base de datos
# % identity = porcentaje de identidad
# Alignment length
# mismatches = número de desajustes
# gap opens = número de apertura "huecos"
system("tblastx -query PPO_tomato_gene.fasta -db chirimoyaCDS_db -out resultados_tomato_tblastx_e-5.txt -evalue 1e-5 -outfmt 6")

# Leer resultados
results <- read.table("resultados_tomato_tblastx_e-5.txt", header = FALSE)
head(results)

# Asignar nombres de columna manualmente
colnames(results) <- c("query_id", "subject_id", "%identity", "alinm_length", "mismatch", 
                          "gapopen", "qstart", "qend", "sstart", "send", 
                          "evalue", "bitscore")

#install.packages("dplyr")
library(dplyr)
# Eliminar todas las filas duplicadas
df_sin_duplicados  <- results %>%
  filter(query_id != subject_id, `%identity` > 80)


# Extraer valores únicos de la columna subject_id
genes_homologos_chirimoya <- data.frame(subject_id = unique(df_sin_duplicados$subject_id))

# Guardar como archivo .txt o .tsv con formato de tabla
write.table(genes_homologos_chirimoya, file = "genes_homolo_chirimoya_POD_tomato.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")
