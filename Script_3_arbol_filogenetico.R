#SCRIPT 3
#Dr. Guillermo Berumen Varela
#Árbol filogenético
####

# Verificar si la paquetería BiocManager esta
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Instalar Biostrings a partir de BiocManager
BiocManager::install("Biostrings")

#Cargar librería
library(Biostrings)

#Establecer ruta de trabajo (es la misma, solo por si se esta este script de 0)
setwd("/home/biocata/Documentos/CONGRESO")

# Leer archivo FASTA
# Modificar la ruta
CDS_chirimoya <- readDNAStringSet("/home/biocata/Documentos/CONGRESO/anche102_CDS_annot.fasta")  

# genes_interes puede tener fragmentos como "ppo_gene"
genes_interes <- c("Anche102Chr1g0114080.1",
                   "Anche102Chr1g0113950.1",
                   "Anche102Chr7g0098410.1",
                   "Anche102Scf0007g0001170.1",
                   "Anche102Chr7g0098120.1")

# Filtrar secuencias cuyo nombre contenga alguno de los términos de interés
#sapply aplicar a--
#grepl es buscar un patrón que haga match 
#functiones crea una funcion
#apply
filtro <- sapply(genes_interes, function(x) grepl(x, names(CDS_chirimoya)))
filtro_comb <- apply(filtro, 1, any)

# Extraer secuencias
fasta_filtrado <- CDS_chirimoya[filtro_comb]

# Guardar archivo
writeXStringSet(fasta_filtrado, filepath = "secuencias_filtradas.fasta")


######### ARBOL FILOGENETICO ##################
BiocManager::install("msa")
BiocManager::install("ggtree")
library(msa)
library(phangorn)
library(ggtree)
library(ggplot2)
library(ape)
library(dplyr)

# Cargar secuencias
sequences <- readDNAStringSet("secuencias_filtradas.fasta")

#Alineamiento con MUSCLE
alignment <- msa(sequences, method = "Muscle")

# Convertir el alineamiento para phangorn
alignment_phangorn <- as.phyDat(alignment)

#p distance
dist_matrix <- dist.p(alignment_phangorn)

#Neighbor-Joining
tree_nj <- NJ(dist_matrix)

tree_ml <- pml(tree_nj, data = alignment_phangorn)
#tree_ml <- optim.pml(tree_ml, model = "JTT", optGamma = TRUE) # Optimización con máxima verosimilitud

# Exportar el árbol en formato Newick
write.tree(tree_ml$tree, file = "arbol.nwk")
tree <- read.tree("arbol.nwk")

# Graficar el árbol en formato rectangular con etiquetas
p <- ggtree(tree, layout = "rectangular") +
  geom_tiplab(size = 3, fontface = "italic") +
  theme_minimal()
print(p)
ggsave("arbol_filo.jpg", p, width = 10, height = 10, units="in", dpi = 300)
dev.off()
#-------------------------------------------------------------------------------