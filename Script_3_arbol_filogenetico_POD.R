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

# genes_interes puede tener fragmentos como "POD_gene"
genes_interes <- c("Anche102Chr1g0115690.1",
                   "Anche102Chr5g0010720.1",
                   "Anche102Chr1g0018570.1",
                   "Anche102Chr1g0019060.1",
                   "Anche102Chr4g0023540.1",
                   "Anche102Chr6g0028230.1",
                   "Anche102Chr3g0056990.1",
                   "Anche102Chr3g0057110.1",
                   "Anche102Chr1g0027890.1",
                   "Anche102Chr3g0054110.1",
                   "Anche102Chr1g0019070.1",
                   "Anche102Chr1g0018580.1",
                   "Anche102Chr5g0013880.1",
                   "Anche102Chr5g0010750.1",
                   "Anche102Chr3g0036090.1",
                   "Anche102Chr5g0010710.1",
                   "Anche102Chr5g0013900.1",
                   "Anche102Chr4g0031590.1",
                   "Anche102Chr5g0047250.1",
                   "Anche102Chr5g0014210.1",
                   "Anche102Chr3g0036400.1",
                   "Anche102Chr1g0027830.1",
                   "Anche102Scf0158g0000010.1",
                   "Anche102Chr1g0027940.1",
                   "Anche102Chr1g0027750.1",
                   "Anche102Chr6g0089010.1",
                   "Anche102Chr7g0026550.1",
                   "Anche102Chr1g0061790.1",
                   "Anche102Chr6g0067920.1",
                   "Anche102Chr5g0054680.1",
                   "Anche102Chr1g0027790.2",
                   "Anche102Chr1g0027790.1",
                   "Anche102Chr6g0070870.1",
                   "Anche102Chr2g0035340.1",
                   "Anche102Chr2g0035480.1",
                   "Anche102Chr1g0014450.1",
                   "Anche102Chr1g0014340.1",
                   "Anche102Chr2g0035470.1",
                   "Anche102Chr7g0037640.1",
                   "Anche102Chr7g0037640.2",
                   "Anche102Chr7g0037930.1")

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
writeXStringSet(fasta_filtrado, filepath = "secuencias_filtradas_POD.fasta")


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
sequences <- readDNAStringSet("secuencias_filtradas_POD.fasta")

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
write.tree(tree_ml$tree, file = "arbol_POD.nwk")
tree <- read.tree("arbol_POD.nwk")

# Graficar el árbol en formato rectangular con etiquetas
p <- ggtree(tree, layout = "rectangular") +
  geom_tiplab(size = 3, fontface = "italic") +
  theme_minimal()
print(p)
ggsave("arbol_filo_POD.jpg", p, width = 10, height = 10, units="in", dpi = 300)
dev.off()
#-------------------------------------------------------------------------------