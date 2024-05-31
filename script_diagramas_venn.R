#!/usr/bin/env Rscript
counts <- read.csv("CountsDE.csv" , row.names = 1)
counts <- round(counts , 0) #Para quitar los decimales
rownames(counts) <- gsub("\\..*" , "" , rownames(counts)) #Para quitar los decimales de los identificadores de los genes
info <- read.csv("sample_treat.csv" , row.names = 1)

Exp_counts <- read.csv("GSE200621_counts.csv" , row.names = 1)
NS_counts <- Exp_counts[4:6]
nt_counts <- Exp_counts[1:3]
Experimental_counts <- cbind(NS_counts,nt_counts)
Experimental_counts <- round(Experimental_counts , 0) #Para quitar los decimales
rownames(Experimental_counts) <- gsub("\\..*" , "" , rownames(Experimental_counts)) #Para quitar los decimales de los identificadores de los genes

library(ExpHunterSuite)

degh_out_one_pack <- main_degenes_Hunter(raw=counts,
                                         target=info,
                                         modules="D")


degh_out_one_pack_2 <- main_degenes_Hunter(raw=Experimental_counts,
                                         target=info,
                                         modules="D")

genes <- (degh_out_one_pack$DE_all_genes)
Prev_genes <- subset(genes, genes_tag == "PREVALENT_DEG")

genes_exp <- (degh_out_one_pack_2$DE_all_genes)
Exp_genes <- subset(genes_exp, genes_tag == "PREVALENT_DEG")

library(biomaRt) #interfaz para el acceso a las bases de datos de BioMart

#####Primera parte de obtención de nombres de nustros genes

Row_names <- rownames(Prev_genes)

# Conectar a la base de datos de Ensembl para el ratón (Mus musculus)
ensembl_mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl") #useMart() nos permite conectar con la base de datos de Ensembl para el ratón

# Obtener información sobre los genes
gene_info <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "description"), #el comando getBM()  obtener información sobre los genes
                   filters = "ensembl_gene_id",
                   values = Row_names,
                   mart = ensembl_mouse)

#función merge() para combinar los datos obtenidos de getBM() con Row_names, asegurándonos así de que los identificadores de Ensembl se mantengan incluso para los genes sin información correspondiente en la base de datos.
merged_info <- merge(data.frame(ensembl_gene_id = Row_names), gene_info, by = "ensembl_gene_id", all.x = TRUE)

nombres_genes <- merged_info$external_gene_name

#####Segunda parte de obtención de nombres de los genes del artículo

Row_names <- rownames(Exp_genes)

# Conectar a la base de datos de Ensembl para el ratón (Mus musculus)
ensembl_mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl") #useMart() nos permite conectar con la base de datos de Ensembl para el ratón

# Obtener información sobre los genes
gene_info <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "description"), #el comando getBM()  obtener información sobre los genes
                   filters = "ensembl_gene_id",
                   values = Row_names,
                   mart = ensembl_mouse)

#función merge() para combinar los datos obtenidos de getBM() con Row_names, asegurándonos así de que los identificadores de Ensembl se mantengan incluso para los genes sin información correspondiente en la base de datos.
merged_info <- merge(data.frame(ensembl_gene_id = Row_names), gene_info, by = "ensembl_gene_id", all.x = TRUE)

nombres_genes_Exp <- merged_info$external_gene_name

#DIAGRAMAS DE VENN

library(VennDiagram)

nombres_genes <- setdiff(nombres_genes, "")
nombres_genes_Exp <- setdiff(nombres_genes_Exp, "")
nombres_genes_Exp <- setdiff(nombres_genes_Exp, NA)

#Ahora obtenemos un diagrama de Venn para los genes DE obtenidos y los genes DE del artículo

venn.diagram(
        x = list(nombres_genes, nombres_genes_Exp),
        category.names = c("Genes DE" , "Genes DE (Artículo)"),
        filename = 'venn_diagramm_genesDE_vs_publicados.png',
        output=TRUE,
        col = c("black", "black"),
        fill = c("red", "green"),
        cex = 2,
        cat.cex = 2,
        lwd = 2,
        imagetype = "png" ,
        height = 3840 , 
        width = 3840 ,
        resolution = 300 ,
        margin = 0.05
)

#Ahora extraemos los nombres de HPO y los pasamos a los genes de ratón para poder compararlos

nombres_hpo <- read.table("genes_for_HP_0000510.txt" , row.names = 1)
hpo_nombres <- nombres_hpo[,1]

# Cargar la librería biomaRt
library(biomaRt)

# Configurar el mart de Ensembl para humanos
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Buscar el ID de Ensembl para un gen HGNC
hgnc_name <- hpo_nombres
filters <- "hgnc_symbol"
attributes <- c("ensembl_gene_id", "hgnc_symbol", "chromosome_name")
result <- getBM(attributes = attributes, filters = filters, values = hgnc_name, mart = ensembl)

# Verificar si se encontró el ID de Ensembl
if (nrow(result) == 0) {
  stop("No Ensembl ID found for HGNC name")
}
ensembl_gene_id <- result$ensembl_gene_id

# Cambiar a un mart de homología para obtener homólogos de ratón
homology_mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "www.ensembl.org")

# Listar los atributos disponibles para la homología
listAttributes(homology_mart)[grepl("mmusculus_homolog", listAttributes(homology_mart)$name),]

# Obtener parálogos en ratón
attributes <- c("ensembl_gene_id", "mmusculus_homolog_ensembl_gene", "mmusculus_homolog_associated_gene_name")
filters <- "ensembl_gene_id"
values <- ensembl_gene_id
mouse_ortologs <- getBM(attributes = attributes, filters = filters, values = values, mart = homology_mart)

ortologos_hpo <- mouse_ortologs$mmusculus_homolog_associated_gene_name

ortologos_hpo <- setdiff(ortologos_hpo, "")

#Ahora obtenemos un diagrama de Venn para los genes DE obtenidos y los genes de HPO

venn.diagram(
        x = list(nombres_genes, ortologos_hpo),
        category.names = c("Genes DE" , "Genes de HPO"),
        filename = 'venn_diagramm_DE_vs_hpo.png',
        output=TRUE,
        col = c("black", "black"),
        fill = c("red", "blue"),
        cex = 2,
        cat.cex = 2,
        lwd = 2,
        imagetype = "png" ,
        height = 3840 , 
        width = 3840 ,
        resolution = 300 ,
        margin = 0.05
)

#Ahora obtenemos un diagrama de Venn para comparación de los tres casos vistos

venn.diagram(
        x = list(nombres_genes, nombres_genes_Exp, ortologos_hpo),
        category.names = c("Genes DE" , "Genes DE (Artículo)" , "Genes de HPO"),
        filename = 'venn_diagramm_completo.png',
        output=TRUE,
        col = c("black", "black", "black"),
        fill = c("red", "green", "blue"),
        cex = 2,
        cat.cex = 2,
        lwd = 2,
        imagetype = "png" ,
        height = 3840 , 
        width = 3840 ,
        resolution = 300 ,
        margin = 0.05
)

##Nivel de significacia de las intersecciones

phyper(1328+37,501+2+1328+37,15382-501-2-1328-37,1328+37+2+152,lower.tail=F) #Para la interseccion de DEG y DEG del articulo

phyper(2+37,501+2+1328+37,15382-501-2-1328-37,2+37+2+194,lower.tail=F) #Para la interseccion de DEG y los genes de HPO

phyper(2+37,1328+37+2+152,15382-(1328+37+2+152),2+37+2+194,lower.tail=F) #Para la interseccion de DEG y los genes de HPO

##Creación de tabla de genes intersecantes del diagrama de Venn 

# Convertir las cadenas en vectores de caracteres
genes1 <- strsplit(nombres_genes, ",")
genes2 <- strsplit(nombres_genes_Exp, ",")
genes3 <- strsplit(ortologos_hpo, ",")

# Encontrar la intersección de los tres vectores
interseccion_genes <- Reduce(intersect, list(genes1, genes2, genes3))
interseccion_genes <- as.character(interseccion_genes)

genes_interseccion_info <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "description"), #el comando getBM()  obtener información sobre los genes
                   filters = "external_gene_name",
                   values = interseccion_genes,
                   mart = ensembl_mouse)

genes_intersección_info

write.table(genes_interseccion_info, file = "genes_interseccion_info.txt")
