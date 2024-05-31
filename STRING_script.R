#!/usr/bin/env Rscript

counts <- read.csv("CountsDE.csv" , row.names = 1)
counts <- round(counts , 0) #Para quitar los decimales
rownames(counts) <- gsub("\\..*" , "" , rownames(counts)) #Para quitar los decimales de los identificadores de los genes
info <- read.csv("sample_treat.csv" , row.names = 1)

library(ExpHunterSuite)

degh_out_one_pack <- main_degenes_Hunter(raw=counts,
                                         target=info,
                                         modules="D")
#Ahora para usar STRING, vamos a realizar una serie de pasos para preparar los datos para introducirlos en el algoritmo
genes <- (degh_out_one_pack$DE_all_genes)

#Vamos a necesitar "traducir" los identificadores Ensembl para que STRING pueda leerlos

library(biomaRt) #interfaz para el acceso a las bases de datos de BioMart

Row_names <- rownames(counts)

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

#Ahora ya podemos usar los datos en STRING

library(STRINGdb)

string_db <- STRINGdb$new(  species=10090,
                           score_threshold=200, network_type="full", input_directory="")

vector_of_pvalues <- genes[,3]
vector_of_logFC <- genes[,1]
vector_of_gene_names <- nombres_genes

test_df <- data.frame(pvalue= vector_of_pvalues, logFC= vector_of_logFC, "gene" = vector_of_gene_names)

example1_mapped <- string_db$map(test_df, "gene", removeUnmappedRows = TRUE )
hits <- example1_mapped$STRING_id[1:200]
string_db$plot_network( hits )  

#Esto nos generará un mapa de genes con sus interacciones

#Lo siguiente es obtener los 4 clusteres mas probables
clustersList <- string_db$get_clusters(example1_mapped$STRING_id[1:200])
par(mfrow=c(2,2))
for(i in seq(1:4)){ string_db$plot_network(clustersList[[i]]) }
