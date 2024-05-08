counts <- read.csv("CountsDE.csv" , row.names = 1)
counts <- round(counts , 0) #Para quitar los decimales
rownames(counts) <- gsub("\\..*" , "" , rownames(counts)) #Para quitar los decimales de los identificadores de los genes
info <- read.csv("sample_treat.csv" , row.names = 1)

library(ExpHunterSuite)

degh_out_one_pack <- main_degenes_Hunter(raw=counts,
                                         target=info,
                                         modules="D")

fh_out_one_pack <- main_functional_hunter( #Perform enrichment analysis
        degh_out_one_pack,
        model_organism = 'Mouse', # Use specified organism database 
        enrich_dbs = c("MF", "BP", "CC", "Kegg", "Reactome"), # Enrichment analysis for GO, KEGG and Reactome
        enrich_methods = "o" # Use overepresentation analysis only
)

write_expression_report(exp_results=degh_out_one_pack) #Reporte html de expresión
write_enrich_files(func_results=fh_out_one_pack)
write_functional_report(hunter_results=degh_out_one_pack, 
                        func_results=fh_out_one_pack) #Reporte html de enriquecimiento main_functional_hunter

