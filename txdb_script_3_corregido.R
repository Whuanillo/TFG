library(GenomicFeatures)
txdb.filename <- "~/db/refanno/gencode.vM34.annotation.sqlite"
txdb <- loadDb(txdb.filename)
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, k, "GENEID", "TXNAME")
files <- file.path("salmon_data", list.files("salmon_data") , "quant.sf")
names(files) <- list.files("salmon_data")
library(tximport)
txi.salmon.g <- tximport(files, type = "salmon", tx2gene = tx2gene)
txi.salmon.t <- tximport(files, type = "salmon", txOut = TRUE)

# Gene-level
head(txi.salmon.g$counts)[,1:6] 

# Transcript-level
head(txi.salmon.t$counts)[,1:6]

# Performs scaling for DTU analysis
txi.salmon.s <- tximport(files, type = "salmon", 
txOut = TRUE, tx2gene = tx2gene, countsFromAbundance = "dtuScaledTPM")

# Transcript-level, dtuScaledTPM
head(txi.salmon.s$counts)[,1:6]

txi.salmon <- tximport(files, type = "none", txIn = TRUE, txOut = TRUE, 
txIdCol = "Name", abundanceCol = "TPM", 
countsCol = "NumReads", lengthCol = "Length",
importer = function(x) readr::read_tsv(x))
head(txi.salmon$counts)[,1:6]