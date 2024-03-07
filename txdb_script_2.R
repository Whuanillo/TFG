library(GenomicFeatures)
txdb.filename <- "~/db/refanno/gencode.vM34.annotation.sqlite"
txdb <- loadDb(txdb.filename)
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, k, "GENEID", "TXNAME")
files <- file.path("salmon_counts" , "quant.sf")
names(files) <- "SRR18716583"
library(tximport)
txi.salmon.g <- tximport(files, type = "salmon", tx2gene = tx2gene)
txi.salmon.t <- tximport(files, type = "salmon", txOut = TRUE)

# Gene-level
head(txi.salmon.g$counts) #[,1:6] porque trabajaremos con 6 sets de datos

# Transcript-level
head(txi.salmon.t$counts) #[,1:6]

# Performs scaling for DTU analysis
txi.salmon.s <- tximport(files, type = "salmon", 
txOut = TRUE, tx2gene = tx2gene, countsFromAbundance = "dtuScaledTPM")

# Transcript-level, dtuScaledTPM
head(txi.salmon.s$counts) #[,1:6]

txi.salmon <- tximport(files, type = "none", txIn = TRUE, txOut = TRUE, 
txIdCol = "Name", abundanceCol = "TPM", 
countsCol = "NumReads", lengthCol = "Length",
importer = function(x) readr::read_tsv(x))
head(txi.salmon$counts) #[,1:6]
