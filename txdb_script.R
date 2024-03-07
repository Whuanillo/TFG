library(GenomicFeatures)
gtf <- "gencode.vM34.annotation.gtf"
txdb.filename <- "gencode.vM34.annotation.sqlite"
txdb <- makeTxDbFromGFF(gtf)
saveDb(txdb, txdb.filename)
txdb <- loadDb(txdb.filename)
genes(txdb)
