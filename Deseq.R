library(tximport)
library(readr)
library(DESeq2)

######To Generate the Transcripts to geneids tx2gene.gencode file #########

#library(tximportData)
#library(GenomicFeatures)
#txdb <- makeTxDbFromGFF(file="gencode.v28.annotation.gff3.gz")
#txdb
#saveDb(x=txdb, file = "gencode.v28.annotation.TxDb")

#k <- keys(txdb, keytype = "TXNAME")
#tx2gene <- select(txdb, k, "GENEID", "TXNAME")
#head(k)
#head(tx2gene)
#dim(tx2gene)
#length(k)
#write.table(tx2gene, "tx2gene.gencode.v28.csv", sep = ",", row.names = FALSE)

#counts file generated from Salmon
dirs = "./salmon_out"
sys_dir <- system.file("extdata", package = "tximportData")
samples = dir(file.path(dirs))
files <- file.path(dirs, samples, "quant.sf")

#####Read in tx2gene#####
tx2gene <- read.csv("tx2gene.gencode.v28.csv",header=TRUE)

txi <- tximport(files, type="salmon", tx2gene=tx2gene)
print(head(txi$counts))
sampleTable <- data.frame(condition = factor(c("Untreated","Dex","Untreated","Dex","Untreated","Dex","Untreated","Dex")))
rownames(sampleTable) <- colnames(txi$counts)


ds <- DESeqDataSetFromTximport(txi,
                                colData = sampleTable,
                                design = ~ condition)
#Execute the Deseq pipeline
dds <- DESeq(ds)

## Get differential expression results
res <- results(dds)
table(res$padj<0.02)

## Order by adjusted p-value
res <- res[order(res$padj), ]

## Merge with normalized count data
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata)[1] <- "Gene ID"
head(resdata)

## Write results
write.table(resdata, file='salmon_deseq.tsv', quote=FALSE, sep='\t', col.names = NA)
