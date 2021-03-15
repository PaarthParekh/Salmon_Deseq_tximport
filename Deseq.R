library(tximport)
library(readr)
library(DESeq2)

######To Generate the Transcripts ids to geneids ids needed for Deseq2 #####
###  We generate the tx2gene.gencode.csv file mapping each transcript id to Gene id #########

#library(tximportData)
#library(GenomicFeatures)
#txdb <- makeTxDbFromGFF(file="gencode.v28.annotation.gff3.gz")
#saveDb(x=txdb, file = "gencode.v28.annotation.TxDb")

#k <- keys(txdb, keytype = "TXNAME")
#tx2gene <- select(txdb, k, "GENEID", "TXNAME")
#write.table(tx2gene, "tx2gene.gencode.v28.csv", sep = ",", row.names = FALSE)

#get all the counts file generated from Salmon for all the samples 
dirs = "./salmon_out"
sys_dir <- system.file("extdata", package = "tximportData")
samples = dir(file.path(dirs))
files <- file.path(dirs, samples, "quant.sf")

#####Read in tx2gene gencode csv file containing transcript and gene ids#####
tx2gene <- read.csv("tx2gene.gencode.v28.csv",header=TRUE)

# using tximport map the transcript ids to gene ids
txi <- tximport(files, type="salmon", tx2gene=tx2gene)
print(head(txi$counts))
# define the sample conditions between which we are trying to evaluate the deseq result
sampleTable <- data.frame(condition = factor(c("Untreated","Dex","Untreated","Dex","Untreated","Dex","Untreated","Dex")))
rownames(sampleTable) <- colnames(txi$counts)

# Run the Deseq pipeline after defining to the Deseq how the data is coming in 
ds <- DESeqDataSetFromTximport(txi,colData = sampleTable,design = ~ condition)
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
