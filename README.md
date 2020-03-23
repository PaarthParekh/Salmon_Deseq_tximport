# Salmon_Deseq_tximport
<br/>
## Files download <br/>
To download the datasets, go to https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE52778 to get the raw reads. Download the files, after downloading the SraAccList.txt from the geo website, using the command:<br/>
prefetch $(<SraAccList.txt)<br/>
<br/>
## Processing the files and getting the count data <br/>
After downloading the sets, use the command:<br/>
fastq-dump --split-files * <br/>
To get the paired reads in the file. After getting the paired reads use salmon_code.sh to process all the split files.<br/>
<br/>
## Deseq output<br/>
Run the Deseq.R to get the Deseq output in the current working directory as salmon_deseq.tsv as a tab seperated file
