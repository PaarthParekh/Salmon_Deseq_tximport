#!/bin/bash
for fn in ./fastq_files/SRR10395{17..21};
do
samp=`basename ${fn}`
echo "Processing sample ${samp}"
salmon quant -i transcript_index -l A \
         -1 ./fastq_files/${samp}_1.fastq \
         -2 ./fastq_files/${samp}_2.fastq \
		 -p 6 --validateMappings --gcBias --numGibbsSample 20 -o ./salmon_out/${samp}_quant
done
