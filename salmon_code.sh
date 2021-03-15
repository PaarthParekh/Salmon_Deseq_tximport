#!/bin/bash
# for all the files (17 to 21) in the directory run the salmon quantification command to obtain the counts for each transcript  
for fn in ./fastq_files/SRR10395{17..21};
do
# this variable is to print to std out, what sample is currently running for salmon quantification
samp=`basename ${fn}`
echo "Processing sample ${samp}"
salmon quant -i transcript_index -l A \
         -1 ./fastq_files/${samp}_1.fastq \
         -2 ./fastq_files/${samp}_2.fastq \
		 -p 6 --validateMappings --gcBias --numGibbsSample 20 -o ./salmon_out/${samp}_quant
done
