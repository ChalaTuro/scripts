#!/usr/bin/env python
# coding: utf-8
# Rahmat data nalysis in ipynb jupyter 
import os
get_ipython().system('ls  # magic command ')


get_ipython().run_cell_magic('bash', '', 'for i in a b c;\ndo \necho $i\ndone ')

#pwd= /resm/Rahmat
## The data fouind to contain contamination from kiwi chloroplast genome
## as such it need to be cleaned from the plant sequences using reference isolate as a baits
## use mirabiat to check if there is a contamination in the reads 
## assembled scaffolds/contigs used for testing 


get_ipython().run_cell_magic('bash', '', 'spades.py -k 21, 33, 55, 77 --isolate -1 1B_1.filtered.fq.gz  -2 1B_2.filtered.fq.gz -s 1B_1.unpaired.fq.gz -o 1BAssembly\nspades.py -k 21, 33, 55, 77 --isolate -1 2B_1.filtered.fq.gz  -2 2B_2.filtered.fq.gz -s 2B_1.unpaired.fq.gz -o 2BAssembly\nspades.py -k 21, 33, 55, 77 --isolate -1 3B_1.filtered.fq.gz  -2 3B_2.filtered.fq.gz -s 3B_1.unpaired.fq.gz -o 3BAssembly')

# The assemblies of the genome  result 
perl ~/Scripts/n50calc_James.pl 1BAssembly/scaffolds.fasta | grep "total length"
perl ~/Scripts/n50calc_James.pl 1BAssembly/scaffolds.fasta | grep "total length"
perl ~/Scripts/n50calc_James.pl 1BAssembly/scaffolds.fasta | grep "total length" 


# This asseblies looks larger than Ptt genome and can not be Ptm
#Disscussed with Lincoln and  confirmed  to be Ptm from  Rahmat
# CONTAMINATION ?
# Use the reference SG1 as Bait to check if there is contamination


get_ipython().run_cell_magic('bash', '', 'cp /resm/Rahmat/GCA_900231935.1_Ptm_SG1_genomic.fasta .\n#cp /resm/Rahmat/1BAssembly/scaffolds.fasta .\n#mirabait -b GCA_900231935.1_Ptm_SG1_genomic.fasta  scaffolds.fasta -I -t 10 \nmirabait -b GCA_900231935.1_Ptm_SG1_genomic.fasta 1BAssembly/scaffolds.fasta -I -t 10\nrm scaffolds.fasta \ncp ../2BAssembly/scaffolds.fasta  .\nmirabait -b GCA_900231935.1_Ptm_SG1_genomic.fasta  scaffolds.fasta -I -t 10 \nrm scaffolds.fasta\ncp .. 3BAssembly/scaffolds.fasta  .\nmirabait  -b GCA_900231935.1_Ptm_SG1_genomic.fasta  scaffolds.fasta -I -t 10')


get_ipython().run_cell_magic('bash', '', '#  Using custem perl scipt determine the size of the comtaminants from the asemblies\nperl ~/Scripts/n50calc_James.pl mirabiat/1B/1B_bait_miss_scaffolds.fasta | grep "total length"\nperl ~/Scripts/n50calc_James.pl mirabiat/2B/2B_bait_miss_scaffolds.fasta | grep "total length"\nperl ~/Scripts/n50calc_James.pl mirabiat/3B/3B_bait_miss_scaffolds.fasta | grep "total length"')


## The longest assembled sequence from contaminat is over 20kb
# Blast to NCBI indicats the contaminat is from Actinidia deliciosa. 
# so need to clean up the sequences and re-run the assemblies 

# install and use bbSplit to classify the genome based on ref isolate SG1

get_ipython().run_cell_magic('bash', '', '/nanodata/apps2/bbmap/bbsplit.sh qin=33 qout=33 in=1B_1.fq  in2=1B_2.fq ref=SG1.fasta basename=out1Bo%_#.fq outu1=rubish1B_1.fq outu2=rubish1B_2.fq \n/nanodata/apps2/bbmap/bbsplit.sh qin=33 qout=33 in=2B_1.fq.gz  in2=2B_2.fq.gz ref=SG1.fasta basename=2B_cleano%_#.fq outu1=2B_rubish_1.fq outu2=2B_rubbish_2.fq\n/nanodata/apps2/bbmap/bbsplit.sh qin=33 qout=33 in=3B_1.fq.gz  in2=3B_2.fq.gz ref=SG1.fasta basename=3B_cleano%_#.fq outu1=3B_rubish_1.fq outu2=3B_rubbish_2.fq')

# Re-analyse the data using fastp:  Filter the reads at Q28 with minimum read length 25 pb and sort the data according to paired reads.  Detected adapters for PE reads  activated.


get_ipython().run_cell_magic('bash', '', 'fastp -i 1B_cleanedSG1_1.fq -o 1B_1.filtered.fq.gz -I 1B_cleanedSG1_2.fq -O 1B_2.filtered.fq.gz --unpaired1 \\ 1B_1.unpaired.fq.gz --unpaired2 1B_2.unpaired.fq.gz --detect_adapter_for_pe  -z 6 -l 25 -q 28 -h 1B_report.html\nfastp -i 2B_cleanoSG1_1.fq -o 2B_1.filtered.fq.gz -I 2B_cleanoSG1_2.fq -O 2B_2.filtered.fq.gz --unpaired1  \\ 2B_1.unpaired.fq.gz --unpaired2 2B_2.unpaired.fq.gz  --detect_adapter_for_pe  -z 6 -l 25 -q 28 -h 2B_report.html\nfastp -i 3B_cleanoSG1_1.fq -o 3B_1.filtered.fq.gz -I 3B_cleanoSG1_2.fq -O 3B_2.filtered.fq.gz --unpaired1  \\ 3B_1.unpaired.fq.gz --unpaired2 3B_2.unpaired.fq.gz --detect_adapter_for_pe  -z 6 -l 25 -q 28 -h 3B_report.html ')

#Assemble the Filtered sequences using spades tool.
# here only the 


# In[ ]:


get_ipython().run_cell_magic('bash', '', 'spades.py -k 21, 33, 55, 77 --isolate -1 1B_cleanedSG1_1.fq  -2 1B_cleanedSG1_2.fq -s 1B_1.unpaired.fq.gz -o 1BAssembly\nspades.py -k 21, 33, 55, 77 --isolate -1 2B_cleanedSG1_1.fq  -2 2B_cleanedSG1_2.fq -s 2B_1.unpaired.fq.gz -o 2BAssembly\nspades.py -k 21, 33, 55, 77 --isolate -1 3B_cleanedSG1_1.fq  -2 3B_cleanedSG1_2.fq -s 3B_1.unpaired.fq.gz -o 3BAssembly\n')


# In[ ]:


get_ipython().run_cell_magic('bash', '', '#  Using custem perl scipt determine the size of the comtaminants from the asemblies\nperl ~/Scripts/n50calc_James.pl mirabiat/1B/1B_bait_miss_scaffolds.fasta | grep "total length"\nperl ~/Scripts/n50calc_James.pl mirabiat/2B/2B_bait_miss_scaffolds.fasta | grep "total length"\nperl ~/Scripts/n50calc_James.pl mirabiat/3B/3B_bait_miss_scaffolds.fasta | grep "total length"')



# this show better concordance with the refrence but slightly lower (2Mb)
## run bash command in shell as below no comment above the %%  in the cell



get_ipython().run_cell_magic('bash', '', '#pull out the chr carrying CYP51 in SG1 and plot if there is any genome rearrangement\n\n/nanodata/apps2/smashpp/smashpp -r SG1_CYP51_Chr_LT934397.1.fasta \\\n    -t ../1BAssembly/1B_cleandedAssembly_scaffolds.fasta -n 12\n\n###\n/nanodata/apps2/smashpp/smashpp -r SG1_CYP51_Chr_LT934397.1.fasta \\\n    -t ../2BAssembly/2B_cleaned_scaffolds.fasta -n 12\n\n/nanodata/apps2/smashpp/smashpp -r SG1_CYP51_Chr_LT934397.1.fasta \\\n    -t ../3BAssembly/3B_cleanedAssembly_scaffolds.fasta -n 12\n# generate the .SVG file from the .pos file \n\n/nanodata/apps2/smashpp/smashpp -viz -rn " SG1" -o "1B"  SG1_1B.svg  \\\n    SG1_CYP51_Chr_LT934397.1.fasta.1B_cleandedAssembly_scaffolds.fasta.pos\n/nanodata/apps2/smashpp/smashpp -viz -rn " SG1" -tn "2B" -o  SG1_2B.svg  \\\n    SG1_CYP51_Chr_LT934397.1.fasta.2B_cleaned_scaffolds.fasta.pos\n/nanodata/apps2/smashpp/smashpp -viz -rn " SG1" -tn  "3B" -o  SG1_3B.svg  \\\n    SG1_CYP51_Chr_LT934397.1.fasta.3B_cleanedAssembly_scaffolds.fasta.pos')

get_ipython().run_cell_magic('bash', '', '#Order the assemblies as per SG1 ref\ncd /resm/Rahmat/compaative\n#python3 ~/RaGOO/ragoo.py ../1BAssembly/1B_cleandedAssembly_scaffolds.fasta  ../SG1.fasta\n#python3 ~/RaGOO/ragoo.py ../2BAssembly/2B_cleaned_scaffolds.fasta  ../SG1.fasta\npython3 ~/RaGOO/ragoo.py ../3BAssembly/3B_cleanedAssembly_scaffolds.fasta  ../SG1.fasta\ncd /resm/Rahmat')

get_ipython().run_cell_magic('bash', '', 'cd /resm/Rahmat/compaative\n/nanodata/apps2/smashpp/smashpp -r ../SG1.fasta \\\n    -t ../1BAssembly/1B_cleandedAssembly_scaffolds.fasta -n 12\ncd /resm/Rahmat')

get_ipython().run_cell_magic('bash', '', 'cd /resm/Rahmat/compaative\n/nanodata/apps2/smashpp/smashpp -viz -rn " SG1" "1B" -o  All_SG1_1B.ordered.svg  \\\n    SG1.fasta.1B_cleandedAssembly_scaffolds.fasta.pos\ncd /resm/Rahmat')

get_ipython().run_cell_magic('bash', '', '## install and run sibeliaz whole genome alignment \n## rename the fasta header to distinguish the isolates\ncd /resm/Rahmat/sibeliaz\n/usr/bin/bin/sibeliaz  *.fasta   # use default parameters to run ')

get_ipython().run_cell_magic('bash', '', '# does not show great re-arangements \n#Maybe useful to use Altive to run comparative analysis simeltanosely\n# Blast analysis CYP51 +- 600bp  and alignment to SG1 region showed no difference!\n# Is the host selection pressure needed for mutation to occur? or Natural environment + host ?\n# https://alitv.readthedocs.io/en/latest/index.html\n\n# cd ~/APPS/AliTV-perl-interface/rahmat\n#../../bin/alitv.pl *.fasta --project 3BSG1\n#visualize the result at: https://alitvteam.github.io/AliTV/d3/AliTV.html')

get_ipython().run_cell_magic('bash', '', '# run  variant analysis \nbwa index  SG1_masked.fasta\nbwa mem SG1_masked.fasta 1B_cleanedSG1_1.fq 1B_cleanedSG1_2.fq | samtools view -Sb - -o 1B.bam\nbwa mem SG1_masked.fasta 2B_cleanedSG1_1.fq 2B_cleanedSG1_2.fq | samtools view -Sb - -o 2B.bam\nbwa mem SG1_masked.fasta 3B_cleanedSG1_1.fq 3B_cleanedSG1_2.fq | samtools view -Sb - -o 3B.bam\n# sort the alignment \nsamtools sort -O BAM -@ 8 -o 1B.sorted.bam 1B.bam\nsamtools sort -O BAM -@ 8 -o 2B.sorted.bam 2B.bam\nsamtools sort -O BAM -@ 8 -o 3B.sorted.bam 3B.bam\n#remove unsorted bam to save space \nrm -rf 1B.bam 2B.bam 3B.bam\n\n')

# mark duplicates; the need the bam to be sorted.
%%bash
# run SNP analysis 
# Pirooznia et al. Human Genomics 2014, 8:14 indicated that GATK is better in accuracy of SNP call

gatk MarkDuplicates -I 1B.bam  -O 1B_markedDuplils *.bam 
cates.bam  -M 1B.marked_dup_metrics.txt
gatk MarkDuplicates -I 2B.bam  -O 2B_markedDuplicates.bam  -M 2B.marked_dup_metrics.txt
gatk MarkDuplicates -I 3B.bam  -O 3B_markedDuplicates.bam  -M 3B.marked_dup_metrics.txt
# remove duplicates 
gatk MarkDuplicates -I 1B.sorted.bam  -O 1B_Noduplicate.bam --REMOVE_DUPLICATES true -M 1B.metrics
gatk MarkDuplicates -I 2B.sorted.bam  -O 2B_Noduplicate.bam --REMOVE_DUPLICATES true -M 2B.metrics
gatk MarkDuplicates -I 3B.sorted.bam  -O 3B_Noduplicate.bam --REMOVE_DUPLICATES true -M 3B.metrics
# recalibrate 
# this did not worked 
gatk BaseRecalibrator -I 1B_markedDuplicates.bam -R SG1_masked.fasta  -O 1B.recalbrated.table
gatk BaseRecalibrator -I 2B_markedDuplicates.bam -R SG1_masked.fasta  -O 2B.recalbrated.table
gatk BaseRecalibrator -I 3B_markedDuplicates.bam -R SG1_masked.fasta  -O 3B.recalbrated.table

#apply recalibration #not worked/used as I have no known indels

gatk ApplyBQSR -R  SG1_masked.fasta -I 1B_markedDuplicates.bam     --bqsr-recal-file 1B.recalbrated.table -O 1BCalibrated.bam

gatk ApplyBQSR -R  SG1_masked.fasta -I 1B_markedDuplicates.bam     --bqsr-recal-file 2B.recalbrated.table -O 2BCalibrated.bam

gatk ApplyBQSR -R  SG1_masked.fasta -I 1B_markedDuplicates.bam     --bqsr-recal-file 3B.recalbrated.table -O 3BCalibrated.bam

# add read group 
gatk AddOrReplaceReadGroups        -I 1B_Noduplicate.bam        -O 1B_Noduplicate-RG.bam        --RGID 4        --RGLB B11        --RGPL ILLUMINA        --RGPU UNIT1        --RGSM B1
gatk AddOrReplaceReadGroups        -I  2B_Noduplicate.bam        -O 2B_Noduplicate-RG.bam        --RGID 4        --RGLB B22        --RGPL ILLUMINA        --RGPU UNIT1        --RGSM B2
gatk AddOrReplaceReadGroups        -I 3B_Noduplicate.bam        -O 3B_Noduplicate-RG.bam        --RGID 4        --RGLB B33        --RGPL ILLUMINA        --RGPU UNIT1        --RGSM B3
# index bamfiles
samtools index 1B_Noduplicate-RG.bam
samtools index  2B_Noduplicate-RG.bam
samtools index  3B_Noduplicate-RG.bam
# create sequence dic of the ref 
gatk CreateSequenceDictionary -R  SG1_masked.fasta -O SG1_masked.dict

# call the SNP 
gatk HaplotypeCaller --java-options "-Xmx4G" -R SG1_masked.fasta -I 1B_Noduplicate-RG.bam -O 1B.vcf.gz 
gatk HaplotypeCaller --java-options "-Xmx4G" -R SG1_masked.fasta -I 2B_Noduplicate-RG.bam -O 2B.vcf.gz 
gatk HaplotypeCaller --java-options "-Xmx6G" -R SG1_masked.fasta -I 3B_Noduplicate-RG.bam -O 3B.vcf.gz 
# not the java option did not worked (i.e number of threads)

# Samtools 
##Call the SNP 
ref=SG1_masked.fasta
for f in *BCalibrated.bam
do
 gatk --java-options "-Xmx8g" HaplotypeCaller     -R $ref     -I $f    -O $f.vcf.gz    -ERC GVCF

echo " Partial analysis completed !"
echo ""
cd /resm/Rahmat


# Structural variATION ANALYSIS 
#USE lUMPY TO RUN TE STRUCTURAL VARIATION ANALYSIS
speedseq  align -o 1B -t 12 -R "@RG\tID:id\tSM:1B\tLB:1B" SG1_masked.fasta 1B_cleanedSG1_1.fq 1B_cleanedSG1_2.fq
lumpyexpress     -B 1B.bam     -S 1B.splitters.bam     -D 1B.discordants.bam     -o 1B.vcf
    
speedseq  align -o 1B -t 12 -R "@RG\tID:id\tSM:1B\tLB:2B" SG1_masked.fasta 2B_cleanedSG1_1.fq  2B_cleanedSG1_2.fq
lumpyexpress     -B 2B.bam     -S 2B.splitters.bam     -D 2B.discordants.bam     -o 2B.vcf
speedseq  align -o 1B -t 12 -R "@RG\tID:id\tSM:1B\tLB:2B" SG1_masked.fasta 3B_cleanedSG1_1.fq   3B_cleanedSG1_2.fq
lumpyexpress     -B 3B.bam     -S 3B.splitters.bam     -D 3B.discordants.bam     -o 3B.vcf
## proceed to analysis of VCF

