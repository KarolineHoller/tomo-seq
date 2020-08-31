#This mapping script has been originally written by Bastiaan Spanjaard.

#Load an environment that can access the Software of STAR
source /environment_path
# Add UMIs and BCs from the R2 to the fastq file of R1
awk 'NR%4==2{printf "_UMI_%s_BC_%s\n\n\n\n", substr($0, 0, 6), substr($0, 7, 8)}' file_R2.fastq > UMIBC.txt
#writes UMIs and barcodes from file_R2.fastq into a txt file
paste -d " " <(paste -d "" <(awk '{print $1}' file_R1.fastq) UMIBC.txt) <(awk '{print $2}' file_R1.fastq) > extended_file_R1.fastq
#adds UMI and barcode next to the paired read and creates a new fastq file
rm UMIBC.txt #removes the txt file

# Map the fastq file
STAR --runThreadN  --genomeDir /genome_path --quantMode TranscriptomeSAM --readFilesIn extended_file_R1.fastq
samtools view Aligned.toTranscriptome.out.bam > Transcr_NCBI_UMIBC.sam
#converts the binary .bam into a human readable .sam file

# Extract transcripts
join -1 1 -2 3 <(sort -k1,1 /gene_name_conversion_path/translator.tsv) <(sort -k3,3 Transcr_NCBI_UMIBC.sam) > Transcr_genes.sam
#converts gene identifier into gene names
awk '{if($4==0){split($3, a, "_"); printf("%s\t%s\t%s\tsense\n", a[3], a[5], $2)}if($4==16){split($3, a, "_"); printf("%s\t%s\t%s\tantisense\n", a[3], a[5], $2)}}' Transcr_genes.sam | sort -k3,3
-k4,4 -k2,2 -k1,1 | uniq | cut -f2,3,4 | uniq -c > genes_transcripts.tsv
#extracts gene counts per section barcode
join -1 2 -2 2 <(sort -k2,2 genes_transcripts.tsv) <(sort -k2,2 /path/section_barcodes.csv) > genes_transcripts_bc.tsv
# adds the section number to the respective section barcode
