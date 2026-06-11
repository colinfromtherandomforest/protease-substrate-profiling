#list all files in directoty and loop over them
FILES=(*gz)
REF=$1

#Argument box
echo "Reference input: $REF";

echo "" > total_reads.txt

#looping through files
for f in ${FILES[@]}
do 

	echo "Aligning ... "$f 

	#----Align fastq files to reference
	bwa index $REF
	bwa mem $REF $f > ${f:0:22}.sam
	#----Convert SAM to BAM files
	samtools view -Sb -o ${f:0:22}.bam ${f:0:22}.sam
	samtools sort -o ${f:0:22}.sorted.bam ${f:0:22}.bam
	#----Filtering reads with MAPQ less than 20, secondary alignments, and CIGAR base length less than 220
	samtools view -b --min-qlen 200 -F2048 -F4 -q20 ${f:0:22}.sorted.bam > ${f:0:22}.mapped_q20.bam  
	rm ${f:0:22}.sorted.bam
	rm ${f:0:22}.sam
	rm ${f:0:22}.bam
	#----extract sequence and MD tag into txt file
	samtools view ${f:0:22}.mapped_q20.bam | awk -v OFS="," '{print substr($13,6,100),$6,$10}' | awk '{gsub(/\^/, "A", $1); print $1}' >> ${f:0:22}_match_seq2.txt 
	MATCH_READS=$(cat ${f:0:22}_match_seq2.txt | wc -l )

	echo "${f:0:22}_match_seq2.txt : "$MATCH_READS >> total_reads.txt
	echo "${f:0:22} - Reads passing MAPQ20 : "$MATCH_READS
	echo "####################" 
done

echo "*** Alignment complete ***"
echo ''
echo "---Counting AA---"
#run python script for translation and counting aa
python3 count_aa_xCFH016.py $(pwd)

echo "---- run completed, tidying files"
mkdir alignment_output
mkdir summary_output
mv *bam alignment_output
mv *.csv summary_output
mv total_reads.txt summary_output


