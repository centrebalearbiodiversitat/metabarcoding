#! /bin/bash


### PREPARE FOLDERS AND STRUCTURE FOR DATA OBTAINED


mkdir quality_control
mkdir reads
mkdir reads/paired_reads
mkdir reads/filtered
mkdir quality_control/fastqc_raw_reads
mkdir quality_control/filtered_reads_fastqc


### PRIMER REMOVAL AND QUALITY CONTROL


ls | grep -E '_R[12].fastq$' | sed -E 's/_R[12].fastq$//'>SLIST


##_PRIMER_F=''
PRIMER_F='GxxxxTYCCGCG'
 
##ArR5_PRIMER_R='GTRATIGCICCIGCIARIACIGG'
PRIMER_R='GTRATxxxGCGARGACGGG'




while read -r S
do
#        fastqc ${S}_R1.fastq ${S}_R2.fastq -o quality_control/fastqc_raw_reads
        cutadapt ${S}_R1.fastq ${S}_R2.fastq -g $PRIMER_F -G $PRIMER_R -o reads/filtered/${S}_R1_filtered.fastq -p reads/filtered/${S}_R2_filtered.fastq --discard-untrimmed -e 0.1 
#	fastqc reads/filtered/${S}_R1_filtered.fastq reads/filtered/${S}_R2_filtered.fastq -o quality_control/filtered_reads_fastqc


done < SLIST


### PAIRED-END READS MERGING


while read -r S
do
	pear -f reads/filtered/${S}_R1_filtered.fastq -r reads/filtered/${S}_R2_filtered.fastq -o reads/paired_reads/${S} -j 12 -t 4 -n 32
done < SLIST




#### QUALITY CONTROL AFTER MERGING
mkdir quality_control/paired_reads


#while read -r S
#do
#	fastqc reads/paired_reads/${S}.assembled.fastq -o quality_control/paired_reads
#done < SLIST


#python -m multiqc quality_control/paired_reads/* -o quality_control/paired_reads


rm SLIST




mkdir reads/trimmed_paired_reads


### MERGED READS TRIMMERING 


cd reads/paired_reads
 
ls | grep -E '.assembled.fastq$' | sed -E 's/.assembled.fastq$//' > SLIST
while read -r S
do
        TrimmomaticSE -threads 8 -phred33 ${S}.assembled.fastq ${S}_paired_trimmed.fastq SLIDINGWINDOW:6:20 MINLEN:250
done < SLIST


rm SLIST
mv *_paired_trimmed.fastq ~/Project_folder/reads/trimmed_paired_reads/


cd ..


### CLUSTERING AND SORTING 


mkdir chimera_removal cdhit_out vsearch_out


cd trimmed_paired_reads


ls| grep -E '_trimmed.fastq$'| sed -E 's/.fastq$//' > SLIST
while read -r S
do
	sample=${S:10:6}
        ## FASTQ TO FASTA
        seqkit fq2fa ${S}.fastq -o ${S}.fasta


        ## CHIMERA FILTERING, SEARCHES ALL OF THE AMPLICON SEQUENCES VARIANT AGAINST ONE ANOTHER.  
        vsearch --threads 4 --uchime_denovo ${S}.fasta \
	--chimeras ~/Project_folder/reads/chimera_removal/${S}_unique_denovo.chimera \
        --nonchimeras ~/Project_folder/reads/chimera_removal/${S}_unique_denovo.good  


	#### CLUSTERING 1


	## DE NOVO VSEARCH CLUSTERING ##0.99 identity, sorting, adding cluster id, size and relabeling. 
	vsearch --threads 4 --cluster_fast ~/Project_folder/reads/chimera_removal/${S}_unique_denovo.good \
	 --sizeout --relabel ${sample}_hap_ --clusterout_id \
	--centroids ~/Project_folder/reads/vsearch_out/${S}_unique_denovo.cluster.centroids \
	--id 0.99 --sizeorder --clusterout_id --consout ~/Ports_IB/Pons3_PrimerSortedDemultiplexed_2018430/ArF5_ArR5/reads/vsearch_out/${S}_consensus_cluster 


	#### CLUSTERING 2
	
        ## DE NOVO BASED OTU CLUSTERING ## -c sequence identity threshold . -T threads -d description len -sf sort bt size -s length difference cutoff
 
        cdhit -i ~/Ports_IB/Pons3_PrimerSortedDemultiplexed_2018430/ArF5_ArR5/reads/chimera_removal/${S}_unique_denovo.good \
	-o ~/Ports_IB/Pons3_PrimerSortedDemultiplexed_2018430/ArF5_ArR5/reads/cdhit_out/${S}.clust \
	-c 0.99 -T 4 -d 35 -sf 1 -sc 0 -s 0.9


done < SLIST


rm SLIST


cd ..
cd vsearch_out


mkdir cluster_size1_removed
mkdir cluster_size2_removed
mkdir cluster_size3_removed
mkdir cluster_size4_removed
mkdir cluster_size5_removed
mkdir cluster_size6_removed


### REMOVE CLUSTERS WITH SMALLER SIZE THAN [1-6].  FROM CLUSTERS FOUND THROUGH VSEARCH 




ls | grep -E '\.centroids$' | sed -E 's/.centroids$//' > SLIST


while read -r S  
do
	cf.pl -w *${S}.centroids Fasta > ${S}.inline.fasta
	awk '/^>/{header=$0; getline seq; if (header !~ /;size=[0-9]*[0-1]/) {print header; print seq}}' ${S}.inline.fasta > ${S}.removedsize1.fasta
        awk '/^>/{header=$0; getline seq; if (header !~ /;size=[0-9]*[0-2]/) {print header; print seq}}' ${S}.inline.fasta > ${S}.removedsize2.fasta
        awk '/^>/{header=$0; getline seq; if (header !~ /;size=[0-9]*[0-3]/) {print header; print seq}}' ${S}.inline.fasta > ${S}.removedsize3.fasta
        awk '/^>/{header=$0; getline seq; if (header !~ /;size=[0-9]*[0-4]/) {print header; print seq}}' ${S}.inline.fasta > ${S}.removedsize4.fasta
        awk '/^>/{header=$0; getline seq; if (header !~ /;size=[0-9]*[0-5]/) {print header; print seq}}' ${S}.inline.fasta > ${S}.removedsize5.fasta
        awk '/^>/{header=$0; getline seq; if (header !~ /;size=[0-9]*[0-6]/) {print header; print seq}}' ${S}.inline.fasta > ${S}.removedsize6.fasta


	sed '/^>/ s/;/_/g; s/=/_/g' ${S}.removedsize1.fasta  > ${S}.good1.fasta
	sed '/^>/ s/;/_/g; s/=/_/g' ${S}.removedsize2.fasta  > ${S}.good2.fasta
        sed '/^>/ s/;/_/g; s/=/_/g' ${S}.removedsize3.fasta  > ${S}.good3.fasta
        sed '/^>/ s/;/_/g; s/=/_/g' ${S}.removedsize4.fasta  > ${S}.good4.fasta
        sed '/^>/ s/;/_/g; s/=/_/g' ${S}.removedsize5.fasta  > ${S}.good5.fasta
        sed '/^>/ s/;/_/g; s/=/_/g' ${S}.removedsize6.fasta  > ${S}.good6.fasta
	rm ${S}.inline.fasta
	rm ${S}.removedsize[0-6].fasta
	mv ${S}.good1.fasta cluster_size1_removed
        mv ${S}.good2.fasta cluster_size2_removed
        mv ${S}.good3.fasta cluster_size3_removed
        mv ${S}.good4.fasta cluster_size4_removed
        mv ${S}.good5.fasta cluster_size5_removed
	mv ${S}.good6.fasta cluster_size6_removed
 
done < SLIST
rm SLIST


cd ../..


cd reads/cdhit_out


# Adding TAGs to the Haplotyes found by cdhit. 


ls | grep -E '.clust.clstr$' > SLIST 
ls | grep -E '.clust.clstr$' | cut -c 10-16 > NLIST
while read -r S
do
        TotalClusters=$(grep -c '>Cluster' "$S")
                for ((i=0; i<${TotalClusters}; i++)); do
                        j=$((i+1))
			while read -r N
			do
				substring=${N}
				if [[ $S == *"$substring"* ]]; then
                        		lines=$(awk -v i="$i" -v j="$j" '/^>Cluster /{flag=0} /^>Cluster '"$i"'/{flag=1} /^>Cluster '"$j"'/ && flag {exit} flag{count++} END {print count}' "$S")
                        		linesExcludingPattern=$((lines-2))
                        		echo " $substring Hap ${i} -- $((linesExcludingPattern+1)) sequences " >> ${S}.txt
				fi
			done < NLIST
                done
done < SLIST
rm SLIST 




rm NLIST
####### CLUSTERS WITH POUETS NAME


cat *.txt > haploid_sample_size.txt


mkdir clust_names
ls | grep -E '\.clust$'| cut -c 10-16 > SLIST 
ls | grep -E '_paired_trimmed.clust$' > MLIST
while read -r M
do 
        while read -r S
        do
                substring="${S}"
                if [[ $M == *"$substring"* ]]; then
                        awk -v text="$substring" '/>/ { $0 = $0 " " text}{print}' "$M" > tmp.txt && mv tmp.txt "$M.cluster_sampleName"
                        awk -v text="$substring" '/Cluster/ { $0 = $0 " " text}{print}' "$M.clstr" > tmp.txt && mv tmp.txt "$M.stats"
                        awk -v text="$substring" '/Cluster/ { $0 = $0 " " text}{print}' "$M.clstr.txt" > tmp.txt && mv tmp.txt "$M.hapsize"
                fi
        done < SLIST
done < MLIST


rm SLIST MLIST


###### CAT al .test1


cp *.cluster_sampleName clust_names
mv *.hapsize clust_names
mv *.stats clust_names


cat *.cluster_sampleName > Ar.aligned
