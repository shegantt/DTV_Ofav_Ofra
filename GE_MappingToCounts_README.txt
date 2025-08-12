
#Tag-based RNA-seq reads processing pipeline
#last revised by Shelby Gantt on 12 Nov 2024
#adapted from Dan's codes (https://github.com/wuitchik/MPCC_2018/blob/main/Host_analyses/Mapping/Astrangia/Astrangia_mapping.sh)

#for all code, run directly in terminal pasting all consecutive rows of code at once.

#More info on the packages
#fastx_toolkit: http://hannonlab.cshl.edu/fastx_toolkit/download.html
#bowtie2: http://bowtie-bio.sourceforge.net/index.shtml 
#pipeline came from https://github.com/z0on/tag-based_RNAseq


#This code is for already cleaned and trimmed sequences
#start by concatenating genomes

#to request a more powerful and longer session than the login node
qrsh -l h_rt=12:00:00 -pe omp 4

cd /projectnb/~directory path~  #directory with genomes
cd /projectnb/~directory path~ #directory with trimmed/cleaned files


#-----Grab data and scripts-----

#concatenate host and symbiont genomes (from Dan)
#cat coral_refgenome.fasta algae_reftranscriptome.fasta > holobiont.fasta

cat /projectnb/~directory path~/Ofav/Genomes/Orbicella_faveolata_gen_17.scaffolds.fa /projectnb/~directory path~/Ofav/Genomes/Dtrenchii_CCMP2556_ASSEMBLY_fasta.fasta >  holobiont.fasta

#count the number of contigs in your fasta file
grep '>' holobiont.fasta | wc -l 
#29188

#concatenate host and symbiont annotation files (from Dan)
#did not do. decided to split host and sym files before annotation/counts
#cat coral.gff algae_ref.gff > holobiont.gff3
cat /projectnb/~directory path~/Orbicella_faveolata_gen_17.gff3 /projectnb/~directory path~/Ofav/Genomes/Dtrenchii_CCMP2556_ANNOT_gff.gff3 > holobiont.gff3

#-----Mapping-----

############################# Bowtie2 ###################

#Make a new job file called BowtieBuild.sh and put the following lines in it.

#create bowtie2 index for reference genome
module load bowtie2
module load htslib/1.19.1
module load samtools/1.19.2  	
module load bbmap		
module load fastx-toolkit	

bowtie2-build holobiont.fasta holobiont.fasta 
samtools faidx holobiont.fasta 
ls

#then call job file BowtieBuild.sh to build reference
qsub BowtieBuild.sh

#create batch file and launch the mapping process for all files at the same time
#run this code directly in terminal
module load bowtie2 		
module load htslib/1.19.1
module load samtools/1.19.2  	
module load bbmap		
module load fastx-toolkit	

for file in *.trim
	do echo "bowtie2 -x holobiont.fasta -U $file --local -p 4 -S ${file/.trim/}.sam">> maps
done
nano maps

# create a batch job (first line) and submit (second line) to execute all commands written to file ‘maps’ 
scc6_qsub_launcher.py -N maps -P davieslab -jobsfile maps
qsub maps_array.qsub


#make directory and move map* files to new directory
mkdir parallel_outputs
mv maps* parallel_outputs

#------- Sort Sam files ---------

#make directory for SAM files
mkdir SamFiles

#transfer all SAM files to this directory
mv *.sam SamFiles

#------Split host and Symbiont reads into separate Sam files---- 
#Run all of splitting code directly in terminal

#from Hannah Aichelman's GitHub for her Oculina_Host_Sym_GE work
#'ofavscaf\|SO:unsorted' and 'Dtrenchii\|SO:unsorted' are specific identifiers in the host and sym genome files.

#run the following in the terminal
cd SamFiles
for x in *.sam; do 
	echo $x; grep 'ofavscaf\|SO:unsorted' $x &> ${x//sam/host.sam}; grep 'Dtrenchii\|SO:unsorted' $x &> ${x//sam/Sym.sam}; done

#make directories for Sym and Ofav files and move files into these directories
mkdir HostSams
mkdir SymSams
mv *host.sam HostSams
mv *Sym.sam SymSams

#now process each in those directories as follows
cd /projectnb/~directory path~/Ofav/SamFiles/HostSams
cd /projectnb/~directory path~/Ofav/SamFiles/SymSams

#------Gene Count Tables ----------

#-----Sort the sam alignment files and convert to bams-----
module load htslib/1.19.1
module load samtools/1.19.2

for file in *.sam
	do echo "samtools sort -n -O bam -o ${file/.sam/}.bam $file" >> sortConvert
done 
nano sortConvert

scc6_qsub_launcher.py -N sortconverting -P davieslab -jobsfile sortConvert
qsub sortconverting_array.qsub

mkdir parallel_outputs
mv sortconverting* parallel_outputs


module load subread

# If running Host & Sym together choose the GFF holobiont.gff3
# If running separate host and sym use: 
# Orbicella_faveolata_gen_17.gff3 for host files
# and Dtrenchii_CCMP2556_ANNOT_gff.gff3 for sym files
#I did both ways to confirm that we get the same count table whether they are combined or not, and we do.

MY_GFF="holobiont.gff3"; GENE_ID="ID"
featureCounts -a $MY_GFF -t gene -g $GENE_ID -o feature_counts_out.txt -T 32 *.bam

#get alignment counts
for file in *.bam
	do echo "samtools flagstat $file > ${file/.bam/}_flagStats.txt" >> getInitialAlignment
done

scc6_qsub_launcher.py -N samtools_alignment_counts -P davieslab -jobsfile getInitialAlignment
qsub samtools_alignment_counts_array.qsub

mv samtools* parallel_outputs

#-----Making raw read count table-----
module purge
module load python2
module load htseq/0.11.0

# If running Host & Sym together choose the GFF holobiont.gff3
# If running separate host and sym use: 
# Orbicella_faveolata_gen_17.gff3 for host files
# and Dtrenchii_CCMP2556_ANNOT_gff.gff3 for sym files

MY_GFF="holobiont.gff3"; GENE_ID="ID"

for file in *.sam
	do echo "htseq-count -t gene -i ID -m intersection-nonempty --stranded=no $file $MY_GFF > ${file/sam/counts.txt}" >> doCounts
done

scc6_qsub_launcher.py -N counts -P davieslab -jobsfile doCounts
qsub counts_array.qsub 


#-----Compile to make big host and Sym specific expression tables-----

#for host & sym counts together
expression_compiler.pl *.counts.txt > holobiont_counts.txt

#For split SAM host counts
expression_compiler.pl *.counts.txt > Ofav_counts.txt

#For split SAM sym counts
expression_compiler.pl *.counts.txt > Dtrenchii_counts.txt


#-----Making a table compiling the number of counts in each step-----
#Dan's code for checking counts - follow the text to see in which directory with the appropriate files 
#each needs to be run, then move all tsv files to the same directory to compile count estimates across steps.

>mapped_count.tsv
for file in *_flagStats.txt
do pp=$(grep "mapped" $file | head -n 1)
 echo -e "$file\t$pp" |\
 awk '{split($1, a, "_flagStats.txt")
 print a[1]"\t"$2"\tmapped"}' >> mapped_count.tsv
 done

# raw
wc -l *.fastq |\
	awk '{split($2, a, ".fastq")
	print a[1]"\t"$1/4"\rawCounts"}' |\
grep -v total > raw_read_counts.tsv 

cd /projectnb/davieslab/shelby/Ofav

# trimmed
wc -l *.trim |\
	awk '{split($2, a, ".trim")
	print a[1]"\t"$1/4"\ttrimmedCounts"}' |\
grep -v total > trimmed_read_counts.tsv 

mv trimmed* DtrenchiiOnlySams

cd /projectnb/davieslab/shelby/Ofav/DtrenchiiOnlySams

# Compile pipeline counts for supps
cat *.tsv > final_pipeline_counts.txt


#export counts tables from scc to desktop for R analysis or run on scc R
