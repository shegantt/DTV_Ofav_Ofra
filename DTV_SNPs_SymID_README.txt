#Nicola's DTV experiment

##first gunzipped the files
#start an interactive session
qrsh -pe omp 12

gunzip *.gz
#all files should unzip and now be .fastq

#these samples were sequenced across two lanes so need to concatenate the reps
mkdir concat
for x in *L001*fastq; do echo $x; o=${x//_*/.fastq}; l=${x//L001/L002}; cat $x $l > concat/$o; done

###first we wants to call SNPs to make sure reps are reps and there are no clones that we were not aware of
#trim all fastq files
for F in *.fastq; do
echo "tagseq_clipper.pl $F | cutadapt - -a AAAAAAAA -a AGATCGG -q 15 -m 25 -o ${F/.fq/}.trim" >>clean;
done

nano clean
#you will see that each file gets a line that looks like this:
tagseq_clipper.pl KG12-2.fastq | cutadapt - -a AAAAAAAA -a AGATCGG -q 15 -m 25 -o KG12-2.fastq.trim

scc6_qsub_launcher.py -N trim -P coral -M daviessw@gmail.com -j y -h_rt 24:00:00 -jobsfile clean

module load cutadapt
qsub trim_array.qsub

#prep genomes/transcriptomes for mapping (doing symbionts here too for future mapping to determine sym type)
module load bowtie2
bowtie2-build Ofav_GCF_002042975.1_ofav_dov_v1_rna.fna Ofav_GCF_002042975.1_ofav_dov_v1_rna.fna
bowtie2-build symABCD.fasta symABCD.fasta

#Now map all .trim files
2bRAD_bowtie2_launch.pl '\.trim$' Ofav_GCF_002042975.1_ofav_dov_v1_rna.fna >bowtie2_forsnp

scc6_qsub_launcher.py -N mapping -P coral -M daviessw@gmail.com -j y -h_rt 24:00:00 -jobsfile bowtie2_forsnp

module load bowtie2
qsub mapping_array.qsub

#example line in script:
bowtie2 --no-unal --score-min L,16,1 --local -L 16 -x Ofav_GCF_002042975.1_ofav_dov_v1_rna.fna -U VE30-1_S76_L001_R1_001.fastq.trim $

#breakdown =
#--no-unal = suppress SAM records for unaligned reads
#--score-min L,16,1 [L = linear --> f(x) = 16 + 1 * x]
#--local -L 16 = Sets the length of the seed substrings to align during multiseed alignment. 
#Smaller values make alignment slower but more sensitive

#this is different from the tagseq mapping script in a few ways. the tagseq method includes these unique flags:
--no-hd = suppress SAM header lines (starting with @)
--no-sq  = suppress @SQ SAM header lines
-k = indicates -k mode. In this mode, bowtie searches for up to N distinct, valid alignments for each read. We have -k 5, meaning bowtie will search for at most 5 distinct alignments. 
The manual says that the -k mode is effective in situations where you care more about whether a read aligns, or aligns a certain number of times, rather than where exactly it originated

##now we have .sam files
#now try making bam files:
ls *bt2.sam > sams
cat sams | wc -l
#40 files

cat sams | perl -pe 's/(\S+)\.sam/samtools import Ofav_GCF_002042975.1_ofav_dov_v1_rna.fna $1\.sam $1\.unsorted\.bam && samtools sort -o $1\.sorted\.bam $1\.unsorted\.bam && picard AddOrReplaceReadGroups INPUT=$1\.sorted\.bam OUTPUT=$1\.bam RGID=group1 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=$1 && samtools index $1\.bam/' >s2b

#Example command for one file:
samtools import Ofav_GCF_002042975.1_ofav_dov_v1_rna.fna VG6-1_S84_L002_R1_001.fastq.trim.bt2.sam VG6-1_S84_L002_R1_001.fastq.trim.bt2.unsorted.bam && samtools sort -o VG6-1_S84_L002_R1_001.fastq.trim.bt2.sorted.bam VG6-1_S84_L002_R1_001.fastq.trim.bt2.unsorted.bam && picard AddOrReplaceReadGroups INPUT=VG6-1_S84_L002_R1_001.fastq.trim.bt2.sorted.bam OUTPUT=VG6-1_S84_L002_R1_001.fastq.trim.bt2.bam RGID=group1 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=VG6-1_S84_L002_R1_001.fastq.trim.bt2 && samtools index VG6-1_S84_L002_R1_001.fastq.trim.bt2.bam

scc6_qsub_launcher.py -N sam2bam -P coral -M daviessw@gmail.com -j y -h_rt 24:00:00 -jobsfile s2b

##note these module loads need to be within the qsub script
module load htslib/1.9
module load samtools/1.9
module load picard

qsub sam2bam_array.qsub

#once job finishes, remove unnecessary intermediate files
rm -f *sorted*

ls *bt2.bam > bams
cat bams | wc -l
#40

### angsd genotyping - initial quality check
module load samtools
module load picard

#first need to create a dictionary and an index (.fai) file for the reference transcriptomes
picard CreateSequenceDictionary R= Ofav_GCF_002042975.1_ofav_dov_v1_rna.fna O= Ofav_GCF_002042975.1_ofav_dov_v1_rna.fna.dict

samtools faidx Ofav_GCF_002042975.1_ofav_dov_v1_rna.fna

#--------------- population structure (based on common polymorphisms)

# F I L T E R S :
# (ALWAYS record your filter settings and explore different combinations to confirm that results are robust. )
# Suggested filters :
# -minMapQ 20 : only highly unique mappings (prob of erroneous mapping = 1%)
# -minQ 25 : only highly confident base calls
#could make this 20 to get more snps
# -minInd 35 : the site must be genotyped in at least 35 individuals (note: set this to ~ 80% of your total number of your individuals) 0.8*44 =35.2
#lower this if you want more snps
# -minIndDepth 5 : depth of at least 5 in non-missing individual
#change to 5, no lower than 5 though
# -snp_pval 1e-5 : high confidence that the SNP is not just sequencing error 
# -minMaf 0.05 : only common SNPs, with allele frequency 0.05 or more.
# Note: the last two filters are very efficient against sequencing errors but introduce bias against true rare alleles. It is OK (and even desirable) - UNLESS we want to do AFS analysis. We will generate data for AFS analysis in the next part.\
# also adding  filters against very badly non-HWE sites (such as, all calls are heterozygotes => lumped paralog situation) and sites with really bad strand bias:\

# T O   D O : \
# -GL 1 : samtools likelihood model\
# -doGlf 2 : output beagle format (for admixture)\
# -doGeno 8 : binary genotype likelihoods format (for ngsCovar => PCA)\44
# -doMajorMinor 4 -ref $GENOME_REF : infer major and minor alleles from reference (in our case it is outgroup taxon)\
# -makeMatrix 1 -doIBS 1 -doCov 1 : identity-by-state and covariance matrices based on single-read resampling (robust to variation in coverage across samples)\
# TO-DO commands can be changed to generate different files needed for different analyses as well (http://www.popgen.dk/angsd/index.php/ANGSD#Overview) 
#skipTrialellic 1 #could add this if you want to remove triallelic snps, usually a sequencing error

#Request an interactive node before you start this analysis, might get kicked off for taking up too much memory on the head node
qrsh -pe omp 12

#use these filters:
FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 25 -dosnpstat 1 -doHWE 1 -sb_pval 1e-5 -hetbias_pval 1e-5 -skipTriallelic 1 -minInd 36 -snp_pval 1e-5 -minMaf 0.05"
TODO="-doMajorMinor 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 8 -doBcf 1 -doPost 1 -doGlf 2"

module load angsd
angsd -b bams -GL 1 $FILTERS $TODO -P 1 -out nicola

Output filenames:
		->"nicola.arg"
		->"nicola.beagle.gz"
		->"nicola.mafs.gz"
		->"nicola.hwe.gz"
		->"nicola.geno.gz"
		->"nicola.snpStat.gz"
		->"nicola.ibs.gz"
		->"nicola.ibsMat"
		->"nicola.covMat"
		->"nicola.bcf"

NSITES=`zcat nicola.beagle.gz | wc -l`

echo $NSITES
#48,194 sites 

# scp *Mat, *qopt and bams files to laptop to look at hclust, PCA etc
#now do hclust in R and check for clones etc ibn R file "nicola_angsd_ibs"
#found two sets of mixups and two sets of clones (one in fav and one in franks)

###Now determine which symbiont they have
cat host.fasta symABCD.fast >holobiont.fasta
module load bowtie2
bowtie2-build holobiont.fasta holobiont.fasta

#map reads to holobiont
2bRAD_bowtie2_launch.pl '\.trim$' holobiont.fasta > bt2-b
scc6_qsub_launcher.py -N bt2 -P coral -M daviessw@gmail.com -j y -h_rt 24:00:00 -jobsfile bt2-b

qsub bt2_array.qsub

zooxtype.pl host ="XM_020744737.1" >zoox_counts.txt
#analyze the sym data in R file "Nicola_sym.R"


###now we know who the clones are and which symbionts they host. Given that all but KF host D we will map to Ofav+D reference
#------------------------------TRIMMING FILES
module load fastx-toolkit

# creating and launching the cleaning process for all files in the same time:
tagseq_trim_launch.pl '\.fastq$' > clean
##look at what is being done in clean

nano clean
#in the next steps you'll see that these commands are clipping poly A tails, Illumina adapters, short reads (<20bp) and low quality (to keep 90% of a read must be quality score >20)

# now execute all commands written to file 'clean', preferably in array format
scc6_qsub_launcher.py -N trim -P coral -jobsfile clean
#this should create a trim.array.qsub and a trim_array_commands.txt files

qsub trim_array.qsub

#let's see what is happening on the cluster
qstat -u daviessw (change to your username)

##when the job is done, have a look in the trim.e* file
cat trim.(tab complete)
##this has all of the info for trimming. You'll see many sequences are PCR duplicates because this is TagSeq data and remember that we incorporated the degenerate bases into the cDNA synthesis. 

#---------------------------------------Getting the reference transcriptomes
wget "https://zenodo.org/records/10151798/files/Orbicella_faveolata_gen_17.mrna-transcripts.fa?download=1"
#this gets the Ofav transcriptome from this paper: https://zenodo.org/records/10151798

#then I ran a series of commands to make a dummy seq2iso file for the new reference
grep '>' Orbicella_faveolata_gen_17.mrna-transcripts.fa | awk 'BEGIN{OFS="\t"}{print $1, $1}' >Orb_fav_seq2iso.tab
cat Orb_fav_seq2iso.tab | sed 's/>//g' >Orb_fav_seq2iso2.tab

##grabbed the newest genome from Dougan https://espace.library.uq.edu.au/view/UQ:27da3e7
gunzip Dtrenchii_CCMP2556_CDS_fasta.gz
grep '>' Dtrenchii_CCMP2556_CDS_fasta | wc -l
#76789
grep '>' Dtrenchii_CCMP2556_CDS_fasta | awk 'BEGIN{OFS="\t"}{print $1, $1}' >Dtrenchii_seq2iso.tab
cat Dtrenchii_seq2iso.tab | sed 's/>//g' >Dtrenchii_seq2iso2.tab

#Dtrenchi transcriptome was collapsing isoforms wrong, so ran this command to properly collapse reads by isogroups
awk '{split($2, a, "_i"); print $1, a[1]}' Dtrenchii_seq2iso2.tab >Dtrenchii_seq2iso3.tab

#then I concatenated the Ofav reference with the D reference
cat Orbicella_faveolata_gen_17.mrna-transcripts.fa Dtrenchii_CCMP2556_CDS_fasta >newholobiont.fasta

#same for the seq2iso files
cat Orb_fav_seq2iso.tab Dtrenchii_seq2iso.tab >newholobiont_seq2iso.tab

#---------------------------------------Making the mapping database for your reference transcriptome
module load bowtie2
# creating bowtie2 index for your transcriptome:
bowtie2-build newholobiont.fasta newholobiont.fasta 

#---------------------------------------Mapping reads to reference transcriptome

# cd where the trimmed read files are (extension "trim")
tagseq_bowtie2map.pl "trim$" holobiont.fasta  > maps
nano maps

scc6_qsub_launcher.py -N maps -P coral -jobsfile maps
#this should create a maps.array.qsub and a maps_array_commands.txt files

qsub maps_array.qsub


###now you have individual sam files for each trimmed files

# alignment rates:
nano maps.o(tab complete)

#---------------------------------------Generating read-counts-per gene 

# NOTE: Must have a tab-delimited file giving correspondence between contigs in the transcriptome fasta file and genes. Typically, each gene is represented by several contigs in the transcriptome. 
head newholobiont_seq2iso.tab

module load samtools
# counting hits per isogroup:
perl samcount_launch_bt2.pl '\.sam' holobiont_seq2iso.tab > sc
nano sc

scc6_qsub_launcher.py -N sc -P coral -jobsfile sc
#this should create a sc.array.qsub and a sc_array_commands.txt files

qsub sc_array.qsub

#nano sc.o(tab complete)
#you will see this: disregarding reads mapping to multiple isogroups
#we do not count reads that map to multiple places in this script, conservative approach.

#now you have individual counts files for each of your samples. Let's compile them into a single table!

# assembling all counts into a single table:
expression_compiler.pl *.sam.counts > nicola_holobiont_counts_new2.txt

head nicola_holobiont_counts_new2.txt

# DONE! Now split the counts bu host and symbiont and then analyze separately


