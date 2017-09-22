#!/bin/bash

# Placeholders for your job name, e-mail, etc
#SBATCH --job-name=star_gatk_variants
#SBATCH --mail-user=kelly.quek@sahmri.com
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --output=/homes/kquek/star/AA_star_gatk.log

# Resources allocation request parameters
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem-per-cpu=4096       # this will request 96 Gb of memory
#SBATCH --time=100:00:00          # Run time in hh:mm:ss


# This is job's working directory
echo "Starting at `date`"
echo "Running on hosts: $SLURM_NODELIST"
echo "Running on $SLURM_NNODES nodes."
echo "Running on $SLURM_NPROCS processors."
echo "Current working directory is `pwd`"


############ parameter to set #################################

#### software version used in this pipeline

## samtools v1.2
## htslib v1.2.1
## bedtools v2.21.0
## STAR v2.4.2a_modified
## java v1.8.0_51
## picard v1.136
## GATK v3.4.46
## ANNOVAR as of 28 Aug 2015
## SNPiR as of July 2013
## bcftools v1.2


############ Total estimated time for running a sample ########
### sample with ~89 million reads run ~28 hours full pipeline from step 1 to 8

### ANNOVAR results - user can use xxx.FINAL.txt for summary of output, or xxx.step10.gene/varlist output for variant prioritisation

## specify output parent and tmp directory
saveDir=/data/kquek/star_output                             ## path to save outputs - will be created once pipeline start
tmpy=/data/kquek/tmp                                        ## tmp directory for java - will be created once pipeline start

## fastq and sample name (id) information
id=AA                                                       ## sample id
read1=/homes/kquek/CNL/AA_161013_S1_combined_R1_001.fastq   ## assume is gz format, so do not need to put .gz at the back
read2=/homes/kquek/CNL/AA_161013_S1_combined_R2_001.fastq   ## assume is gz format, so do not need to put .gz at the back
rgsm=${id}                                                  ## RG header for RGSM, by default will be similiar to sample id

## path for references annotation for STAR
genomeDir=/homes/kquek/star/annotation/GRCh37_Gencode19     ## star genome directory 
GTF_ref=$genomeDir/gencode.v19.annotation.gtf               ## gtf file - gene annotation for star
ref_fasta=$genomeDir/hg19genome.fa                          ## assuming the reference genome fasta of .dict has been created by picard and samtools faidx (chr order ChrM,Chr1..X,Y)

## Software - path to software or scripts directory
software=/apps/packages/star/bin/Linux_x86_64
picardDir=/apps/packages/picard
GATK=/apps/packages/gatk
snpirDir=/homes/kquek/SNPiR
annovarDir=/homes/kquek/annovar

## assume intersectBed (from bedtools) is in the $PATH as SNPiR will use intersectBed (/usr/bin/intersectBed) for the analysis

## path for BaseRecalibration annotation downloaded from GATK bundle FTP
mills=/homes/kquek/gatk_resources/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf
dbsnp=/homes/kquek/gatk_resources/dbsnp_138.hg19.vcf
phase1indel=/homes/kquek/gatk_resources/1000G_phase1.indels.hg19.sites.vcf

## path for SNPiR annotation created according to SNPiR instructions
repeatRef=$snpirDir/hg19_repeat_masker_annotation.bed
snpirGeneRef=$snpirDir/hg19_overlap_candidates.sorted.txt
RNAeditRef=$snpirDir/Human_AG_all_hg19.bed


######################## parameter do not need to be set by user ############################################################

runDir=$saveDir/${id}_result_output/01_1pass                    
genomeDir2=$saveDir/${id}_result_output/02_hg19_2pass           
runDir2=$saveDir/${id}_result_output/03_2pass                  
output=$saveDir/${id}_result_output/04_MarkDuplicates
output1=$saveDir/${id}_result_output/05_BaseRecalibration
output2=$saveDir/${id}_result_output/06_VariantCalling
output3=$saveDir/${id}_result_output/07_SNPIR_Filtering
output4=$saveDir/${id}_result_output/08_ANNOVAR_Annotation
output5=$saveDir/${id}_result_output/08_ANNOVAR_Annotation/INDEL
output6=$saveDir/${id}_result_output/08_ANNOVAR_Annotation/SNP
output7=$saveDir/${id}_result_output/08_ANNOVAR_Annotation/SNP_noSNPiR

mkdir -pv $saveDir
mkdir -pv $tmpy
mkdir -pv $runDir
mkdir -pv $runDir2
mkdir -pv $genomeDir2
mkdir -pv $output
mkdir -pv $output1
mkdir -pv $output2
mkdir -pv $output3
mkdir -pv $output4
mkdir -pv $output5
mkdir -pv $output6
mkdir -pv $output7


##################### ANALYSIS START ###############################

#################### Unzip the FASTQ in gz format ##############################

gunzip ${read1}.gz
gunzip ${read2}.gz

#################### Execute STAR alignment ########################################
# 1-pass STAR 

cd $runDir

$software/STAR --genomeDir $genomeDir --readFilesIn $read1 $read2 --runThreadN 24 --sjdbGTFfile $GTF_ref --sjdbOverhang 100 --chimSegmentMin 15 --chimJunctionOverhangMin 15 --outReadsUnmapped Fastx

# 2-pass STAR

$software/STAR --runMode genomeGenerate --genomeDir $genomeDir2 --genomeFastaFiles $ref_fasta --sjdbFileChrStartEnd $runDir/SJ.out.tab --sjdbOverhang 100 --runThreadN 24

# Final alignment

cd $runDir2

$software/STAR --genomeDir $genomeDir2 --readFilesIn $read1 $read2 --runThreadN 24 --chimSegmentMin 15 --chimJunctionOverhangMin 15 --outReadsUnmapped Fastx

### alternate way for adding prefix and add RG
### --outFileNamePrefix ${id}. --outSAMattrRGline ID:${id} LB:Truseq PU:unkn PL:ILLUMINA SM:${id}
### out file name will be ${id}.Aligned.out.sam by star output

########## Rename the aligned reads

cd $runDir2
mv Aligned.out.sam ${id}.aligned.out.sam

################## Add read groups, sort, mark duplicates, and create index ########################

## Sort
java -Djava.io.tmpdir=$tmpy -jar $picardDir/picard.jar AddOrReplaceReadGroups I=$runDir2/${id}.aligned.out.sam O=$output/${id}_rg_added_sorted.bam SO=coordinate RGID=${id} RGLB=library RGPL=illumina RGPU=unit1 RGSM=${rgsm}

## Mark duplicates
java -Djava.io.tmpdir=$tmpy -jar $picardDir/picard.jar MarkDuplicates I=$output/${id}_rg_added_sorted.bam O=$output/${id}.dedupped.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=$output/${id}.output.metrics

################## Split'N'Trim and reassign mapping qualities ###################################

java -Djava.io.tmpdir=$tmpy -jar $GATK/GenomeAnalysisTK.jar -T SplitNCigarReads -R ${ref_fasta} -I $output/${id}.dedupped.bam -o $output/${id}.split.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS


################## Base recalibration ###################################################

java -Djava.io.tmpdir=$tmpy -jar $GATK/GenomeAnalysisTK.jar -T BaseRecalibrator -R ${ref_fasta} -I $output/${id}.split.bam -knownSites ${dbsnp} -knownSites ${mills} -knownSites ${phase1indel} -o $output1/recal_data.table

## second pass
java -Djava.io.tmpdir=$tmpy -jar $GATK/GenomeAnalysisTK.jar -T BaseRecalibrator -R ${ref_fasta} -I $output/${id}.split.bam -knownSites ${dbsnp} -knownSites ${mills} -knownSites ${phase1indel} -BQSR $output1/recal_data.table -o $output1/post_recal_data.table

## plot before/after plots
java -Djava.io.tmpdir=$tmpy -jar $GATK/GenomeAnalysisTK.jar -T AnalyzeCovariates -R ${ref_fasta} -before $output1/recal_data.table -after $output1/post_recal_data.table -plots $output1/recalibration_plots.pdf

## apply recalibration to the seq data
java -Djava.io.tmpdir=$tmpy -jar $GATK/GenomeAnalysisTK.jar -T PrintReads -R ${ref_fasta} -I $output/${id}.split.bam -BQSR $output1/recal_data.table -o $output1/${id}_recal_reads.bam

################# GATK HaplotypeCaller variant calling, variant tagged and filtering ##################################################

java -Djava.io.tmpdir=$tmpy -jar $GATK/GenomeAnalysisTK.jar -T HaplotypeCaller -R ${ref_fasta} -I $output1/${id}_recal_reads.bam -dontUseSoftClippedBases -stand_call_conf 20.0 -stand_emit_conf 20.0 -o $output2/${id}.HC.output.vcf

## variant filtering
java -Djava.io.tmpdir=$tmpy -jar $GATK/GenomeAnalysisTK.jar -T VariantFiltration -R ${ref_fasta} -V $output2/${id}.HC.output.vcf -window 35 -cluster 3 --filterName FS -filter "FS > 30.0" --filterName QD -filter "QD < 2.0" --filterName "LowCoverage" -filter "DP < 8 " -o $output2/${id}.HC.tagged.output.vcf

java -Djava.io.tmpdir=$tmpy -jar $GATK/GenomeAnalysisTK.jar -T SelectVariants -R ${ref_fasta} --variant $output2/${id}.HC.tagged.output.vcf -select 'vc.isNotFiltered()' -o $output2/${id}.HC.filtered.output.vcf


## perform QC and stats on filtered variant call

bgzip -c $output2/${id}.HC.filtered.output.vcf > $output2/${id}.HC.filtered.output.vcf.gz && tabix -p vcf $output2/${id}.HC.filtered.output.vcf.gz

bcftools stats $output2/${id}.HC.filtered.output.vcf.gz > $output2/${id}.HC.filtered.output.vcf.gz.check
 
plot-vcfstats -p $output2/${id}.HC.filtered.output.check.plots/ $output2/${id}.HC.filtered.output.vcf.gz.check

################## split SNP and INDEL out to separate VCF ##############################################

## INDEL VCF only
java -Djava.io.tmpdir=$tmpy -jar $GATK/GenomeAnalysisTK.jar -T SelectVariants -R ${ref_fasta} -V $output2/${id}.HC.filtered.output.vcf -o $output2/${id}.HC.filtered.INDEL.vcf -selectType INDEL

## SNP only
java -Djava.io.tmpdir=$tmpy -jar $GATK/GenomeAnalysisTK.jar -T SelectVariants -R ${ref_fasta} -V $output2/${id}.HC.filtered.output.vcf -o $output2/${id}.HC.filtered.SNP.vcf -xlSelectType INDEL

################# SNPiR filtering ####################################################################

## have to be in snpir directory

cd $snpirDir

## Convert VCF format to our custom variant format and filter variants with low quality

bash $snpirDir/convertVCF.sh $output2/${id}.HC.filtered.SNP.vcf $output3/${id}.HC.filtered.SNP.converted.txt 20

## Remove mismatches in first 6 bp of reads

perl $snpirDir/filter_mismatch_first6bp.pl -infile $output3/${id}.HC.filtered.SNP.converted.txt -outfile $output3/${id}.rmhex.txt -bamfile $output1/${id}_recal_reads.bam 

## filter variants in repetitive regions

awk '{OFS="\t";$2=$2-1"\t"$2;print $0}' $output3/${id}.rmhex.txt | intersectBed -a stdin -b $repeatRef -v | cut -f1,3-7 > $output3/${id}.rmhex.rmsk.txt

## filter intronic sites that are within 4bp of splicing junctions

perl $snpirDir/filter_intron_near_splicejuncts.pl -infile $output3/${id}.rmhex.rmsk.txt -outfile $output3/${id}.rmhex.rmsk.rmintron.txt -genefile $snpirGeneRef

## filter variants in homopolymers

perl $snpirDir/filter_homopolymer_nucleotides.pl -infile $output3/${id}.rmhex.rmsk.rmintron.txt -outfile $output3/${id}.rmhex.rmsk.rmintron.rmhom.txt -refgenome ${ref_fasta}

## filter variants that were caused by mismapped reads

perl $snpirDir/BLAT_candidates.pl -infile $output3/${id}.rmhex.rmsk.rmintron.txt -outfile $output3/${id}.rmhex.rmsk.rmintron.rmhom.rmblat.txt -bamfile $output1/${id}_recal_reads.bam -refgenome ${ref_fasta}

## remove known RNA editing sites

awk '{OFS="\t";$2=$2-1"\t"$2;print $0}' $output3/${id}.rmhex.rmsk.rmintron.rmhom.rmblat.txt | intersectBed -a stdin -b $RNAeditRef -v > $output3/${id}.rmhex.rmsk.rmintron.rmhom.rmblat.rmedit.bed


## Output file format:
## 1.) chromosome
## 2.) position-1
## 3.) position
## 4.) coverage, number of non-reference nucleotides
## 5.) reference nucleotide
## 6.) non-reference nucleotide
## 7.) frequency of non-reference nucleotide

######################### ANNOVAR ANNOTATION ##################################################

## change to annovar directory

cd $annovarDir

################################ INDEL ##############################################

perl $annovarDir/convert2annovar.pl --format vcf4 $output2/${id}.HC.filtered.INDEL.vcf --includeinfo --withfreq > $output5/${id}.HC.filtered.INDEL.avinput

## summarise variants

## INDEL
## if want to remove intermediate files, use --remove

perl $annovarDir/table_annovar.pl $output5/${id}.HC.filtered.INDEL.avinput humandb/ --buildver hg19 --out $output5/${id}.INDEL.anno --otherinfo --protocol refGene,knownGene,ensGene,wgEncodeGencodeBasicV19,genomicSuperDups,esp6500siv2_all,exac03,1000g2015aug_all,1000g2015aug_afr,1000g2015aug_amr,1000g2015aug_eas,1000g2015aug_eur,1000g2015aug_sas,cg69,snp138NonFlagged,avsnp138,snp138,ljb26_all,clinvar_20150629,nci60,cosmic70 -operation g,g,g,g,r,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f -nastring .

head -1 $output5/${id}.INDEL.anno.hg19_multianno.txt | sed -e 's/^/SampleID\t/' - > $output5/header.INDEL.txt
sed '1d' $output5/${id}.INDEL.anno.hg19_multianno.txt | sed -e "s/^/${id}$(printf '\t')/" - | cat $output5/header.INDEL.txt - > $output5/${id}.INDEL.anno.hg19_multianno.FINAL.txt

perl $annovarDir/variants_reduction.pl $output5/${id}.HC.filtered.INDEL.avinput humandb/ --genetype wgEncodeGencodeBasicV19 --remove --protocol nonsyn_splicing,genomicSuperDups,exac03,1000g2015aug_all,esp6500siv2_all,avsnp138,snp138,snp138NonFlagged,cg69,dominant --operation g,rr,f,f,f,f,f,f,f,m --outfile $output5/${id}.INDEL.reduce --buildver hg19


################################# SNP with SNPiR filtering ################################################

## convert to annovar format
awk '{print $1,$3,$3,$5,$6,$4"|"$7}' $output3/${id}.rmhex.rmsk.rmintron.rmhom.rmblat.rmedit.bed > $output6/${id}.SNP.avinput

## summarise variants

perl $annovarDir/table_annovar.pl $output6/${id}.SNP.avinput humandb/ --buildver hg19 -out $output6/${id}.SNP.anno --otherinfo --protocol refGene,knownGene,ensGene,wgEncodeGencodeBasicV19,genomicSuperDups,esp6500siv2_all,exac03,1000g2015aug_all,1000g2015aug_afr,1000g2015aug_amr,1000g2015aug_eas,1000g2015aug_eur,1000g2015aug_sas,cg69,snp138NonFlagged,avsnp138,snp138,ljb26_all,clinvar_20150629,nci60,cosmic70 -operation g,g,g,g,r,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f -nastring . 

## annovar output {out}.hg19_multianno.txt format
## add new 1st column for sample name

head -1 $output6/${id}.SNP.anno.hg19_multianno.txt | sed -e 's/^/SampleID\t/' - > $output6/header.SNP.txt
sed '1d' $output6/${id}.SNP.anno.hg19_multianno.txt | sed -e "s/^/${id}$(printf '\t')/" - | cat $output6/header.SNP.txt - > $output6/${id}.SNP.anno.hg19_multianno.FINAL.txt


## variant reduction to prioritise casual variants

perl $annovarDir/variants_reduction.pl $output6/${id}.SNP.avinput humandb/ --genetype wgEncodeGencodeBasicV19 --protocol nonsyn_splicing,genomicSuperDups,exac03,1000g2015aug_all,esp6500siv2_all,avsnp138,snp138,snp138NonFlagged,cg69,dominant --operation g,rr,f,f,f,f,f,f,f,m --outfile $output6/${id}.SNP.reduce --buildver hg19


################################# SNP without SNPiR ###########################################

## convert to annovar input

perl $annovarDir/convert2annovar.pl --format vcf4 $output2/${id}.HC.filtered.SNP.vcf --includeinfo --withfreq > $output7/${id}.HC.filtered.SNP.noSNPiR.avinput

## summarised variants

perl $annovarDir/table_annovar.pl $output7/${id}.HC.filtered.SNP.noSNPiR.avinput humandb/ --buildver hg19 -out $output7/${id}.SNP.noSNPiR.anno --otherinfo --protocol refGene,knownGene,ensGene,wgEncodeGencodeBasicV19,genomicSuperDups,esp6500siv2_all,exac03,1000g2015aug_all,1000g2015aug_afr,1000g2015aug_amr,1000g2015aug_eas,1000g2015aug_eur,1000g2015aug_sas,cg69,snp138NonFlagged,avsnp138,snp138,ljb26_all,clinvar_20150629,nci60,cosmic70 -operation g,g,g,g,r,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f -nastring .

head -1 $output7/${id}.SNP.noSNPiR.anno.hg19_multianno.txt | sed -e 's/^/SampleID\t/' - > $output7/header.SNP.noSNPiR.txt
sed '1d' $output7/${id}.SNP.noSNPiR.anno.hg19_multianno.txt | sed -e "s/^/${id}$(printf '\t')/" - | cat $output7/header.SNP.noSNPiR.txt - > $output7/${id}.SNP.noSNPiR.anno.hg19_multianno.FINAL.txt

## variant reduction to prioritise casual variants

perl $annovarDir/variants_reduction.pl $output7/${id}.HC.filtered.SNP.noSNPiR.avinput humandb/ --genetype wgEncodeGencodeBasicV19 --remove --protocol nonsyn_splicing,genomicSuperDups,exac03,1000g2015aug_all,esp6500siv2_all,avsnp138,snp138,snp138NonFlagged,cg69,dominant --operation g,rr,f,f,f,f,f,f,f,m --outfile $output7/${id}.SNP.noSNPiR.reduce --buildver hg19


################################ compressed back the fastq files #############################

gzip ${read1}
gzip ${read2}

################################ END ####################################


