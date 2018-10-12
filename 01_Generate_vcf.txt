###################################################################################################################################
#
# Copyright 2017 IRD
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/> or
# write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston,
# MA 02110-1301, USA.
#
# You should have received a copy of the CeCILL-C license with this program.
# If not see <http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.txt>
#
# Intellectual property belongs to IRD 
#
# Written by Nora Scarcelli
#
###################################################################################################################################


###################################################################################################################################
# This script generate a vcf file from raw .fastq data
# Needs files adapteurs.conf, Filter.pl and Compare.pl provided in folder 01_Generate_vcf
####################################################################################################################################


##############
# File preparation
##############

#Demultiplexing
for i in *R1.fastq; do python demultadapt.py -f $i -p Demul-R1 FILE_ADAPT; done

#Remove adapter
#Need file adapteurs.conf
for i in Demul-*; do cutadapt $(<adapteurs.conf) -o Cutadapt-$i $i; done
for i in *R2.fastq; do cutadapt $(<adapteurs.conf) -o Cutadapt-$i $i; done

#Quality filter
#Need file Filter.pl
for i in Cutadapt*; do perl Filter.pl -f $i -o Filter-$i; done

#R1 and R2 matching
#Need file Compare.pl
for i in Filter-Cutadapt-Demul*.fastq; do perl Compare.pl -f $i -r Filter-Cutadapt-R2.fastq -of Paired-$i -or R2-Paired-$i -os Single-$i; done



##############
# Mapping
##############

#Mapping
for i in Paired-*.fastq; do bwa mem ref.fasta -M -B 4 $i R2-$i> Sam-$i.sam; done

#Delete non properly paired reads
for i in *.sam; do samtools view -h -b -S -f 0x2 $i -o $i.bam; done

#Sort reads
for i in *.bam; do samtools sort $i Sort-$i; done



##############
# SNP calling
##############

#Index reference
samtools faidx ref.fasta
java -jar CreateSequenceDictionary.jar REFERENCE=ref.fasta OUTPUT= ref.dict

#Add read group
for i in *.bam; do java -jar AddOrReplaceReadGroups.jar INPUT= $i OUTPUT=RG-$i RGLB=WGS RGPL=Illumina RGPU=RunHi RGSM=$i VALIDATION_STRINGENCY=SILENT; done

#Index bam files
for i in RG*; do samtools index $i; done

#Local realignment for indel
for i in *.bam; do java -jar GenomeAnalysisTK.jar -T RealignerTargetCreator -R ref.fasta -I $i -o INDEL-$i.intervals; done
for i in *.bam; do java -jar GenomeAnalysisTK.jar -T IndelRealigner --targetIntervals INDEL-$i.intervals -o Realign-$i -I $i -R ref.fasta --maxReadsForRealignment 10000000; done

#Calling
java -jar GenomeAnalysisTK.jar -T UnifiedGenotyper -I bam.list -R ref.fasta -o Chr_Name_out.vcf --output_mode EMIT_VARIANTS_ONLY -glm SNP -stand_emit_conf 10 -mbq 30 -rfBadCigar -L Chr_name



##############
# SNP filtration
##############

for i in *.vcf; do java -jar GenomeAnalysisTK.jar -T VariantFiltration -R ref.fasta -V $i --filterExpression 'QUAL<200' --filterName 'LOW-QUAL' --filterExpression 'MQ0>=4 && ((MQ0/(1.0*DP))>0.1)' --filterName 'LOW-MQ0' --filterExpression 'DP<10' --filterName 'LOW-DP' --clusterSize 3 --clusterWindowSize 10 --filterExpression 'DP>20000' --filterName 'HIGH-DP' -o TAG_$i; vcftools --vcf TAG_$i --remove-filtered-all --recode --recode-INFO-all --out Filter_$i; done
for i in *.vcf; do vcftools --vcf $i --max-missing 0.75 --recode --out Missing_$i; vcftools --vcfMissing_$i --min-alleles 2 --max-alleles 2 --recode --out Final_$i; done

