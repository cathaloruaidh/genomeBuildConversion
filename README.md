# 1&nbsp; Genome Build Conversion
Code for identifying positions in the human genome that are unstable when converting between GRCh37 and GRCh38, using either `liftOver` or `CrossMap`. 
The data in mind for this procedure are single nucleotide variants (SNVs) obtained from whole genome sequencing (WGS). 
Previous work has highlighted unusual behaviour in build conversion, such as SNVs mapping to a different chromosome. 
BED files describing these conversion-unstable positions (CUPs) are provided, as well as a combined file of all regions ("novel_CUPs"). 
These positions are determined by the chain file, and are the same regardless of the conversion tool. 
Pre-excluding SNVs in these regions before converting between builds removes all unstable behaviour. 

If you have any queries or feedback, please contact the [author](mailto:cathalormond@gmail.com). If you use the exclude files or this algorithm in your publication, please cite the following paper:

Navigation: 
- [Set Up](#2-set-up) 

- [Full Genome Data](#3-full-genome-data)
- [Real WGS Example](#4-real-wgs-example)


# 2&nbsp; Set Up
## 2.1&nbsp; Notes
- Prerequisites: `liftOver`, `CrossMap`, `bedtools` and `GNU parallel`. Additionally, `bwa` is required for the WGS example. 
- Reference FASTA files are not included due to file size, but are required for the application to the real WGS data. Code to download and index these reference files is provided below. 
- This process assumes `chr1, chr2, ..., chrX, chrY` nomenclature. 
- The input BED files for the full-genome search for one build are ~150GB in size. Once the algorithm is applied, all files can take up to 1.5TB in size. 
- For a single run of the algorithm, a source and target such as GRCh37 and GRCh38 must be chosen as well as a tool (`liftOver` or `CrossMap`). 
- For the WGS example, JAR files for `picard` and `GATK` are included in the REFERENCE directory. 



## 2.2&nbsp; Installation
Download the resource material and initialise:
```
git clone https://github.com/cathaloruaidh/genomeBuildConversion.git
cd genomeBuildConversion ;

MAIN_DIR=${PWD}
REF_DIR=${MAIN_DIR}/REFERENCE
SCRIPT_DIR=${MAIN_DIR}/SCRIPTS
export REF_DIR
export SCRIPT_DIR
mkdir CHR COMBINE WGS_DATA ;
chmod +x ${SCRIPT_DIR}/*sh

```

Run for GRCh37: 
```
SOURCE=GRCh37
TARGET=GRCh38
export SOURCE TARGET
REGIONS=${REF_DIR}/GRCh37.region.Standard.bed

```

Run for GRCh38: 
```
SOURCE=GRCh38
TARGET=GRCh37
export SOURCE TARGET
REGIONS=${REF_DIR}/GRCh38.regions.Standard.bed 

```


Run for liftOver: 
```
TOOL=liftOver
export TOOL
LOOP_BED=${SCRIPT_DIR}/loop_${TOOL}.FullGenome.BED.sh

```

Run for CrossMap: 
```
TOOL=CrossMap
export TOOL
LOOP_BED=${SCRIPT_DIR}/loop_${TOOL}.FullGenome.BED.sh

```





# 3&nbsp; Full Genome Data
## 3.1&nbsp; Create Input BED

Generate the input BED files for the conversion process. 
Every individual base-pair position in the genome will have a BED entry, based on the lengths of the standard 23 pairs of chromosomes.
This is parallelised for speed, as it can take up to 90 minutes on the largest chromosome. 

```
cd CHR 
date ;
parallel --plus -j 12 --colsep '\t' "${SCRIPT_DIR}/createInputBed.sh {1} {2} {3} ${SOURCE}" :::: ${REGIONS}

```


## 3.2&nbsp; Apply Algorithm
Run the main script to identify unstable positions. 
The loop script takes as arguments the input filename, the start iteration, the end iteration and the source build. 
Both scripts will add the tool name to the file output, so there should be no over-writing of output files. 
This is parallelised using GNU parallel, with 12 CPUs available. 


```
parallel --plus -j12  ". ${LOOP_BED} {} 1 2 ${SOURCE}" ::: $( ls -d *_${SOURCE} | sort -V )

```

The script was set up so that iterations could be interrupted and restarted if neccessary. 
Two iterations were run above to determine if the algorithm was stable, or if new sites would be identified at each step (the former was expected and observed). 



## 3.3&nbsp; Sanity Check
Check if there are entries in the files of unstable positions for the second iteration by counting the lines (regardless of the source/target builds). 
The first two commands should return zeroes for all files, and the third command should return nothing. 

```
for FILE in $( find . -iname '*GRCh37_2.reject.extract.bed' ) ; do wc -l ${FILE} ; done
for FILE in $( find . -iname '*GRCh38_2.reject.extract.bed' ) ; do wc -l ${FILE} ; done
for FILE in $( find . -iname '*_2.jump*' ) ; do wc -l ${FILE} ; done

```


## 3.4&nbsp; Combine Sites
For each of the five CUP files, combine all the individual base-pair positions and collapse into multi-site regions. 
Additionally, combine the four novel CUPs into one file. 

```
cd ${MAIN_DIR}/COMBINE ; 

cat $( find ../CHR/ -iname "*${TARGET}_*jump_CHR.bed" ) | bedtools sort -i - | bedtools merge -i - > FASTA_BED.${TOOL}.ALL_${SOURCE}.CHR_jump_1.bed
cat $( find ../CHR/ -iname "*${SOURCE}_*jump_CHR.bed" ) | bedtools sort -i - | bedtools merge -i - > FASTA_BED.${TOOL}.ALL_${SOURCE}.CHR_jump_2.bed
cat $( find ../CHR/ -iname "*${SOURCE}_*jump_POS.bed" ) | bedtools sort -i - | bedtools merge -i - > FASTA_BED.${TOOL}.ALL_${SOURCE}.POS_jump.bed
cat $( find ../CHR/ -iname "*${SOURCE}_*.reject.extract.bed" ) | bedtools sort -i - | bedtools merge -i - > FASTA_BED.${TOOL}.ALL_${SOURCE}.reject_1.bed
cat $( find ../CHR/ -iname "*${TARGET}_*.reject.extract.bed" ) | bedtools sort -i - | bedtools merge -i - > FASTA_BED.${TOOL}.ALL_${SOURCE}.reject_2.bed
cat FASTA_BED.${TOOL}.ALL_${SOURCE}*jump*.bed FASTA_BED.${TOOL}.ALL_${SOURCE}*reject_2.bed | bedtools sort -i - | bedtools merge -i - > FASTA_BED.${TOOL}.ALL_${SOURCE}.novel_CUPs.bed

```



# 4&nbsp; Real WGS Example
As a proof of principle, a modified version of the above algorithm can be applied to WGS data. 
`liftOver` is implemented by the `LiftoverVCF` module from `picard`, since `liftOver` cannot directly process VCF files. 

## 4.1&nbsp; Download
Download VCF files for NA12877 and NA12878 from the [Illumina Platinum Genomes project](https://www.illumina.com/platinumgenomes.html): 

```
cd ${MAIN_DIR}/WGS_DATA
mkdir GRCh37 GRCh38

wget https://s3.eu-central-1.amazonaws.com/platinum-genomes/2017-1.0/hg19/small_variants/ConfidentRegions.bed.gz -O GRCh37/ConfidentRegions.GRCh37.bed.gz
wget https://s3.eu-central-1.amazonaws.com/platinum-genomes/2017-1.0/hg19/small_variants/ConfidentRegions.bed.gz.tbi -O GRCh37/ConfidentRegions.GRCh37.bed.gz.tbi
wget https://s3.eu-central-1.amazonaws.com/platinum-genomes/2017-1.0/hg19/small_variants/NA12877/NA12877.vcf.gz -O GRCh37/NA12877.GRCh37.vcf.gz
wget https://s3.eu-central-1.amazonaws.com/platinum-genomes/2017-1.0/hg19/small_variants/NA12877/NA12877.vcf.gz.tbi -O GRCh37/NA12877.GRCh37.vcf.gz.tbi
wget https://s3.eu-central-1.amazonaws.com/platinum-genomes/2017-1.0/hg19/small_variants/NA12878/NA12878.vcf.gz -O GRCh37/NA12878.GRCh37.vcf.gz
wget https://s3.eu-central-1.amazonaws.com/platinum-genomes/2017-1.0/hg19/small_variants/NA12878/NA12878.vcf.gz.tbi -O GRCh37/NA12878.GRCh37.vcf.gz.tbi

wget https://s3.eu-central-1.amazonaws.com/platinum-genomes/2017-1.0/hg38/small_variants/ConfidentRegions.bed.gz -O GRCh38/ConfidentRegions.GRCh38.bed.gz
wget https://s3.eu-central-1.amazonaws.com/platinum-genomes/2017-1.0/hg38/small_variants/ConfidentRegions.bed.gz.tbi -O GRCh38/ConfidentRegions.GRCh38.bed.gz.tbi
wget https://s3.eu-central-1.amazonaws.com/platinum-genomes/2017-1.0/hg38/small_variants/NA12877/NA12877.vcf.gz -O GRCh38/NA12877.GRCh38.vcf.gz
wget https://s3.eu-central-1.amazonaws.com/platinum-genomes/2017-1.0/hg38/small_variants/NA12877/NA12877.vcf.gz.tbi  -O GRCh38/NA12877.GRCh38.vcf.gz.tbi
wget https://s3.eu-central-1.amazonaws.com/platinum-genomes/2017-1.0/hg38/small_variants/NA12878/NA12878.vcf.gz -O GRCh38/NA12878.GRCh38.vcf.gz
wget https://s3.eu-central-1.amazonaws.com/platinum-genomes/2017-1.0/hg38/small_variants/NA12878/NA12878.vcf.gz.tbi -O GRCh38/NA12878.GRCh38.vcf.gz.tbi

```

Annotate the variants with unique identifier and extract bi-allelic SNVs subset to confident regions: 
```
for SAMPLE in NA12877 NA12878
do 
    for BUILD in GRCh37 GRCh38
    do
        echo "${SAMPLE} on ${BUILD}"
        bcftools query -f "%CHROM\t%POS\t%REF\t%ALT\n" ${BUILD}/${SAMPLE}.${BUILD}.vcf.gz | awk -v OFS="\t" -v SOURCE=${BUILD} '{print $1,$2-1,$2,$1 "_" $2 "_" $3 "_" $4 "_" SOURCE }' | gzip -c - > ${BUILD}/${SAMPLE}.${BUILD}.annotate.bed.gz
        bcftools annotate -c CHROM,FROM,TO,ID -a ${BUILD}/${SAMPLE}.${BUILD}.annotate.bed.gz -Oz -o ${BUILD}/${SAMPLE}.${BUILD}.annotate.vcf.gz ${BUILD}/${SAMPLE}.${BUILD}.vcf.gz
        tabix ${BUILD}/${SAMPLE}.${BUILD}.annotate.vcf.gz
        bcftools view -Ov --max-alleles 2 --types snps -R ${BUILD}/ConfidentRegions.${BUILD}.bed.gz ${BUILD}/${SAMPLE}.${BUILD}.annotate.vcf.gz > ${BUILD}/${SAMPLE}.${BUILD}.annotate.bi_SNV.original.vcf
        
        echo -e "\tExtracting SNVs at stable positions"
        vcftools --vcf ${BUILD}/${SAMPLE}.${BUILD}.annotate.bi_SNV.original.vcf --exclude-bed ${MAIN_DIR}/CUP_FILES/FASTA_BED.ALL_${BUILD}.novel_CUPs.bed --recode --recode-INFO-all --out ${BUILD}/${SAMPLE}.${BUILD}.annotate.bi_SNV 
        mv ${BUILD}/${SAMPLE}.${BUILD}.annotate.bi_SNV.recode.vcf ${BUILD}/${SAMPLE}.${BUILD}.annotate.bi_SNV.stable.vcf
        
        echo -e "\n\n"
    done
done

```

Create the directories: 
```
for BUILD in GRCh37 GRCh38
    do for CAT in original stable
        do for SAMPLE in NA12877 NA12878
            do for TYPE in VCF BED
                do mkdir -p ${BUILD}/${CAT}/${SAMPLE}/${TYPE}
            done
        done
    done
done

```

Finally, download reference genomes for GRCh37 and GRCh38, and index:
```
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/ucsc.hg19.fasta.gz -O ${REF_DIR}/GRCh37.fa.gz
gunzip ${REF_DIR}/GRCh37.fa.gz
samtools faidx ${REF_DIR}/GRCh37.fa
java -jar ${REF_DIR}/picard.jar CreateSequenceDictionary REFERENCE=${REF_DIR}/GRCh37.fa OUTPUT=${REF_DIR}/GRCh37.dict

wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa -O ${REF_DIR}/GRCh38.fa
samtools faidx ${REF_DIR}/GRCh38.fa
java -jar ${REF_DIR}/picard.jar CreateSequenceDictionary REFERENCE=${REF_DIR}/GRCh38.fa OUTPUT=${REF_DIR}/GRCh38.dict

```

Select the original or the stable data (i.e. SNVs at CUPs pre-excluded):
```
CATEGORY=original
# or
CATEGORY=stable

```

Additionally, select the WGS sample: 
```
SAMPLE=NA12877
# or
SAMPLE=NA12878

```



## 4.2&nbsp; VCF Data
Apply the algorithm: 
```
. ${SCRIPT_DIR}/loop_${TOOL}.WGS.VCF.sh ${SOURCE}/${SAMPLE}.${SOURCE}.annotate.bi_SNV.${CATEGORY}.vcf 1 2 ${SOURCE} ${SOURCE}/${CATEGORY}/${SAMPLE}/VCF

```

Check that there are no entries in the files of CUPs for the second iteration. 
All commands should return zeroes for the files. 
```
for FILE in $( find ${SOURCE}/${CATEGORY}/${SAMPLE}/VCF -iname "*${SOURCE}_2.reject.extract.bed" ) ; do wc -l ${FILE} ; done
for FILE in $( find ${SOURCE}/${CATEGORY}/${SAMPLE}/VCF -iname "*${TARGET}_2.reject.extract.bed" ) ; do wc -l ${FILE} ; done
for FILE in $( find ${SOURCE}/${CATEGORY}/${SAMPLE}/VCF -iname "*_2.*jump*" ) ; do wc -l ${FILE} ; done

```

Combine the results. 
```
FILES=$( find ${SOURCE}/${CATEGORY}/${SAMPLE}/VCF/ -iname "${SAMPLE}*_${TARGET}_*CHR_jump.bed" ) ; if [[ ${#FILES[@]} > 0 ]] ; then cat ${FILES[@]} | bedtools sort -i - | bedtools merge -i - > ${SOURCE}/${CATEGORY}/${SAMPLE}/${SAMPLE}.${TOOL}.${SOURCE}.${CATEGORY}.VCF.CHR_jump_1.bed ; else touch ${SOURCE}/${CATEGORY}/${SAMPLE}/${SAMPLE}.${TOOL}.${SOURCE}.${CATEGORY}.VCF.CHR_jump_1.bed ; fi

FILES=$( find ${SOURCE}/${CATEGORY}/${SAMPLE}/VCF/ -iname "${SAMPLE}*_${SOURCE}_*CHR_jump.bed" ) ; if [[ ${#FILES[@]} > 0 ]] ; then cat ${FILES[@]} | bedtools sort -i - | bedtools merge -i - > ${SOURCE}/${CATEGORY}/${SAMPLE}/${SAMPLE}.${TOOL}.${SOURCE}.${CATEGORY}.VCF.CHR_jump_2.bed ; else touch ${SOURCE}/${CATEGORY}/${SAMPLE}/${SAMPLE}.${TOOL}.${SOURCE}.${CATEGORY}.VCF.CHR_jump_2.bed ; fi

FILES=$( find ${SOURCE}/${CATEGORY}/${SAMPLE}/VCF/ -iname "${SAMPLE}*_${SOURCE}_*POS_jump.bed" ) ; if [[ ${#FILES[@]} > 0 ]] ; then cat ${FILES[@]} | bedtools sort -i - | bedtools merge -i - > ${SOURCE}/${CATEGORY}/${SAMPLE}/${SAMPLE}.${TOOL}.${SOURCE}.${CATEGORY}.VCF.POS_jump.bed ; else touch ${SOURCE}/${CATEGORY}/${SAMPLE}/${SAMPLE}.${TOOL}.${SOURCE}.${CATEGORY}.VCF.POS_jump.bed ; fi

FILES=$( find ${SOURCE}/${CATEGORY}/${SAMPLE}/VCF/ -iname "${SAMPLE}*_${SOURCE}_*.reject.extract.bed" ) ; if [[ ${#FILES[@]} > 0 ]] ; then cat ${FILES[@]} | bedtools sort -i - | bedtools merge -i - > ${SOURCE}/${CATEGORY}/${SAMPLE}/${SAMPLE}.${TOOL}.${SOURCE}.${CATEGORY}.VCF.reject_1.bed ; else touch ${SOURCE}/${CATEGORY}/${SAMPLE}/${SAMPLE}.${TOOL}.${SOURCE}.${CATEGORY}.VCF.reject_1.bed ; fi

FILES=$( find ${SOURCE}/${CATEGORY}/${SAMPLE}/VCF/ -iname "${SAMPLE}*_${TARGET}_*.reject.extract.bed" ) ; if [[ ${#FILES[@]} > 0 ]] ; then cat ${FILES[@]} | bedtools sort -i - | bedtools merge -i - > ${SOURCE}/${CATEGORY}/${SAMPLE}/${SAMPLE}.${TOOL}.${SOURCE}.${CATEGORY}.VCF.reject_2.bed ; else touch ${SOURCE}/${CATEGORY}/${SAMPLE}/${SAMPLE}.${TOOL}.${SOURCE}.${CATEGORY}.VCF.reject_2.bed ; fi

FILES=$( find ${SOURCE}/${CATEGORY}/${SAMPLE}/VCF/ -iname "${SAMPLE}*_${SOURCE}_*.mismatch.bed" ) ; if [[ ${#FILES[@]} > 0 ]] ; then cat ${FILES[@]} | bedtools sort -i - | bedtools merge -i - > ${SOURCE}/${CATEGORY}/${SAMPLE}/${SAMPLE}.${TOOL}.${SOURCE}.${CATEGORY}.VCF.mismatch_1.bed ; else touch ${SOURCE}/${CATEGORY}/${SAMPLE}/${SAMPLE}.${TOOL}.${SOURCE}.${CATEGORY}.VCF.mismatch_1.bed ; fi

FILES=$( find ${SOURCE}/${CATEGORY}/${SAMPLE}/VCF/ -iname "${SAMPLE}*_${TARGET}_*.mismatch.bed" ) ; if [[ ${#FILES[@]} > 0 ]] ; then cat ${FILES[@]} | bedtools sort -i - | bedtools merge -i - > ${SOURCE}/${CATEGORY}/${SAMPLE}/${SAMPLE}.${TOOL}.${SOURCE}.${CATEGORY}.VCF.mismatch_2.bed ; else touch ${SOURCE}/${CATEGORY}/${SAMPLE}/${SAMPLE}.${TOOL}.${SOURCE}.${CATEGORY}.VCF.mismatch_2.bed ; fi 

```


## 4.3&nbsp; BED Data
Convert VCF data to BED data
```
bcftools query -f "%CHROM\t%POS\t%REF\t%ALT\n" ${SOURCE}/${SAMPLE}.${SOURCE}.annotate.bi_SNV.${CATEGORY}.vcf.gz | awk -v OFS="\t" -v REF=${SOURCE} '{print $1,$2-1,$2,$1 "_" $2-1 "_" $3 "_" $4 "_" REF }'  > ${SOURCE}/${SAMPLE}.${SOURCE}.annotate.bi_SNV.${CATEGORY}.bed

```

Apply the algorithm: 
```
. ${SCRIPT_DIR}/loop_${TOOL}.WGS.BED.sh ${SOURCE}/${SAMPLE}.${SOURCE}.annotate.bi_SNV.${CATEGORY}.bed 1 2 ${SOURCE} ${SOURCE}/${CATEGORY}/${SAMPLE}/BED

```

Check that there are no entries in the files of CUPs for the second iteration. 
All commands should return zeroes for the files. 
```
for FILE in $( find ${SOURCE}/${CATEGORY}/${SAMPLE}/BED -iname "*${SOURCE}_2.reject.extract.bed" ) ; do wc -l ${FILE} ; done
for FILE in $( find ${SOURCE}/${CATEGORY}/${SAMPLE}/BED -iname "*${TARGET}_2.reject.extract.bed" ) ; do wc -l ${FILE} ; done
for FILE in $( find ${SOURCE}/${CATEGORY}/${SAMPLE}/BED -iname "*_2.*jump*" ) ; do wc -l ${FILE} ; done

```

Combine the results. 
```
FILES=$( find ${SOURCE}/${CATEGORY}/${SAMPLE}/BED/ -iname "${SAMPLE}*_${TARGET}_*CHR_jump.bed" ) ; if [[ ${#FILES[@]} > 0 ]] ; then cat ${FILES[@]} | bedtools sort -i - | bedtools merge -i - > ${SOURCE}/${CATEGORY}/${SAMPLE}/${SAMPLE}.${TOOL}.${SOURCE}.${CATEGORY}.BED.CHR_jump_1.bed ; else touch ${SOURCE}/${CATEGORY}/${SAMPLE}/${SAMPLE}.${TOOL}.${SOURCE}.${CATEGORY}.BED.CHR_jump_1.bed ; fi

FILES=$( find ${SOURCE}/${CATEGORY}/${SAMPLE}/BED/ -iname "${SAMPLE}*_${SOURCE}_*CHR_jump.bed" ) ; if [[ ${#FILES[@]} > 0 ]] ; then cat ${FILES[@]} | bedtools sort -i - | bedtools merge -i - > ${SOURCE}/${CATEGORY}/${SAMPLE}/${SAMPLE}.${TOOL}.${SOURCE}.${CATEGORY}.BED.CHR_jump_2.bed ; else touch ${SOURCE}/${CATEGORY}/${SAMPLE}/${SAMPLE}.${TOOL}.${SOURCE}.${CATEGORY}.BED.CHR_jump_2.bed ; fi

FILES=$( find ${SOURCE}/${CATEGORY}/${SAMPLE}/BED/ -iname "${SAMPLE}*_${SOURCE}_*POS_jump.bed" ) ; if [[ ${#FILES[@]} > 0 ]] ; then cat ${FILES[@]} | bedtools sort -i - | bedtools merge -i - > ${SOURCE}/${CATEGORY}/${SAMPLE}/${SAMPLE}.${TOOL}.${SOURCE}.${CATEGORY}.BED.POS_jump.bed ; else touch ${SOURCE}/${CATEGORY}/${SAMPLE}/${SAMPLE}.${TOOL}.${SOURCE}.${CATEGORY}.BED.POS_jump.bed ; fi

FILES=$( find ${SOURCE}/${CATEGORY}/${SAMPLE}/BED/ -iname "${SAMPLE}*_${SOURCE}_*.reject.extract.bed" ) ; if [[ ${#FILES[@]} > 0 ]] ; then cat ${FILES[@]} | bedtools sort -i - | bedtools merge -i - > ${SOURCE}/${CATEGORY}/${SAMPLE}/${SAMPLE}.${TOOL}.${SOURCE}.${CATEGORY}.BED.reject_1.bed ; else touch ${SOURCE}/${CATEGORY}/${SAMPLE}/${SAMPLE}.${TOOL}.${SOURCE}.${CATEGORY}.BED.reject_1.bed ; fi

FILES=$( find ${SOURCE}/${CATEGORY}/${SAMPLE}/BED/ -iname "${SAMPLE}*_${TARGET}_*.reject.extract.bed" ) ; if [[ ${#FILES[@]} > 0 ]] ; then cat ${FILES[@]} | bedtools sort -i - | bedtools merge -i - > ${SOURCE}/${CATEGORY}/${SAMPLE}/${SAMPLE}.${TOOL}.${SOURCE}.${CATEGORY}.BED.reject_2.bed ; else touch ${SOURCE}/${CATEGORY}/${SAMPLE}/${SAMPLE}.${TOOL}.${SOURCE}.${CATEGORY}.BED.reject_2.bed ; fi

```


## 4.4&nbsp; Compare Data
Get the jaccard index between the BED and VCF data. 
```
for BED in ${SOURCE}/${CATEGORY}/${SAMPLE}/*BED*.bed ; do VCF=$( echo ${BED} | sed -e 's/BED/VCF/g' ) ; echo -e "${BED}\t$(bedtools jaccard -a ${BED} -b ${VCF} | cut -f3 | tail -1 )" ; done | column -t 

```






