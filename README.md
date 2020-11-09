# 1&nbsp; Genome Build Conversion
Code for identifying positions in the human genome that are unstable when converting between GRCh37 and GRCh38, using either `liftOver` or `CrossMap`. 
The data in mind for this procedure are single nucleotide variants (SNVs) obtained from whole genome sequencing (WGS). 
Previous work has highlighted unusual behaviour in build conversion, such as SNVs mapping to a different chromosome. 
BED files describing these conversion-unstable positions (CUPs) are provided, as well as a combined file of all regions ("novel_CUPs"). 
These positions are determined by the chain file, and are the same regardless of the conversion tool. 
Pre-excluding SNVs in these regions before converting between builds removes all unstable behaviour. 
BED files for these positions may be downloaded in Section 2 below. 

If you have any queries or feedback, please contact the [author](mailto:cathalormond@gmail.com). 
<!---
If you use the exclude files or this algorithm in your publication, please cite the following paper:
--->

Navigation: 
- [CUP Bed Files](#2-novel-cup-bed-files)
- [Algorithm](#3-algorithm) 
- [Full Genome Data](#4-full-genome-data)
- [WGS Example](#5-wgs-example)

Note that sections 3-5 below need only be applied if a user wishes to re-generate the CUP files for validation, or to apply this process to genomes other than GRCh37 and GRCh38. 


# 2&nbsp; Novel CUP Bed Files
BED file containing novel CUP positions for either GRCh37 or GRCh38 can be downloaded directly here: 
```
# GRCh37
wget https://raw.githubusercontent.com/cathaloruaidh/genomeBuildConversion/master/CUP_FILES/FASTA_BED.ALL_GRCh37.novel_CUPs.bed

```

```
# GRCh38
wget https://raw.githubusercontent.com/cathaloruaidh/genomeBuildConversion/master/CUP_FILES/FASTA_BED.ALL_GRCh38.novel_CUPs.bed

```


SNVs at these positions can be removed from a VCF file using, for example, `vcftools` (substitute the appropriate VCF file name and BED genome source build):
```
vcftools --vcf INPUT.vcf --exclude-bed FASTA_BED.ALL_GRCh3N.novel_CUPs.bed --recode --recode-INFO-all --out INPUT.stable

```

The variants in the resulting VCF file are now stable to genome build conversion. 
Note that some variants may still fail the conversion process, due to their positions not being present in the target build. 
This is expected behaviour and is handled by the conversion tools. 


# 3&nbsp; Algorithm
## 3.1&nbsp; Notes
- Prerequisites: 
    * `liftOver`
    * `CrossMap`
    * `bedtools`
    * `bcftools`
    * `vcftools` 
    * `GNU parallel` (optional)
- Additionally, the following are required for the WGS example
    * `bwa` 
    * `picard` (supplied in REFERENCE directory)
    * `GATK` (supplied in REFERENCE directory)
- Reference FASTA files are not included due to file size, but are required for the application to the real WGS data. Code to download and index these reference files is provided below. 
- This process assumes `chr1, chr2, ..., chrX, chrY` nomenclature. 
- The input BED files for the full-genome search for one build are ~150GB in size. Once the algorithm is applied, all files can take up to 1.5TB in size. 
- For a single run of the algorithm, a source and target such as GRCh37 and GRCh38 must be chosen as well as a tool (`liftOver` or `CrossMap`). 


## 3.2&nbsp; Installation
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





# 4&nbsp; Full Genome Data
## 4.1&nbsp; Create Input BED

Generate the input BED files for the conversion process. 
Every individual base-pair position in the genome will have a BED entry, based on the lengths of the standard 23 pairs of chromosomes.
This is parallelised for speed, as it can take up to 90 minutes on the largest chromosome. 

```
cd CHR 
date ;
parallel --plus -j 12 --colsep '\t' "${SCRIPT_DIR}/createInputBed.sh {1} {2} {3} ${SOURCE}" :::: ${REGIONS}

```
If you are working on a different genome than GRCh37 or GRCh38, you can supply alternate REGION files in the REFERENCE directory. 


## 4.2&nbsp; Apply Algorithm
Run the main script to identify unstable positions. 
The loop script takes as arguments the input filename, the start iteration, the end iteration and the source build. 
Both scripts will add the tool name to the file output, so there should be no over-writing of output files. 
This is parallelised using GNU parallel, with 12 CPUs available. 


```
parallel --plus -j12  ". ${LOOP_BED} {} 1 2 ${SOURCE}" ::: $( ls -d *_${SOURCE} | sort -V )

```

The script was set up so that iterations could be interrupted and restarted if neccessary. 
Two iterations were run above to determine if the algorithm was stable, or if new sites would be identified at each step (the former was expected and observed). 



## 4.3&nbsp; Sanity Check
Check if there are entries in the files of unstable positions for the second iteration by counting the lines (regardless of the source/target builds). 
The first two commands should return zeroes for all files, and the third command should return nothing. 

```
for FILE in $( find *_${SOURCE} -iname "*${TOOL}_GRCh37_2.reject.extract.bed" ) ; do wc -l ${FILE} ; done
for FILE in $( find *_${SOURCE} -iname "*${TOOL}_GRCh38_2.reject.extract.bed" ) ; do wc -l ${FILE} ; done
for FILE in $( find *_${SOURCE} -iname '*${TOOL}*_2.jump*' ) ; do wc -l ${FILE} ; done

```


## 4.4&nbsp; Combine Sites
For each of the five CUP files, combine all the individual base-pair positions and collapse into multi-site regions. 
Additionally, combine the four novel CUPs into one file. 

```
cd ${MAIN_DIR}/COMBINE ; 

cat $( find ../CHR/*_${SOURCE} -iname "FASTA_BED.chr*_${SOURCE}.${TOOL}_${TARGET}_*.jump_CHR.bed" ) | bedtools sort -i - | bedtools merge -i - > FASTA_BED.${TOOL}.ALL_${SOURCE}.CHR_jump_1.bed
cat $( find ../CHR/*_${SOURCE} -iname "FASTA_BED.chr*_${SOURCE}.${TOOL}_${SOURCE}_*.jump_CHR.bed" ) | bedtools sort -i - | bedtools merge -i - > FASTA_BED.${TOOL}.ALL_${SOURCE}.CHR_jump_2.bed
cat $( find ../CHR/*_${SOURCE} -iname "FASTA_BED.chr*_${SOURCE}.${TOOL}_${SOURCE}_*.jump_POS.bed" ) | bedtools sort -i - | bedtools merge -i - > FASTA_BED.${TOOL}.ALL_${SOURCE}.POS_jump.bed
cat $( find ../CHR/*_${SOURCE} -iname "FASTA_BED.chr*_${SOURCE}.${TOOL}_${SOURCE}_*.reject.extract.bed" ) | bedtools sort -i - | bedtools merge -i - > FASTA_BED.${TOOL}.ALL_${SOURCE}.reject_1.bed
cat $( find ../CHR/*_${SOURCE} -iname "FASTA_BED.chr*_${SOURCE}.${TOOL}_${TARGET}_*.reject.extract.bed" ) | bedtools sort -i - | bedtools merge -i - > FASTA_BED.${TOOL}.ALL_${SOURCE}.reject_2.bed
cat FASTA_BED.${TOOL}.ALL_${SOURCE}*jump*.bed FASTA_BED.${TOOL}.ALL_${SOURCE}*reject_2.bed | bedtools sort -i - | bedtools merge -i - > FASTA_BED.${TOOL}.ALL_${SOURCE}.novel_CUPs.bed

```

## 4.5&nbsp; Compare 
To confirm that both `liftOver` and `CrossMap` give identical output for each of the CUP caregoties, we calculate the jaccard indices between the files:
```
cd ${MAIN_DIR}/COMBINE
for LIFT in FASTA_BED.liftOver.ALL_${SOURCE}.* ; do CROSS=$( echo ${LIFT} | sed -e 's/liftOver/CrossMap/g' ) ; echo -e "${LIFT}\t$( bedtools jaccard -a ${LIFT} -b ${CROSS} | cut -f3 | tail -1 )" ; echo  ; done | column -t 

```



# 5&nbsp; WGS Example
Note, the code in this section is currently experimental. 
As a proof of principle, a modified version of the above algorithm can be applied to WGS data. 
`liftOver` is implemented by the `LiftoverVCF` module from `picard`, since `liftOver` cannot directly process VCF files. 
The VCF information can be collapsed down to position information only as BED files, and the original algorithm run. 
The VCF-based data should be contained within the BED-based data, which in turn should be contained within the full-genome data. 

## 5.1&nbsp; Download
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



## 5.2&nbsp; VCF Data
Apply the algorithm: 
```
. ${SCRIPT_DIR}/loop_${TOOL}.WGS.VCF.sh ${SOURCE}/${SAMPLE}.${SOURCE}.annotate.bi_SNV.${CATEGORY}.vcf 1 2 ${SOURCE} ${SOURCE}/${CATEGORY}/${SAMPLE}/VCF

```

Check that there are no entries in the files of CUPs for the second iteration. 
All commands should return zeroes for the files. 
```
for FILE in $( find ${SOURCE}/${CATEGORY}/${SAMPLE}/VCF -iname "${SAMPLE}.${TOOL}.${SOURCE}.annotate.bi_SNV.${CATEGORY}_${SOURCE}_2.reject.extract.bed" ) ; do wc -l ${FILE} ; done
for FILE in $( find ${SOURCE}/${CATEGORY}/${SAMPLE}/VCF -iname "${SAMPLE}.${TOOL}.${SOURCE}.annotate.bi_SNV.${CATEGORY}_${TARGET}_2.reject.extract.bed" ) ; do wc -l ${FILE} ; done
for FILE in $( find ${SOURCE}/${CATEGORY}/${SAMPLE}/VCF -iname "${SAMPLE}.${TOOL}.${SOURCE}.annotate.bi_SNV.${CATEGORY}_*_2.*jump*" ) ; do wc -l ${FILE} ; done

```

Combine the results. 
```
FILES=( $( find ${SOURCE}/${CATEGORY}/${SAMPLE}/VCF/ -iname "${SAMPLE}.${SOURCE}*${TOOL}_${TARGET}_*jump_CHR.bed" ) ) ; if [[ ${#FILES[@]} > 0 ]] ; then cat ${FILES[@]} | bedtools sort -i - | bedtools merge -i - > ${SOURCE}/${CATEGORY}/${SAMPLE}/${SAMPLE}.${TOOL}.${SOURCE}.${CATEGORY}.VCF.CHR_jump_1.bed ; else touch ${SOURCE}/${CATEGORY}/${SAMPLE}/${SAMPLE}.${TOOL}.${SOURCE}.${CATEGORY}.VCF.CHR_jump_1.bed ; fi

FILES=( $( find ${SOURCE}/${CATEGORY}/${SAMPLE}/VCF/ -iname "${SAMPLE}.${SOURCE}*${TOOL}_${SOURCE}_*jump_CHR.bed" ) ) ; if [[ ${#FILES[@]} > 0 ]] ; then cat ${FILES[@]} | bedtools sort -i - | bedtools merge -i - > ${SOURCE}/${CATEGORY}/${SAMPLE}/${SAMPLE}.${TOOL}.${SOURCE}.${CATEGORY}.VCF.CHR_jump_2.bed ; else touch ${SOURCE}/${CATEGORY}/${SAMPLE}/${SAMPLE}.${TOOL}.${SOURCE}.${CATEGORY}.VCF.CHR_jump_2.bed ; fi

FILES=( $( find ${SOURCE}/${CATEGORY}/${SAMPLE}/VCF/ -iname "${SAMPLE}.${SOURCE}*${TOOL}_${SOURCE}_*jump_POS.bed" ) ) ; if [[ ${#FILES[@]} > 0 ]] ; then cat ${FILES[@]} | bedtools sort -i - | bedtools merge -i - > ${SOURCE}/${CATEGORY}/${SAMPLE}/${SAMPLE}.${TOOL}.${SOURCE}.${CATEGORY}.VCF.POS_jump.bed ; else touch ${SOURCE}/${CATEGORY}/${SAMPLE}/${SAMPLE}.${TOOL}.${SOURCE}.${CATEGORY}.VCF.POS_jump.bed ; fi

FILES=( $( find ${SOURCE}/${CATEGORY}/${SAMPLE}/VCF/ -iname "${SAMPLE}.${SOURCE}*${TOOL}_${SOURCE}_*.reject.extract.bed" ) ) ; if [[ ${#FILES[@]} > 0 ]] ; then cat ${FILES[@]} | bedtools sort -i - | bedtools merge -i - > ${SOURCE}/${CATEGORY}/${SAMPLE}/${SAMPLE}.${TOOL}.${SOURCE}.${CATEGORY}.VCF.reject_1.bed ; else touch ${SOURCE}/${CATEGORY}/${SAMPLE}/${SAMPLE}.${TOOL}.${SOURCE}.${CATEGORY}.VCF.reject_1.bed ; fi

FILES=( $( find ${SOURCE}/${CATEGORY}/${SAMPLE}/VCF/ -iname "${SAMPLE}.${SOURCE}*${TOOL}_${TARGET}_*.reject.extract.bed" ) ) ; if [[ ${#FILES[@]} > 0 ]] ; then cat ${FILES[@]} | bedtools sort -i - | bedtools merge -i - > ${SOURCE}/${CATEGORY}/${SAMPLE}/${SAMPLE}.${TOOL}.${SOURCE}.${CATEGORY}.VCF.reject_2.bed ; else touch ${SOURCE}/${CATEGORY}/${SAMPLE}/${SAMPLE}.${TOOL}.${SOURCE}.${CATEGORY}.VCF.reject_2.bed ; fi

FILES=( $( find ${SOURCE}/${CATEGORY}/${SAMPLE}/VCF/ -iname "${SAMPLE}.${SOURCE}*${TOOL}_${SOURCE}_*.mismatch.bed" ) ) ; if [[ ${#FILES[@]} > 0 ]] ; then cat ${FILES[@]} | bedtools sort -i - | bedtools merge -i - > ${SOURCE}/${CATEGORY}/${SAMPLE}/${SAMPLE}.${TOOL}.${SOURCE}.${CATEGORY}.VCF.mismatch_1.bed ; else touch ${SOURCE}/${CATEGORY}/${SAMPLE}/${SAMPLE}.${TOOL}.${SOURCE}.${CATEGORY}.VCF.mismatch_1.bed ; fi

FILES=( $( find ${SOURCE}/${CATEGORY}/${SAMPLE}/VCF/ -iname "${SAMPLE}.${SOURCE}*${TOOL}_${TARGET}_*.mismatch.bed" ) ) ; if [[ ${#FILES[@]} > 0 ]] ; then cat ${FILES[@]} | bedtools sort -i - | bedtools merge -i - > ${SOURCE}/${CATEGORY}/${SAMPLE}/${SAMPLE}.${TOOL}.${SOURCE}.${CATEGORY}.VCF.mismatch_2.bed ; else touch ${SOURCE}/${CATEGORY}/${SAMPLE}/${SAMPLE}.${TOOL}.${SOURCE}.${CATEGORY}.VCF.mismatch_2.bed ; fi 

```


## 5.3&nbsp; BED Data
Convert VCF data to BED data
```
bcftools query -f "%CHROM\t%POS\t%REF\t%ALT\n" ${SOURCE}/${SAMPLE}.${SOURCE}.annotate.bi_SNV.${CATEGORY}.vcf | awk -v OFS="\t" -v REF=${SOURCE} '{print $1,$2-1,$2,$1 "_" $2-1 "_" $3 "_" $4 "_" REF }'  > ${SOURCE}/${SAMPLE}.${SOURCE}.annotate.bi_SNV.${CATEGORY}.bed

```

Apply the algorithm: 
```
. ${SCRIPT_DIR}/loop_${TOOL}.WGS.BED.sh ${SOURCE}/${SAMPLE}.${SOURCE}.annotate.bi_SNV.${CATEGORY}.bed 1 2 ${SOURCE} ${SOURCE}/${CATEGORY}/${SAMPLE}/BED

```

Check that there are no entries in the files of CUPs for the second iteration. 
The first two commands should return zeroes for all files, and the third command should return nothing. 
```
for FILE in $( find ${SOURCE}/${CATEGORY}/${SAMPLE}/BED -iname "${SAMPLE}.${TOOL}.${SOURCE}.annotate.bi_SNV.${CATEGORY}_${SOURCE}_2.reject.extract.bed" ) ; do wc -l ${FILE} ; done
for FILE in $( find ${SOURCE}/${CATEGORY}/${SAMPLE}/BED -iname "${SAMPLE}.${TOOL}.${SOURCE}.annotate.bi_SNV.${CATEGORY}_${TARGET}_2.reject.extract.bed" ) ; do wc -l ${FILE} ; done
for FILE in $( find ${SOURCE}/${CATEGORY}/${SAMPLE}/BED -iname "${SAMPLE}.${TOOL}.${SOURCE}.annotate.bi_SNV.${CATEGORY}__2.*jump*" ) ; do wc -l ${FILE} ; done

```

Combine the results. 
```
FILES=( $( find ${SOURCE}/${CATEGORY}/${SAMPLE}/BED/ -iname "${SAMPLE}.${SOURCE}*${TOOL}_${TARGET}_*jump_CHR.bed" ) ) ; if [[ ${#FILES[@]} > 0 ]] ; then cat ${FILES[@]} | bedtools sort -i - | bedtools merge -i - > ${SOURCE}/${CATEGORY}/${SAMPLE}/${SAMPLE}.${TOOL}.${SOURCE}.${CATEGORY}.BED.CHR_jump_1.bed ; else touch ${SOURCE}/${CATEGORY}/${SAMPLE}/${SAMPLE}.${TOOL}.${SOURCE}.${CATEGORY}.BED.CHR_jump_1.bed ; fi

FILES=( $( find ${SOURCE}/${CATEGORY}/${SAMPLE}/BED/ -iname "${SAMPLE}.${SOURCE}*${TOOL}_${SOURCE}_*jump_CHR.bed" ) ) ; if [[ ${#FILES[@]} > 0 ]] ; then cat ${FILES[@]} | bedtools sort -i - | bedtools merge -i - > ${SOURCE}/${CATEGORY}/${SAMPLE}/${SAMPLE}.${TOOL}.${SOURCE}.${CATEGORY}.BED.CHR_jump_2.bed ; else touch ${SOURCE}/${CATEGORY}/${SAMPLE}/${SAMPLE}.${TOOL}.${SOURCE}.${CATEGORY}.BED.CHR_jump_2.bed ; fi

FILES=( $( find ${SOURCE}/${CATEGORY}/${SAMPLE}/BED/ -iname "${SAMPLE}.${SOURCE}*${TOOL}_${SOURCE}_*jump_POS.bed" ) ) ; if [[ ${#FILES[@]} > 0 ]] ; then cat ${FILES[@]} | bedtools sort -i - | bedtools merge -i - > ${SOURCE}/${CATEGORY}/${SAMPLE}/${SAMPLE}.${TOOL}.${SOURCE}.${CATEGORY}.BED.POS_jump.bed ; else touch ${SOURCE}/${CATEGORY}/${SAMPLE}/${SAMPLE}.${TOOL}.${SOURCE}.${CATEGORY}.BED.POS_jump.bed ; fi

FILES=( $( find ${SOURCE}/${CATEGORY}/${SAMPLE}/BED/ -iname "${SAMPLE}.${SOURCE}*${TOOL}_${SOURCE}_*.reject.extract.bed" ) ) ; if [[ ${#FILES[@]} > 0 ]] ; then cat ${FILES[@]} | bedtools sort -i - | bedtools merge -i - > ${SOURCE}/${CATEGORY}/${SAMPLE}/${SAMPLE}.${TOOL}.${SOURCE}.${CATEGORY}.BED.reject_1.bed ; else touch ${SOURCE}/${CATEGORY}/${SAMPLE}/${SAMPLE}.${TOOL}.${SOURCE}.${CATEGORY}.BED.reject_1.bed ; fi

FILES=( $( find ${SOURCE}/${CATEGORY}/${SAMPLE}/BED/ -iname "${SAMPLE}.${SOURCE}*${TOOL}_${TARGET}_*.reject.extract.bed" ) ) ; if [[ ${#FILES[@]} > 0 ]] ; then cat ${FILES[@]} | bedtools sort -i - | bedtools merge -i - > ${SOURCE}/${CATEGORY}/${SAMPLE}/${SAMPLE}.${TOOL}.${SOURCE}.${CATEGORY}.BED.reject_2.bed ; else touch ${SOURCE}/${CATEGORY}/${SAMPLE}/${SAMPLE}.${TOOL}.${SOURCE}.${CATEGORY}.BED.reject_2.bed ; fi

```


## 5.4&nbsp; Compare Data
The VCF-based data should be contained within the BED-based data, with differences arising only due to mis-matching reference alleles. 
To confirm this, we take the intersection of each of the categories from the BED-based data with the VCF-based data, and take the jaccard index of this intersection with the VCF-based data.
This proportion overlap should be equal to 1 for all files, indicating that the VCF-based data is a subset of the BED-based data. 
```
for BED in ${SOURCE}/${CATEGORY}/${SAMPLE}/${SAMPLE}.${TOOL}.${SOURCE}.${CATEGORY}.BED.*.bed ; do VCF=$( echo ${BED} | sed -e 's/BED/VCF/g' ) ; echo -e "${BED}\t$(bedtools intersect -a ${VCF} -b ${BED} | bedtools sort -i - | bedtools jaccard -a - -b ${VCF} | cut -f3 | tail -1 )" ; done | column -t

```

Similarly, we can confirm that the BED-based data is contained within the full-genome data. 
Note that for the `stable` data, the novel CUP files should be empty, so the jaccard index will be `-nan`. 
Otherwise, these should all be equal to 1. 
```
for BED in ${SOURCE}/${CATEGORY}/${SAMPLE}/${SAMPLE}.${TOOL}.${SOURCE}.${CATEGORY}.BED.*.bed ; do CAT=$( basename ${BED} | cut -f6 -d. ) ; echo -e "${BED}\t$(bedtools intersect -a ${BED} -b ${MAIN_DIR}/COMBINE/FASTA_BED.${TOOL}.ALL_${SOURCE}.${CAT}.bed | bedtools sort -i - | bedtools jaccard -a - -b ${BED} | cut -f3 | tail -1 )" ; done | column -t

```

We can confirm that pre-excluding variants at novel CUP positions (i.e. the `stable` data) prior to conversion results in the same list of variants as applying the algorithm to the original, unfiltered data, and removing variants at novel CUPs. 
In the output for the following, note that the Venn-Diagram numbers should indeicate all variants are shared by both files. 
Alternatively, in the Genotype Comparison Summary, the non-reference discordance rate should be zero. 
```
bgzip -c ${SOURCE}/original/${SAMPLE}/VCF/${SAMPLE}.${SOURCE}.annotate.bi_SNV.original_${TOOL}_${TARGET}_2.pass.vcf > ${SOURCE}/original/${SAMPLE}/VCF/${SAMPLE}.${SOURCE}.annotate.bi_SNV.original_${TOOL}_${TARGET}_2.pass.vcf.gz
tabix -f ${SOURCE}/original/${SAMPLE}/VCF/${SAMPLE}.${SOURCE}.annotate.bi_SNV.original_${TOOL}_${TARGET}_2.pass.vcf.gz

bgzip -c ${SOURCE}/stable/${SAMPLE}/VCF/${SAMPLE}.${SOURCE}.annotate.bi_SNV.stable_${TOOL}_${TARGET}_1.pass.vcf > ${SOURCE}/stable/${SAMPLE}/VCF/${SAMPLE}.${SOURCE}.annotate.bi_SNV.stable_${TOOL}_${TARGET}_1.pass.vcf.gz
tabix -f ${SOURCE}/stable/${SAMPLE}/VCF/${SAMPLE}.${SOURCE}.annotate.bi_SNV.stable_${TOOL}_${TARGET}_1.pass.vcf.gz

vcf-compare -g ${SOURCE}/original/${SAMPLE}/VCF/${SAMPLE}.${SOURCE}.annotate.bi_SNV.original_${TOOL}_${TARGET}_2.pass.vcf.gz ${SOURCE}/stable/${SAMPLE}/VCF/${SAMPLE}.${SOURCE}.annotate.bi_SNV.stable_${TOOL}_${TARGET}_1.pass.vcf.gz | grep ^SN | cut -f 2-

```

