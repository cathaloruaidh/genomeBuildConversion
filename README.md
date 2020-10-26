# 1&nbsp; Genome Build Conversion
Code for identifying positions in the human genome that are unstable when converting between GRCh37 and GRCh38, using either liftOver or CrossMap. 
Previous work has highlighted unusual behaviour in build conversion, such as SNVs mapping to a different chromosome. 
BED files describing these conversion-unstable positions (CUPs) are provided, as well as a combined file of all regions ("novel_CUPs"). 
Pre-excluding SNVs in these regions before converting between builds removes all unstable behaviour. 

If you have any queries or feedback, please contact the [author](mailto:cathalormond@gmail.com). If you use the exclude files or this algorithm in your publication, please cite the following paper:

Navigation: 
- [Set Up](#2-set-up) 
- [Full Genome Data](#3-full-genome-data)
- [Real WGS Example](#4-real-wgs-example)


# 2&nbsp; Set Up
## 2.1&nbsp; Notes
- Prerequisites: `liftOver`, `CrossMap`, `bedtools` and `GNU parallel`. 
- Reference FASTA files are not included due to file size, but are required for the application to the real WGS data and should be stored in the REFERENCE directory. 
- This process assumes `chr1, chr2, ..., chrX, chrY` nomenclature. 
- The input BED files for the full-genome search for one build are ~150GB in size. Once the algorithm is applied, all files can take up to 1.5TB in size. 
- For a single run of the algorithm, a source and target such as GRCh37 and GRCh38 must be chosen as well as a tool (liftOver or CrossMap). 



## 2.2&nbsp; Installation
Download the resource material and initialise:
```
git clone https://github.com/cathaloruaidh/genomeBuildConversion.git
cd genomeBuildConversion ;

REF_DIR=$( pwd )/REFERENCE
SCRIPT_DIR=$( pwd )/SCRIPTS
export REF_DIR
export SCRIPT_DIR
mkdir CHR COMBINE ;
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


Run for liftOver
```
TOOL=liftOver
export TOOL
LOOP_BED=${SCRIPT_DIR}/loopLift_BED.sh
```

Run for CrossMap
```
TOOL=CrossMap
export TOOL
LOOP_BED=${SCRIPT_DIR}/loopCrossMap_BED.sh
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
Run the main script to identify unstable regions. 
The loop script takes as arguments the input filename, the start iteration, the end iteration and the source build. 
Both scripts will add the tool name to the file output, so there should be no over-writing of output files. 
This is parallelised using GNU parallel, with 12 CPUs available. 


```
parallel --plus -j12  ". ${LOOP_BED} {} 1 2 ${SOURCE}" ::: $( ls | sort -V)
```

The script was set up so that iterations could be interrupted and restarted if neccessary. 
Two iterations were run above to determine if the algorithm was stable, or if new sites would be identified at each step (the former was expected and observed). 



## 3.3&nbsp; Sanity Check
Check if there are entries in the files of unstable regions for the second iteration by counting the lines (regardless of the source/target builds). 
The first two commands should return zeroes for all files, and the third command should return nothing. 

```
for FILE in $( find . -iname '*GRCh37_2.reject.extract.bed' ) ; do wc -l ${FILE} ; done
for FILE in $( find . -iname '*GRCh38_2.reject.extract.bed' ) ; do wc -l ${FILE} ; done
for FILE in $( find . -iname '*_2.jump*' ) ; do wc -l ${FILE} ; done
cd ../ ; 
```


## 3.4&nbsp; Combine Sites
For each of the five unstable regions, combine all the individual base-pair positions and collapse into multi-site regions. 
Additionally, combine the four novel unstable regions into one file. 

```
cd ../COMBINE ; 

cat $( find ../CHR/ -iname "*${TARGET}_*jump_CHR.bed" ) | bedtools sort -i - | bedtools merge -i - > FASTA_BED.${TOOL}.ALL_${SOURCE}.CHR_jump_1.bed
cat $( find ../CHR/ -iname "*${SOURCE}_*jump_CHR.bed" ) | bedtools sort -i - | bedtools merge -i - > FASTA_BED.${TOOL}.ALL_${SOURCE}.CHR_jump_2.bed
cat $( find ../CHR/ -iname "*${SOURCE}_*jump_POS.bed" ) | bedtools sort -i - | bedtools merge -i - > FASTA_BED.${TOOL}.ALL_${SOURCE}.POS_jump.bed
cat $( find ../CHR/ -iname "*${SOURCE}_*.reject.extract.bed" ) | bedtools sort -i - | bedtools merge -i - > FASTA_BED.${TOOL}.ALL_${SOURCE}.reject_1.bed
cat $( find ../CHR/ -iname "*${TARGET}_*.reject.extract.bed" ) | bedtools sort -i - | bedtools merge -i - > FASTA_BED.${TOOL}.ALL_${SOURCE}.reject_2.bed
cat FASTA_BED.${TOOL}.ALL_${SOURCE}*jump*.bed FASTA_BED.${TOOL}.ALL_${SOURCE}*reject_2.bed | bedtools sort -i - | bedtools merge -i - > FASTA_BED.${TOOL}.ALL_${SOURCE}.novel_CUPs.bed

```




# 4 &nbsp;Real WGS Example
Get the VCF files for NA12877 and NA12878 from the Illumina Platinum Genomes project. 

```
mkdir WGS_Data
cd WGS_Data
mkdir GRCh37 GRCh38

wget https://s3.eu-central-1.amazonaws.com/platinum-genomes/2017-1.0/hg19/small_variants/ConfidentRegions.bed.gz -O GRCh37/ConfidentRegions.bed.gz
wget https://s3.eu-central-1.amazonaws.com/platinum-genomes/2017-1.0/hg19/small_variants/ConfidentRegions.bed.gz.tbi -O GRCh37/ConfidentRegions.bed.gz.tbi
wget https://s3.eu-central-1.amazonaws.com/platinum-genomes/2017-1.0/hg19/small_variants/NA12877/NA12877.vcf.gz -O GRCh37/NA12877.GRCh37.vcf.gz
wget https://s3.eu-central-1.amazonaws.com/platinum-genomes/2017-1.0/hg19/small_variants/NA12877/NA12877.vcf.gz.tbi -O GRCh37/NA12877.GRCh37.vcf.gz.tbi
wget https://s3.eu-central-1.amazonaws.com/platinum-genomes/2017-1.0/hg19/small_variants/NA12878/NA12878.vcf.gz -O GRCh37/NA12878.GRCh37.vcf.gz
wget https://s3.eu-central-1.amazonaws.com/platinum-genomes/2017-1.0/hg19/small_variants/NA12878/NA12878.vcf.gz.tbi -O GRCh37/NA12878.GRCh37.vcf.gz.tbi

wget https://s3.eu-central-1.amazonaws.com/platinum-genomes/2017-1.0/hg38/small_variants/ConfidentRegions.bed.gz -O GRCh38/ConfidentRegions.bed.gz
wget https://s3.eu-central-1.amazonaws.com/platinum-genomes/2017-1.0/hg38/small_variants/ConfidentRegions.bed.gz.tbi -O GRCh38/ConfidentRegions.bed.gz.tbi
wget https://s3.eu-central-1.amazonaws.com/platinum-genomes/2017-1.0/hg38/small_variants/NA12877/NA12877.vcf.gz -O GRCh38/NA12877.GRCh37.vcf.gz
wget https://s3.eu-central-1.amazonaws.com/platinum-genomes/2017-1.0/hg38/small_variants/NA12877/NA12877.vcf.gz.tbi  -O GRCh38/NA12877.GRCh37.vcf.gz.tbi
wget https://s3.eu-central-1.amazonaws.com/platinum-genomes/2017-1.0/hg38/small_variants/NA12878/NA12878.vcf.gz -O GRCh38/NA12878.GRCh38.vcf.gz
wget https://s3.eu-central-1.amazonaws.com/platinum-genomes/2017-1.0/hg38/small_variants/NA12878/NA12878.vcf.gz.tbi -O GRCh38/NA12878.GRCh38.vcf.gz.tbi
```






