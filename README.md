# 1&nbsp; Genome Build Conversion
Code for identifying regions of the genome that are unstable when converting between hg19 and GRCh38, using liftOver or CrossMap. 



# 2&nbsp; Set Up
## 2.1&nbsp; Notes
- Prerequisites: `liftOver`, `CrossMap`, `bedtools` and `GNU parallel`. Binary files for `picard` are supplied. 
- Reference FASTA files are not included due to file size, but are required for the application to the real WGS data and should be stored in the REFERENCE directory. 
- This process assumes `chr1, chr2, ..., chrX, chrY, chrM` nomenclature. 
- The input BED files for the full-genome search for one build are ~150GB in size. Once the algorithm is applied, all files can take up to 1.5TB in size. 
- The code below needs to be run separately for both builds (hg19 and GRCh38) as well as using both tools (liftOver and CrossMap), so one of each should be selected. 



## 2.2&nbsp; Installation
Download the resource material and initialise:
```
git clone https://github.com/cathaloruaidh/genomeBuildConversion.git
cd genomeBuildConversion ;
REF_DIR=$( pwd )/REFERENCE
mkdir CHR COMBINE ;
chmod +x ${REF_DIR}/*sh
```

Run for hg19: 
```
SOURCE=hg19
TARGET=GRCh38
REGIONS=${REF_DIR}/ucsc.hg19.region.Standard.sort.bed
```

Run for GRCh38: 
```
SOURCE=GRCh38
TARGET=hg19
REGIONS=${REF_DIR}/GRCh38_full_analysis_set_plus_decoy_hla.regions.Standard.bed 
```


Run for liftOver
```
TOOL=liftOver
LOOP_BED=${REF_DIR}/loopLift_BED.sh
```

Run for CrossMap
```
TOOL=CrossMap
LOOP_BED=${REF_DIR}/loopCrossMap_BED.sh
```





# 3&nbsp; Full Genome Data
## 3.1&nbsp; Create Input BED

Generate the input BED files for the conversion process. 
Every individual base-pair position in the genome will have a BED entry, based on the lengths of the standard 23 pairs of chromosomes, including the mitochondrial chromosome.
This is parallelised for speed, as it can take up to 90 minutes per chromosome, depending on the size. 

```
cd CHR 
date ;
parallel --plus -j 12 --colsep '\t' "${REF_DIR}/createInputBed.sh {1} {2} {3}" :::: ${REGIONS}
```


## 3.2&nbsp; Apply Algorithm
Run the main script to identify unstable regions. 
The loop script takes as arguments the input filename, the start iteration, the end iteration and the source build. 
Both scripts will add the tool name to the file output, so there should be no over-writing of output files. 
This is parallelised using GNU parallel, with 12 CPUs available. 


```
parallel --plus -j12  "${LOOP_BED} {} 1 2 ${SOURCE}" ::: $( ls )
```

The script was set up so that iterations could be interrupted and restarted if neccessary. 
Two iterations were run abvove to determine if the algorithm was stable, or if new sites would be identified at each step (the former was expected and observed). 



## 3.3&nbsp; Sanity Check
Check if there are entries in the files of unstable regions for the second iteration by counting the lines (regardless of the source/target builds). 
The first two commands should return zeroes for all files, and the third command should return nothing. 

```
for FILE in $( find . -iname '*hg19_2.reject.extract.bed' ) ; do wc -l ${FILE} ; done
for FILE in $( find . -iname '*GRCh38_2.reject.extract.bed' ) ; do wc -l ${FILE} ; done
for FILE in $( find . -iname '*_2.jump*' ) ; do wc -l ${FILE} ; done
cd ../ ; 
```


## 3.4&nbsp; Combine Sites
For each of the five unstable regions, combine all the individual base-pair sites and collapse into multi-site regions. 
Additionally, combine the four novel unstable regions into one file. 

```
cd ../COMBINE ; 

cat $( find ../CHR/ -iname "*${TARGET}_*jump_CHR.bed" ) | bedtools sort -i - | bedtools merge -i - > FASTA_BED.${TOOL}.ALL_${SOURCE}.CHR_jump_1.bed
cat $( find ../CHR/ -iname "*${SOURCE}_*jump_CHR.bed" ) | bedtools sort -i - | bedtools merge -i - > FASTA_BED.${TOOL}.ALL_${SOURCE}.CHR_jump_2.bed
cat $( find ../CHR/ -iname "*${SOURCE}_*jump_POS.bed" ) | bedtools sort -i - | bedtools merge -i - > FASTA_BED.${TOOL}.ALL_${SOURCE}.POS_jump.bed
cat $( find ../CHR/ -iname "*${SOURCE}_*.reject.extract.bed" ) | bedtools sort -i - | bedtools merge -i - > FASTA_BED.${TOOL}.ALL_${SOURCE}.reject_1.bed
cat $( find ../CHR/ -iname "*${TARGET}_*.reject.extract.bed" ) | bedtools sort -i - | bedtools merge -i - > FASTA_BED.${TOOL}.ALL_${SOURCE}.reject_2.bed
cat FASTA_BED.${TOOL}.ALL_${SOURCE}*jump*.bed FASTA_BED.${TOOL}.ALL_${SOURCE}*reject_2.bed | bedtools sort -i - | bedtools merge -i - > FASTA_BED.${TOOL}.ALL_${SOURCE}.novel_exclude.bed

```




