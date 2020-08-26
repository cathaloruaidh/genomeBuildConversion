# Genome Build Conversion
Code for identifying regions of the genome that are unstable when converting between hg19 and GRCh38, using liftOver or CrossMap. 

## Installation
Download the resource material 
`git clone https://github.com/cathaloruaidh/genomeBuildConversion.git`


## Set Up
Directory and script variables to be set prior to running the code. 

```
MAIN_DIR=$( pwd )
REF_DIR=${MAIN_DIR}/REFERENCE

REGIONS_38=${REF_DIR}/GRCh38_full_analysis_set_plus_decoy_hla.regions.Standard.bed 
REGIONS_19=${REF_DIR}/ucsc.hg19.region.Standard.sort.bed
```

Note: Need to select this depending on whether liftOver or CrossMap is used. 
```
TOOL=liftOver
LOOP_BED=${REF_DIR}/loopLift_BED.sh

TOOL=CrossMap
LOOP_BED=${REF_DIR}/loopCrossMap_BED.sh
```

Create the directories
```
mkdir CHR COMBINE ;
mkdir CHR/hg19 CHR/GRCh38; 
```


## Resrequisites/Notes
- To be installed: `liftOver`, `CrossMap` and `bedtools`. 
- Reference FASTA files are not included due to file size, but are required for the 
- This process assumes `chr1, chr2, ..., chrX, chrY, chrM` nomenclature. 
- The input BED files for the full-genome search are ~150GB in size. 


# Code Main
## Create Input BED

Generate the input BED files for the conversion process. 
Every individual base-pair position in the genome will have a BED entry, based on the lengths of the standard 23 pairs of chromosomes, including the mitochondrial chromosome. 

Run for hg19: 
```
cd hg19
date 
while IFS="" read -r LINE || [[ -n "${LINE}" ]]
do
	CHR=$( echo ${LINE} | cut -f1 -d ' '  ) ;
	START=$( echo ${LINE} | cut -f2 -d ' '  ) ;  
	END=$( echo ${LINE} | cut -f3 -d ' '  ) ; 

	mkdir FASTA_BED.${CHR} ; 

	time for (( i=${START}; i<${END}; i++ ))
	do 
		echo -e "${CHR}\t${i}\t$((i+1))\tFASTA_BED_${CHR}_${i}" ; 
	done > FASTA_BED.${CHR}/FASTA_BED.${CHR}.bed ; 

	echo ${CHR} ; 
done < ${REGIONS_19}
cd ../
```

Run for GRCh38: 
```
cd GRCh38 
date 
while IFS="" read -r LINE || [[ -n "${LINE}" ]]
do
	CHR=$( echo ${LINE} | cut -f1 -d ' '  ) ;  
	START=$( echo ${LINE} | cut -f2 -d ' '  ) ;  
	END=$( echo ${LINE} | cut -f3 -d ' '  ) ; 

	mkdir FASTA_BED.${CHR} ; 

	time for (( i=${START}; i<${END}; i++ )) ; do 
		echo -e "${CHR}\t${i}\t$((i+1))\tFASTA_BED_${CHR}_${i}" ; 
	done > FASTA_BED.${CHR}/FASTA_BED.${CHR}.bed ; 

	echo ${CHR} ; 
done < ${REGIONS_38}
cd ../ ; 
```


## Apply Algorithm
Run the main script to identify unstable regions. 
The loop script takes as arguments the input filename, the start iteration, the end iteration and the source build. 
Both scripts will add the tool name to the file output, so there should be no over-writing of output files. 
This is parallelised using GNU parallel, with 12 CPUs available. 
Two iterations were run to determine if the algorithm was stable, or if new sites would be identified at each step (the former was expected). 


Run for hg19: 
```
cd hg19 
parallel --plus -j12  "${LOOP_BED} {} 1 2 hg19" ::: $( ls )
cd ../ 
```

Run for GRCh38: 
```
cd GRCh38 
parallel --plus -j12  "${LOOP_BED} {} 1 2 GRCh38" ::: $( ls )
cd ../ 
```


## Sanity Check
Check if there are entries in the files of unstable regions for the second iteration by counting the lines (regardless of the source/target builds). 
The first two commands should return zeroes for all files, and the final line should return nothing. 

```
for FILE in $( find . -iname '*hg19_2.reject.extract.bed' ) ; do wc -l ${FILE} ; done
for FILE in $( find . -iname '*GRCh38_2.reject.extract.bed' ) ; do wc -l ${FILE} ; done
for FILE in $( find . -iname '*_2.jump*' ) ; do wc -l ${FILE} ; done
cd ../ ; 
```


## Combine
Run for hg19: 
```
SOURCE=GRCh38
TARGET=hg19
```

Run for GRCh38:
```
SOURCE=hg19
TARGET=GRCh38
```


Combine all the individual base-pair sites and collapse into multi-site regions. 
```
cd ../COMBINE ; 

cat $( find ../CHR/ -iname "*${TARGET}_*jump_CHR.bed" ) | bedtools sort -i - | bedtools merge -i - > FASTA_BED.${TOOL}.ALL_${SOURCE}.CHR_jump_1.bed
cat $( find ../CHR/ -iname "*${SOURCE}_*jump_CHR.bed" ) | bedtools sort -i - | bedtools merge -i - > FASTA_BED.${TOOL}.ALL_${SOURCE}.CHR_jump_2.bed
cat $( find ../CHR/ -iname "*${SOURCE}_*jump_POS.bed" ) | bedtools sort -i - | bedtools merge -i - > FASTA_BED.${TOOL}.ALL_${SOURCE}.POS_jump.bed
cat $( find ../CHR/ -iname "*${SOURCE}_*.reject.extract.bed" ) | bedtools sort -i - | bedtools merge -i - > FASTA_BED.${TOOL}.ALL_${SOURCE}.reject_1.bed
cat $( find ../CHR/ -iname "*${TARGET}_*.reject.extract.bed" ) | bedtools sort -i - | bedtools merge -i - > FASTA_BED.${TOOL}.ALL_${SOURCE}.reject_2.bed
cat FASTA_BED.${TOOL}.ALL_${SOURCE}*jump*.bed FASTA_BED.${TOOL}.ALL_${SOURCE}*reject_2.bed | bedtools sort -i - | bedtools merge -i - > FASTA_BED.${TOOL}.ALL_${SOURCE}.novel_exclude.bed

```




