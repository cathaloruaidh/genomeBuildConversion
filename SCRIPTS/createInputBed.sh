#! /bin/bash

CHR=${1};
START=${2};  
END=${3}; 
SOURCE=${4}

mkdir FASTA_BED.${CHR}_${SOURCE} ; 

time for (( i=${START}; i<${END}; i++ ))
do 
	echo -e "${CHR}\t${i}\t$((i+1))\tFASTA_BED_${CHR}_${i}" ; 
done > FASTA_BED.${CHR}_${SOURCE}/FASTA_BED.${CHR}_${SOURCE}.bed ; 

echo ${CHR} ; 