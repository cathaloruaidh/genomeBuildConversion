#! /bin/bash

CHR=${1};
START=${2};  
END=${3}; 

while IFS="" read -r LINE || [[ -n "${LINE}" ]]
do
	mkdir FASTA_BED.${CHR} ; 

	time for (( i=${START}; i<${END}; i++ ))
	do 
		echo -e "${CHR}\t${i}\t$((i+1))\tFASTA_BED_${CHR}_${i}" ; 
	done > FASTA_BED.${CHR}/FASTA_BED.${CHR}.bed ; 

	echo ${CHR} ; 
done 