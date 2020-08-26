#! /bin/bash

#FILE=NA12878.vcf
#PREFIX=NA12878




FILE=${1}
PREFIX=${FILE}/${FILE}.liftOver


if [[ ! -z "${2}" ]]
then
	FROM=${2}
else
	FROM=1
fi


if [[ ! -z "${3}" ]]
then
	TO=${3}
else
	TO=100
fi



if [[ "${4}" = "hg19" ]]
then
	SOURCE="hg19"
	TARGET="GRCh38"

	CHAIN_SOURCE_TO_TARGET=/home/shared/cathal/reference/chainFiles/GRCh37_to_GRCH38.chain.gz
	CHAIN_TARGET_TO_SOURCE=/home/shared/cathal/reference/chainFiles/hg38ToHg19.over.chain.gz

	REF_SOURCE=/home/shared/reference/ReferenceGenome/hg19/ucsc.hg19.fasta
	REF_TARGET=/home/shared/reference/ReferenceGenome/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.fa
else
	SOURCE="GRCh38"
	TARGET="hg19"

	CHAIN_SOURCE_TO_TARGET=/home/shared/cathal/reference/chainFiles/hg38ToHg19.over.chain.gz
	CHAIN_TARGET_TO_SOURCE=/home/shared/cathal/reference/chainFiles/GRCh37_to_GRCH38.chain.gz

	REF_SOURCE=/home/shared/reference/ReferenceGenome/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.fa
	REF_TARGET=/home/shared/reference/ReferenceGenome/hg19/ucsc.hg19.fasta
fi



cp ${PREFIX}.bed ${PREFIX}_${SOURCE}_0.pass.bed


for i in $(seq ${FROM} ${TO} )
do
	echo -e "\n\n"
	date
	echo "Beginning iteration $i"
	echo -e "\n\n"


	p=$(( ${i} - 1 ))	

	SOURCE_PASS_PREV=${PREFIX}_${SOURCE}_${p}.pass.bed
	SOURCE_OUT=${PREFIX}_${SOURCE}_${i}.out.bed
	SOURCE_PASS=${PREFIX}_${SOURCE}_${i}.pass.bed
	SOURCE_REJECT=${PREFIX}_${SOURCE}_${i}.reject.bed
	SOURCE_EXTRACT=${PREFIX}_${SOURCE}_${i}.reject.extract.bed
	SOURCE_JUMP_CHR=${PREFIX}_${SOURCE}_${i}.jump_CHR.bed
	SOURCE_JUMP_POS=${PREFIX}_${SOURCE}_${i}.jump_POS.bed

	TARGET_OUT=${PREFIX}_${TARGET}_${i}.out.bed
	TARGET_PASS=${PREFIX}_${TARGET}_${i}.pass.bed
	TARGET_REJECT=${PREFIX}_${TARGET}_${i}.reject.bed
	TARGET_EXTRACT=${PREFIX}_${TARGET}_${i}.reject.extract.bed
	TARGET_JUMP_CHR=${PREFIX}_${TARGET}_${i}.jump_CHR.bed



	# lift the VCF from SOURCE to TARGET
	liftOver ${SOURCE_PASS_PREV} ${CHAIN_SOURCE_TO_TARGET} ${TARGET_OUT} ${SOURCE_REJECT} 

	# Convert the REJECT to an EXTRACT 
	echo "Get Reject"
	grep -v '^#' ${SOURCE_REJECT} | \
	awk -v OFS="\t" '{tmp=$4 ; gsub(/_/, "\t", $4) ; print $1,$2,$3,tmp,$4}' | \
	awk -v OFS="\t" '{print $7,$8,$8+1,$4}' > ${SOURCE_EXTRACT}

	echo "Get Jump"
	awk -v OFS="\t" '{tmp=$4 ; gsub(/_/, "\t", $4) ; print $1,$2,$3,tmp,$4}' ${TARGET_OUT} | \
	awk -v OFS="\t" -v PASS=${TARGET_PASS} -v JUMP=${TARGET_JUMP_CHR} '{if($1 != $7) { print $7,$8,$8+1,$4  > JUMP} else { print $1,$2,$3,$4 > PASS }}'
	

	echo -e "\nFinished with ${SOURCE} to ${TARGET}. \n"



	# list the VCF back from TARGET to SOURCE
	liftOver ${TARGET_PASS} ${CHAIN_TARGET_TO_SOURCE} ${SOURCE_OUT} ${TARGET_REJECT}

	# Get the PASS and JUMP file
	echo "Get Jump"
	awk -v OFS="\t" '{tmp=$4 ; gsub(/_/, "\t", $4) ; print $1,$2,$3,tmp,$4}' ${SOURCE_OUT} | \
	awk -v OFS="\t" -v PASS=${SOURCE_PASS} -v JUMPCHR=${SOURCE_JUMP_CHR} -v JUMPPOS=${SOURCE_JUMP_POS} '{if($1 != $7) {print $7,$8,$8+1,$4 > JUMPCHR} else if($2 != $8) {print $7,$8,$8+1,$4 > JUMPPOS} else {print $1,$2,$3,$4 > PASS} }'

	# Convert the REJECT to an EXTRACT 
	echo "Get Reject"
	grep -v '^#' ${TARGET_REJECT} | \
	awk -v OFS="\t" '{tmp=$4 ; gsub(/_/, "\t", $4) ; print $1,$2,$3,tmp,$4}' | \
	awk -v OFS="\t" '{print $7,$8,$8+1,$4}' > ${TARGET_EXTRACT}


	echo -e "\nFinished with ${TARGET} to ${SOURCE}. \n"
	

done
