#! /bin/bash


FILE=${1}
INPUT=${FILE}/${FILE}
PREFIX=${INPUT}.CrossMap

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
	TO=10
fi



if [[ "${4}" = "GRCh37" ]]
then
	SOURCE="GRCh37"
	TARGET="GRCh38"

	CHAIN_SOURCE_TO_TARGET=${REF_DIR}/hg19ToHg38.over.chain.gz
	CHAIN_TARGET_TO_SOURCE=${REF_DIR}/hg38ToHg19.over.chain.gz
else
	SOURCE="GRCh38"
	TARGET="GRCh37"

	CHAIN_SOURCE_TO_TARGET=${REF_DIR}/hg38ToHg19.over.chain.gz
	CHAIN_TARGET_TO_SOURCE=${REF_DIR}/hg19ToHg38.over.chain.gz
fi


if [[ ! -z "${5}" ]]
then
    OUT_DIR=${5}
else
    OUT_DIR=${PWD}
fi


cp ${INPUT}.bed ${OUT_DIR}/${PREFIX}_${SOURCE}_0.pass.bed

echo "Beginning for ${FILE}"

for i in $(seq ${FROM} ${TO} )
do
	echo -e "\n"
	date
	echo "Beginning iteration $i"
	echo -e "\n"


	p=$(( ${i} - 1 ))	

	SOURCE_PASS_PREV=${OUT_DIR}/${PREFIX}_${SOURCE}_${p}.pass.bed
	SOURCE_OUT=${OUT_DIR}/${PREFIX}_${SOURCE}_${i}.out.bed
	SOURCE_PASS=${OUT_DIR}/${PREFIX}_${SOURCE}_${i}.pass.bed
	SOURCE_REJECT=${OUT_DIR}/${PREFIX}_${SOURCE}_${i}.reject.bed
	SOURCE_EXTRACT=${OUT_DIR}/${PREFIX}_${SOURCE}_${i}.reject.extract.bed
	SOURCE_JUMP_CHR=${OUT_DIR}/${PREFIX}_${SOURCE}_${i}.jump_CHR.bed
	SOURCE_JUMP_POS=${OUT_DIR}/${PREFIX}_${SOURCE}_${i}.jump_POS.bed

	TARGET_OUT=${OUT_DIR}/${PREFIX}_${TARGET}_${i}.out.bed
	TARGET_PASS=${OUT_DIR}/${PREFIX}_${TARGET}_${i}.pass.bed
	TARGET_REJECT=${OUT_DIR}/${PREFIX}_${TARGET}_${i}.reject.bed
	TARGET_EXTRACT=${OUT_DIR}/${PREFIX}_${TARGET}_${i}.reject.extract.bed
	TARGET_JUMP_CHR=${OUT_DIR}/${PREFIX}_${TARGET}_${i}.jump_CHR.bed



	# lift the BED from SOURCE to TARGET
	CrossMap.py bed ${CHAIN_SOURCE_TO_TARGET} ${SOURCE_PASS_PREV} > ${TARGET_OUT}

	# Get the PASS file from the OUT file
	echo "Get PASS and JUMPS"
	grep -v 'Unmap' ${TARGET_OUT} | \
	awk -v OFS="\t" -v PASS=${TARGET_PASS} -v JUMPCHR=${TARGET_JUMP_CHR} '{ if($1 != $6) {print $1,$2,$3,$4 > JUMPCHR} else {print $6,$7,$8,$9 > PASS} }'

	# Get the EXTRACT file from the OUT file
	echo "Get REJECT"
	grep 'Unmap' ${TARGET_OUT} | \
	awk -v OFS="\t" '{tmp=$4 ; gsub(/_/, "\t", $4) ; print $1,$2,$3,tmp,$4}' | \
	awk -v OFS="\t" '{print $7,$8,$8+1,$4}' > ${SOURCE_EXTRACT}


	echo -e "\nFinished with ${SOURCE} to ${TARGET}. \n"



	# lift the BED back from TARGET to SOURCE
	CrossMap.py bed ${CHAIN_TARGET_TO_SOURCE} ${TARGET_PASS} > ${SOURCE_OUT}

	# Get the PASS and JUMP files from the OUT file
	grep -v 'Unmap' ${SOURCE_OUT} | \
	cut -f6-9 | \
	awk -v OFS="\t" '{tmp=$4 ; gsub(/_/, "\t", $4) ; print $1,$2,$3,tmp,$4}' | \
	awk -v OFS="\t" -v PASS=${SOURCE_PASS} -v JUMPCHR=${SOURCE_JUMP_CHR} -v JUMPPOS=${SOURCE_JUMP_POS}  '{ if($1 != $7) {print $7,$8,$8+1,$4 > JUMPCHR}  else if($2 != $8) {print $7,$8,$8+1,$4 > JUMPPOS} else {print $7,$8,$8+1,$4 > PASS} }'

	# Get the EXTRACT file from the OUT file
	grep 'Unmap' ${SOURCE_OUT} | \
	awk -v OFS="\t" '{tmp=$4 ; gsub(/_/, "\t", $4) ; print $1,$2,$3,tmp,$4}' | \
	awk -v OFS="\t" -v N=${i} -v REF=${TARGET} '{print $7,$8,$8+1,$4 "_" REF "_-" N}' > ${TARGET_EXTRACT}


	echo -e "\nFinished with ${TARGET} to ${SOURCE}. \n"
	

done


echo -e "\n\n\n"