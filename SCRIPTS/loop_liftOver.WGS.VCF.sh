#! /bin/bash



FILE=${1}
PREFIX=$( basename ${FILE%.vcf} )


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



if [[ "${4}" = "GRCh37" ]]
then
	SOURCE="GRCh37"
	TARGET="GRCh38"

	CHAIN_SOURCE_TO_TARGET=${REF_DIR}/hg19ToHg38.over.chain.gz
	CHAIN_TARGET_TO_SOURCE=${REF_DIR}/hg38ToHg19.over.chain.gz

	REF_SOURCE=${REF_DIR}/GRCh37.fa
	REF_TARGET=${REF_DIR}/GRCh38.fa
else
	SOURCE="GRCh38"
	TARGET="GRCh37"

	CHAIN_SOURCE_TO_TARGET=${REF_DIR}/hg38ToHg19.over.chain.gz
	CHAIN_TARGET_TO_SOURCE=${REF_DIR}/hg19ToHg38.over.chain.gz

	REF_SOURCE=${REF_DIR}/GRCh38.fa
	REF_TARGET=${REF_DIR}/GRCh37.fa
fi


if [[ ! -z "${5}" ]]
then
    OUT_DIR=${5}
else
    OUT_DIR=${PWD}
fi


cp ${FILE} ${OUT_DIR}/${PREFIX}_liftOver_${SOURCE}_0.pass.vcf


for i in $(seq ${FROM} ${TO} )
do
	echo -e "\n\n"
	date
	echo "Beginning iteration $i"
	echo -e "\n\n"


	p=$(( ${i} - 1 ))	

	SOURCE_PASS_PREV=${OUT_DIR}/${PREFIX}_liftOver_${SOURCE}_${p}.pass.vcf
	SOURCE_OUT=${OUT_DIR}/${PREFIX}_liftOver_${SOURCE}_${i}.out.vcf
	SOURCE_PASS=${OUT_DIR}/${PREFIX}_liftOver_${SOURCE}_${i}.pass.vcf
	SOURCE_PASS_BED=${OUT_DIR}/${PREFIX}_liftOver_${SOURCE}_${i}.pass.bed
	SOURCE_REJECT=${OUT_DIR}/${PREFIX}_liftOver_${SOURCE}_${i}.reject.vcf
	SOURCE_EXTRACT=${OUT_DIR}/${PREFIX}_liftOver_${SOURCE}_${i}.reject.extract.bed
	SOURCE_JUMP_CHR=${OUT_DIR}/${PREFIX}_liftOver_${SOURCE}_${i}.jump_CHR.bed
	SOURCE_JUMP_POS=${OUT_DIR}/${PREFIX}_liftOver_${SOURCE}_${i}.jump_POS.bed
	SOURCE_MISMATCH=${OUT_DIR}/${PREFIX}_liftOver_${SOURCE}_${i}.mismatch.bed
	SOURCE_REMOVE=${OUT_DIR}/${PREFIX}_liftOver_${SOURCE}_${i}.remove.txt

	TARGET_OUT=${OUT_DIR}/${PREFIX}_liftOver_${TARGET}_${i}.out.vcf
	TARGET_PASS=${OUT_DIR}/${PREFIX}_liftOver_${TARGET}_${i}.pass.vcf
	TARGET_PASS_BED=${OUT_DIR}/${PREFIX}_liftOver_${TARGET}_${i}.pass.bed
	TARGET_REJECT=${OUT_DIR}/${PREFIX}_liftOver_${TARGET}_${i}.reject.vcf
	TARGET_EXTRACT=${OUT_DIR}/${PREFIX}_liftOver_${TARGET}_${i}.reject.extract.bed
	TARGET_JUMP_CHR=${OUT_DIR}/${PREFIX}_liftOver_${TARGET}_${i}.jump_CHR.bed
	TARGET_MISMATCH=${OUT_DIR}/${PREFIX}_liftOver_${TARGET}_${i}.mismatch.bed
	TARGET_REMOVE=${OUT_DIR}/${PREFIX}_liftOver_${TARGET}_${i}.remove.txt



	# lift the VCF from SOURCE to TARGET

	java -Djava.io.tmpdir=${OUT_DIR} -jar ${SCRIPT_DIR}/picard.jar LiftoverVcf \
		I=${SOURCE_PASS_PREV} \
		O=${TARGET_OUT} \
		C=${CHAIN_SOURCE_TO_TARGET} \
		REJECT=${SOURCE_REJECT} \
		R=${REF_TARGET}
	
	echo -e "\n"

	# Convert the REJECT file to EXTRACT and MISMATCH files 
	echo -e "Get Reject\n\n"
    
    if [[ -f ${SOURCE_REJECT} ]]
    then
        grep -v "^#" ${SOURCE_REJECT} | \
        grep 'NoTarget' | \
        awk -v OFS="\t" '{print $1,$2-1,$2,$3}' | \
        awk -v OFS="\t" '{tmp=$4 ; gsub(/_/, "\t", $4) ; print $1,$2,$3,tmp,$4}' | \
        awk -v OFS="\t" '{print $5,$6-1,$6,$4}' > ${SOURCE_EXTRACT}
    else
        touch ${SOURCE_EXTRACT}
    fi


	echo -e "Get MisMatch\n\n"

    if [[ -f ${SOURCE_REJECT} ]]
    then
        grep -v "^#" ${SOURCE_REJECT} | \
        grep 'MismatchedRefAllele' | \
        awk -v OFS="\t" '{print $1,$2-1,$2,$3}' | \
        awk -v OFS="\t" '{tmp=$4 ; gsub(/_/, "\t", $4) ; print $1,$2,$3,tmp,$4}' | \
        awk -v OFS="\t" '{print $5,$6-1,$6,$4}' > ${SOURCE_MISMATCH}
    else
        ${SOURCE_MISMATCH}
    fi

	echo -e "Get ChrJump1\n\n"
	bcftools query -f "%CHROM\t%POS\t%ID\n" ${TARGET_OUT} | \
	awk -v OFS="\t" '{tmp=$3 ; gsub(/_/, "\t", $3) ; print $1,$2,tmp,$3}' | \
	awk -v OFS="\t" '{if($1 != $4) print $4,$5-1,$5,$3}' > ${TARGET_JUMP_CHR}

	cut -f4 ${TARGET_JUMP_CHR}  > ${TARGET_REMOVE}
	
	if [[ -f ${TARGET_REMOVE} && $( wc -l ${TARGET_REMOVE} | cut -f1 ) > 0 ]]
	then
		java -Djava.io.tmpdir=${OUT_DIR} -jar ${SCRIPT_DIR}/GenomeAnalysisTK_3.8.jar \
			-T SelectVariants \
			-V ${TARGET_OUT} \
			-R ${REF_TARGET} \
			--excludeIDs ${TARGET_REMOVE} \
			-o ${TARGET_PASS}
	else
		cp ${TARGET_OUT} ${TARGET_PASS}
	fi

	echo -e "\nFinished with ${SOURCE} to ${TARGET}. \n"



	# list the VCF back from TARGET to SOURCE
	java -Djava.io.tmpdir=${OUT_DIR} -jar ${SCRIPT_DIR}/picard.jar LiftoverVcf \
		I=${TARGET_PASS} \
		O=${SOURCE_OUT} \
		C=${CHAIN_TARGET_TO_SOURCE} \
		REJECT=${TARGET_REJECT} \
		R=${REF_SOURCE}
		   

	# Get the PASS and JUMP file
	echo -e "\n\nGet ChrJump2 and PosJump\n\n"
	grep -v "^#" ${SOURCE_OUT} | \
	awk -v OFS="\t" '{tmp=$3 ; gsub(/_/, "\t", $3) ; print $1,$2-1,$2,tmp,$3}' | \
	awk -v OFS="\t" -v PASS=${SOURCE_PASS_BED} -v JUMPCHR=${SOURCE_JUMP_CHR} -v JUMPPOS=${SOURCE_JUMP_POS} '{ if($1 != $5) {print $5,$6-1,$6,$4 > JUMPCHR}  else if($3 != $6) {print $5,$6-1,$6,$4 > JUMPPOS} else {print $1,$2,$3,$4 > PASS} }'

	cat ${SOURCE_JUMP_CHR} ${SOURCE_JUMP_POS} | cut -f4 > ${SOURCE_REMOVE}


	echo -e "Get Pass \n\n"
    
	if [[ -f ${SOURCE_REMOVE} && $( wc -l ${SOURCE_REMOVE} | cut -f1 ) > 0 ]]
	then
		java -Djava.io.tmpdir=${OUT_DIR} -jar ${SCRIPT_DIR}/GenomeAnalysisTK_3.8.jar \
			-T SelectVariants \
			-V ${SOURCE_OUT} \
			-R ${REF_SOURCE} \
			--excludeIDs ${SOURCE_REMOVE} \
			-o ${SOURCE_PASS}
	else
		cp ${SOURCE_OUT} ${SOURCE_PASS}
	fi


	# Convert the REJECT to an EXTRACT 
	echo -e "Get Reject\n\n"
    
    if [[ -f ${TARGET_REJECT} ]]
    then
        grep -v "^#" ${TARGET_REJECT} | grep 'NoTarget' | \
        awk -v OFS="\t" '{print $1,$2-1,$2,$3}' | \
        awk -v OFS="\t" '{tmp=$4 ; gsub(/_/, "\t", $4) ; print $1,$2,$3,tmp,$4}' | \
        awk -v OFS="\t" '{print $5,$6-1,$6,$4}' > ${TARGET_EXTRACT}
    else
        touch ${TARGET_EXTRACT}
    fi
    

	echo -e "Get Mismatch\n\n"
    
    if [[ -f ${TARGET_REJECT} ]]
    then
        grep -v "^#" ${TARGET_REJECT} | grep 'MismatchedRefAllele' | \
        awk -v OFS="\t" '{print $1,$2-1,$2,$3}' | \
        awk -v OFS="\t" '{tmp=$4 ; gsub(/_/, "\t", $4) ; print $1,$2,$3,tmp,$4}' | \
        awk -v OFS="\t" '{print $5,$6-1,$6,$4}' > ${TARGET_MISMATCH}
    else
        touch ${TARGET_MISMATCH}
    fi


	echo -e "\nFinished with ${TARGET} to ${SOURCE}. \n"
	

done
