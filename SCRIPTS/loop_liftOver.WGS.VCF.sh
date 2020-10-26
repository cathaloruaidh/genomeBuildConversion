#! /bin/bash

#FILE=NA12878.vcf
#PREFIX=NA12878




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



if [[ "${4}" = "hg19" ]]
then
	SOURCE="hg19"
	TARGET="GRCh38"

	CHAIN_SOURCE_TO_TARGET=/home/shared/cathal/reference/chainFiles/GRCh37_to_GRCH38.chain.gz
	CHAIN_TARGET_TO_SOURCE=/home/shared/cathal/reference/chainFiles/hg38ToHg19.over.chain.gz

	REF_SOURCE=/home/shared/cathal/reference/ReferenceGenome/hg19/ucsc.hg19.fasta
	REF_TARGET=/home/shared/cathal/reference/ReferenceGenome/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.fa
else
	SOURCE="GRCh38"
	TARGET="hg19"

	CHAIN_SOURCE_TO_TARGET=/home/shared/cathal/reference/chainFiles/hg38ToHg19.over.chain.gz
	CHAIN_TARGET_TO_SOURCE=/home/shared/cathal/reference/chainFiles/GRCh37_to_GRCH38.chain.gz

	REF_SOURCE=/home/shared/cathal/reference/ReferenceGenome/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.fa
	REF_TARGET=/home/shared/cathal/reference/ReferenceGenome/hg19/ucsc.hg19.fasta
fi



cp ${PREFIX}.vcf ${PREFIX}_${SOURCE}_0.pass.vcf


for i in $(seq ${FROM} ${TO} )
do
	echo -e "\n\n"
	date
	echo "Beginning iteration $i"
	echo -e "\n\n"


	p=$(( ${i} - 1 ))	

	SOURCE_PASS_PREV=${PREFIX}_${SOURCE}_${p}.pass.vcf
	SOURCE_OUT=${PREFIX}_${SOURCE}_${i}.out.vcf
	SOURCE_PASS=${PREFIX}_${SOURCE}_${i}.pass.vcf
	SOURCE_PASS_BED=${PREFIX}_${SOURCE}_${i}.pass.bed
	SOURCE_REJECT=${PREFIX}_${SOURCE}_${i}.reject.vcf
	SOURCE_EXTRACT=${PREFIX}_${SOURCE}_${i}.reject.extract.bed
	SOURCE_JUMP_CHR=${PREFIX}_${SOURCE}_${i}.CHR_jump.bed
	SOURCE_JUMP_POS=${PREFIX}_${SOURCE}_${i}.POS_jump.bed
	SOURCE_MISMATCH=${PREFIX}_${SOURCE}_${i}.mismatch.bed
	SOURCE_REMOVE=${PREFIX}_${SOURCE}_${i}.remove.txt

	TARGET_OUT=${PREFIX}_${TARGET}_${i}.out.vcf
	TARGET_PASS=${PREFIX}_${TARGET}_${i}.pass.vcf
	TARGET_PASS_BED=${PREFIX}_${TARGET}_${i}.pass.bed
	TARGET_REJECT=${PREFIX}_${TARGET}_${i}.reject.vcf
	TARGET_EXTRACT=${PREFIX}_${TARGET}_${i}.reject.extract.bed
	TARGET_JUMP_CHR=${PREFIX}_${TARGET}_${i}.CHR_jump.bed
	TARGET_MISMATCH=${PREFIX}_${TARGET}_${i}.mismatch.bed
	TARGET_REMOVE=${PREFIX}_${TARGET}_${i}.remove.txt



	# lift the VCF from SOURCE to TARGET

	java -Djava.io.tmpdir=/home/shared/cathal/tmp -jar /home/shared/cathal/tools/picard.jar LiftoverVcf \
		I=${SOURCE_PASS_PREV} \
		O=${TARGET_OUT} \
		C=${CHAIN_SOURCE_TO_TARGET} \
		REJECT=${SOURCE_REJECT} \
		R=${REF_TARGET}
	
	echo -e "\n"

	# Convert the REJECT file to EXTRACT and MISMATCH files 
	echo -e "Get Reject\n\n"
	grep -v "^#" ${SOURCE_REJECT} | \
	grep 'NoTarget' | \
	awk -v OFS="\t" '{print $1,$2-1,$2,$3}' | \
	awk -v OFS="\t" '{tmp=$4 ; gsub(/_/, "\t", $4) ; print $1,$2,$3,tmp,$4}' | \
	awk -v OFS="\t" '{print $5,$6-1,$6,$4}' > ${SOURCE_EXTRACT}

	echo -e "Get MisMatch\n\n"
	grep -v "^#" ${SOURCE_REJECT} | \
	grep 'MismatchedRefAllele' | \
	awk -v OFS="\t" '{print $1,$2-1,$2,$3}' | \
	awk -v OFS="\t" '{tmp=$4 ; gsub(/_/, "\t", $4) ; print $1,$2,$3,tmp,$4}' | \
	awk -v OFS="\t" '{print $5,$6-1,$6,$4}' > ${SOURCE_MISMATCH}

	echo -e "Get ChrJump1\n\n"
	bcftools query -f "%CHROM\t%POS\t%ID\n" ${TARGET_OUT} | \
	awk -v OFS="\t" '{tmp=$3 ; gsub(/_/, "\t", $3) ; print $1,$2,tmp,$3}' | \
	awk -v OFS="\t" '{if($1 != $4) print $4,$5-1,$5,$3}' > ${TARGET_JUMP_CHR}

	cut -f4 ${TARGET_JUMP_CHR}  > ${TARGET_REMOVE}
	
	if [[ -f ${TARGET_REMOVE} && $( wc -l ${TARGET_REMOVE} | cut -f1 ) > 0 ]]
	then
		java -Djava.io.tmpdir=/home/shared/cathal/tmp -jar /home/shared/cathal/tools/gatk/GenomeAnalysisTK_3.8.jar \
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
	java -Djava.io.tmpdir=/home/shared/cathal/tmp -jar /home/shared/cathal/tools/picard.jar LiftoverVcf \
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
		java -Djava.io.tmpdir=/home/shared/cathal/tmp -jar /home/shared/cathal/tools/gatk/GenomeAnalysisTK_3.8.jar \
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
	grep -v "^#" ${TARGET_REJECT} | grep 'NoTarget' | \
	awk -v OFS="\t" '{print $1,$2-1,$2,$3}' | \
	awk -v OFS="\t" '{tmp=$4 ; gsub(/_/, "\t", $4) ; print $1,$2,$3,tmp,$4}' | \
	awk -v OFS="\t" '{print $5,$6-1,$6,$4}' > ${TARGET_EXTRACT}

	echo -e "Get Mismatch\n\n"
	grep -v "^#" ${TARGET_REJECT} | grep 'MismatchedRefAllele' | \
	awk -v OFS="\t" '{print $1,$2-1,$2,$3}' | \
	awk -v OFS="\t" '{tmp=$4 ; gsub(/_/, "\t", $4) ; print $1,$2,$3,tmp,$4}' | \
	awk -v OFS="\t" '{print $5,$6-1,$6,$4}' > ${TARGET_MISMATCH}


	echo -e "\nFinished with ${TARGET} to ${SOURCE}. \n"
	

done
