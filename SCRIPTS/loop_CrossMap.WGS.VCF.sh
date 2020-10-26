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
	TO=10
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
	SOURCE_SORT=${PREFIX}_${SOURCE}_${i}.out.sort.vcf
	SOURCE_PASS=${PREFIX}_${SOURCE}_${i}.pass.vcf
	SOURCE_PASS_BED=${PREFIX}_${SOURCE}_${i}.pass.bed
	SOURCE_REJECT=${PREFIX}_${SOURCE}_${i}.reject.vcf
	SOURCE_EXTRACT=${PREFIX}_${SOURCE}_${i}.reject.extract.bed
	SOURCE_JUMP_CHR=${PREFIX}_${SOURCE}_${i}.CHR_jump.bed
	SOURCE_JUMP_POS=${PREFIX}_${SOURCE}_${i}.POS_jump.bed
	SOURCE_MISMATCH=${PREFIX}_${SOURCE}_${i}.mismatch.bed
	SOURCE_REMOVE=${PREFIX}_${SOURCE}_${i}.remove.txt

	TARGET_OUT=${PREFIX}_${TARGET}_${i}.out.vcf
	TARGET_SORT=${PREFIX}_${TARGET}_${i}.out.sort.vcf
	TARGET_PASS=${PREFIX}_${TARGET}_${i}.pass.vcf
	TARGET_PASS_BED=${PREFIX}_${TARGET}_${i}.pass.bed
	TARGET_REJECT=${PREFIX}_${TARGET}_${i}.reject.vcf
	TARGET_EXTRACT=${PREFIX}_${TARGET}_${i}.reject.extract.bed
	TARGET_JUMP_CHR=${PREFIX}_${TARGET}_${i}.CHR_jump.bed
	TARGET_MISMATCH=${PREFIX}_${TARGET}_${i}.mismatch.bed
	TARGET_REMOVE=${PREFIX}_${TARGET}_${i}.remove.txt



	# lift the VCF from SOURCE to TARGET
	echo -e "Conversion\n\n"
	CrossMap.py vcf ${CHAIN_SOURCE_TO_TARGET} ${SOURCE_PASS_PREV} ${REF_TARGET} ${TARGET_OUT}

	echo -e "\n\nStandardise reference allele, and sort\n\n"
	cat ${TARGET_OUT} | \
	awk '$1 ~ /^#/ {print $0; next} {if($4 == "A" || $4 == "C" || $4 == "G" || $4 == "T") {print $0 | "sort -k1,1 -k2,2n"}}' > ${TARGET_SORT}
	
	# extract variants with ambiguous reference alleles
	grep -v "^#" ${TARGET_OUT} | 
	awk '$4 != "A" && $4 != "C" && $4 != "G" && $4 != "T"' > ${TARGET_OUT}.refAllele

	# Convert the REJECT file to EXTRACT and MISMATCH files 
	echo -e "Get Reject\n\n"
	grep -v "^#" ${TARGET_OUT}.unmap | \
	grep 'Unmap' | \
	awk -v OFS="\t" '{print $1,$2-1,$2,$3}' | \
	awk -v OFS="\t" '{tmp=$4 ; gsub(/_/, "\t", $4) ; print $1,$2,$3,tmp,$4}' | \
	awk -v OFS="\t" '{print $5,$6-1,$6,$4}' | \
	bedtools sort -i - > ${SOURCE_EXTRACT}

	# include variants with ambiguous reference alleles
	echo -e "Get MisMatch\n\n"
	grep -v "^#" ${TARGET_OUT}.unmap | \
	grep "REF==ALT" | \
	cat ${TARGET_OUT}.refAllele - | \
	awk -v OFS="\t" '{print $1,$2-1,$2,$3}' | \
	awk -v OFS="\t" '{tmp=$4 ; gsub(/_/, "\t", $4) ; print $1,$2,$3,tmp,$4}' | \
	awk -v OFS="\t" '{print $5,$6-1,$6,$4}' | \
	bedtools sort -i - > ${SOURCE_MISMATCH}
	
	echo -e "Get ChrJump1\n\n"
	bcftools query -f "%CHROM\t%POS\t%ID\n" ${TARGET_SORT} | \
	awk -v OFS="\t" '{tmp=$3 ; gsub(/_/, "\t", $3) ; print $1,$2,tmp,$3}' | \
	awk -v OFS="\t" '{if($1 != $4) print $4,$5-1,$5,$3}' | \
	bedtools sort -i - > ${TARGET_JUMP_CHR}

	cut -f4 ${TARGET_JUMP_CHR}  > ${TARGET_REMOVE}
	
	if [[ -f ${TARGET_REMOVE} ]]
	then
		cat <( grep ^# ${TARGET_SORT} ) <(grep -v ^# ${TARGET_SORT} | ~/bin/grep -v -w -f ${TARGET_REMOVE}  ) > ${TARGET_PASS}
	else
		cp ${TARGET_SORT} ${TARGET_PASS}
	fi

	echo -e "\nFinished with ${SOURCE} to ${TARGET}. \n"



	# lift the VCF back from TARGET to SOURCE
	echo -e "Conversion\n\n"
	CrossMap.py vcf ${CHAIN_TARGET_TO_SOURCE} ${TARGET_PASS} ${REF_SOURCE} ${SOURCE_OUT}

	echo -e "\n\nSort\n\n"
	cat ${SOURCE_OUT} | \
	awk '$1 ~ /^#/ {print $0;next} {if($4 == "A" || $4 == "C" || $4 == "G" || $4 == "T") { print $0 | "sort -k1,1 -k2,2n"}}' > ${SOURCE_SORT}

	# extract variants with ambiguous reference alleles
	grep -v "^#" ${SOURCE_OUT} | \
	awk '$4 != "A" && $4 != "C" && $4 != "G" && $4 != "T"' > ${SOURCE_OUT}.refAllele
	
	# Get the PASS and JUMP file
	echo -e "\n\nGet ChrJump2 and PosJump\n\n"
	grep -v "^#" ${SOURCE_SORT} | \
	awk -v OFS="\t" '{tmp=$3 ; gsub(/_/, "\t", $3) ; print $1,$2-1,$2,tmp,$3}' | \
	awk -v OFS="\t" -v PASS=${SOURCE_PASS_BED} -v JUMPCHR=${SOURCE_JUMP_CHR} -v JUMPPOS=${SOURCE_JUMP_POS} '{ if($1 != $5) {print $5,$6-1,$6,$4 > JUMPCHR}  else if($3 != $6) {print $5,$6-1,$6,$4 > JUMPPOS} else {print $1,$2,$3,$4 > PASS} }'

	cat ${SOURCE_JUMP_CHR} ${SOURCE_JUMP_POS} | cut -f4 > ${SOURCE_REMOVE}


	echo -e "Get Pass \n\n"
	if [[ -f ${SOURCE_REMOVE} ]]
	then
		cat <( grep ^# ${SOURCE_SORT} ) <(grep -v ^# ${SOURCE_SORT} | ~/bin/grep -v -w -f ${SOURCE_REMOVE}  ) > ${SOURCE_PASS}
	else
		cp ${SOURCE_SORT} ${SOURCE_PASS}
	fi


	# Convert the REJECT to an EXTRACT 
	echo -e "Get Reject\n\n"
	grep -v "^#" ${SOURCE_OUT}.unmap | \
	grep 'Unmap' | \
	awk -v OFS="\t" '{print $1,$2-1,$2,$3}' | \
	awk -v OFS="\t" '{tmp=$4 ; gsub(/_/, "\t", $4) ; print $1,$2,$3,tmp,$4}' | \
	awk -v OFS="\t" '{print $5,$6-1,$6,$4}' > ${TARGET_EXTRACT}

	echo -e "Get Mismatch\n\n"
	grep -v "^#" ${SOURCE_OUT}.unmap | \
	grep "REF==ALT" | \
	cat ${SOURCE_OUT}.refAllele - | \
	awk -v OFS="\t" '{print $1,$2-1,$2,$3}' | \
	awk -v OFS="\t" '{tmp=$4 ; gsub(/_/, "\t", $4) ; print $1,$2,$3,tmp,$4}' | \
	awk -v OFS="\t" '{print $5,$6-1,$6,$4}' | \
	bedtools sort -i - > ${TARGET_MISMATCH}


	echo -e "\nFinished with ${TARGET} to ${SOURCE}. \n"
	

done
