
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::     ANNOTATION     ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

ANNOTATION () {

	printf "\n\n"
	cat $LOGHI/logo_annotation.txt
	printf $"\n\n\n"
	
	if [[ "$2" == "CANONICAL" ]]
		then
		perl $VEPANN -i $1 \
		-o ${1%.*}.ANN.nosort.vcf \
		--stats_file ${1%.*}.ANN.html \
		--cache \
		--assembly GRCh37 \
		--offline \
		--force_overwrite \
		-v \
		--fork 2 \
		--variant_class \
		--sift b \
		--poly b \
		--vcf_info_field ANN \
		--hgvs \
		--protein \
		--canonical \
		--check_existing \
		--gmaf \
		--pubmed \
		--species homo_sapiens \
		--failed 1 \
		--vcf \
		--pick

	else

		perl $VEPANN -i $1 \
		-o ${1%.*}.ANN.nosort.vcf \
		--stats_file ${1%.*}.ANN.html \
		--cache \
		--assembly GRCh37 \
		--offline \
		--force_overwrite \
		-v \
		--fork 2 \
		--variant_class \
		--sift b \
		--poly b \
		--vcf_info_field ANN \
		--hgvs \
		--protein \
		--canonical \
		--check_existing \
		--gmaf \
		--pubmed \
		--species homo_sapiens \
		--failed 1 \
		--vcf
	fi

	vcf-sort -c ${1%.*}.ANN.nosort.vcf > ${1%.*}.ANN.vcf
	mv $1 $DELETE
	mv ${1%.*}.ANN.nosort.vcf $DELETE
	INPUT=${1%.*}.ANN.vcf
	printf $'\n =========>	ANNOTATION: DONE\n'

}


ADD_ANNOTATION() {

	printf $'\n =========>	ADD ANNOTATION\n'
	
	python $SCRIPT_PIPELINE/Estrai_Annotazione.py \
  	--vcf $1 \
  	--tsv $2 \
  	--tag_list $ANN_LIST \
  	-o  ${2%.*}.ANN.tsv
 	

 	mv $2 $DELETE
 	
 	sort -V ${2%.*}.ANN.tsv > $2
  	 	
 	cp $2 $STORAGE
 	mv $2 $OUT

  	mv ${2%.*}.ANN.tsv $DELETE
  	printf $'\n =========>	ADD ANNOTATION: DONE\n'
}


SPLIT_TRANSCRIPTS () {

	printf $'\n =========>	SPLIT TRANSCRIPTS\n'

	python $SCRIPT_PIPELINE/Estrai_Annotazione.py \
  	--vcf $1 \
  	--trs_list $2 \
  	--split \
  	-o  ${1%.*.*}.Trans

  	INPUT1=${1%.*.*}.Trans.vcf
  	INPUT2=${1%.*.*}.Trans.other.vcf

  	printf $'\n =========>	SPLIT TRANSCRIPTS: DONE\n'

}

MERGE_2VCF () {

	printf $'\n =========>	MERGING VCF\n'

	bgzip $1 && tabix $1.gz
	bgzip $2 && tabix $2.gz

	vcf-concat $1.gz $2.gz > ${1%.*.*}.Merge.vcf
	vcf-sort -c ${1%.*.*}.Merge.vcf > ${1%.*.*}.vcf

	mv $1.gz $DELETE
	mv $2.gz $DELETE
	mv $1.gz.tbi $DELETE
	mv $2.gz.tbi $DELETE
	mv ${1%.*.*}.Merge.vcf $DELETE

	INPUT=${1%.*.*}.vcf

	printf $'\n =========>	MERGING VCF: DONE\n'

}

# ANNOTATION_somatic () {

# 	printf "\n\n"
# 	cat $LOGHI/logo_annotation.txt
# 	printf $"\n\n\n"
	
# 	perl $VEPANN -i $1 \
#  	-o ${1%.*}.ANN.vcf \
#  	--stats_file ${1%.*}.ANN.html \
#  	--cache \
#  	--assembly GRCh37 \
#  	--offline \
#  	--force_overwrite \
#  	-v \
#  	--fork 10 \
#  	--variant_class \
#  	--sift b \
#  	--poly b \
#  	--vcf_info_field ANN \
#  	--hgvs \
#  	--protein \
#  	--canonical \
#  	--check_existing \
#  	--gmaf \
#  	--pubmed \
#  	--species homo_sapiens \
#  	--failed 1 \
#  	--vcf

#  	mv $1 $DELETE
#  	INPUT=${1%.*}.ANN.vcf
#  	printf $'\n =========>	ANNOTATION: DONE\n'
# }


ADD_ANNOTATION_somatic () {
	
	python $SCRIPT_PIPELINE/Estrai_Annotazione_Somatic.py \
  	-i $1 \
  	-f ${1%.*.*}.features.tsv \
  	-l $ANN_LIST_SOMATIC \
 	-t $TRANSCR_LIST \
  	-o ${1%.*.*}.ANN
 	
 	sort -V ${1%.*.*}.ANN.tsv > ${1%.*.*}.ANN.sort.tsv
  	sort -V ${1%.*.*}.ANN.Other_transcripts.tsv > ${1%.*.*}.ANN.Other_transcripts.sort.tsv
 	
 	cp ${1%.*.*}.ANN.sort.tsv $STORAGE
 	cp ${1%.*.*}.ANN.Other_transcripts.sort.tsv $STORAGE
 	mv ${1%.*.*}.ANN.sort.tsv $OUT
 	mv ${1%.*.*}.ANN.Other_transcripts.sort.tsv $OUT
 	
  	mv ${2%.*}.ANN.tsv $DELETE
  	mv ${2%.*}.ANN.Other_transcripts.tsv $DELETE
}