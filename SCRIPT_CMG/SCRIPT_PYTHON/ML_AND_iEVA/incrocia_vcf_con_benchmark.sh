#!bin/bash


cat $1 | while read line
	do
		python /home/jarvis/git/CMG/SCRIPT_CMG/SCRIPT_PYTHON/ML_AND_iEVA/incrocia_benchmark_con_vcf_IEVA.py \
		--vcf $line \
		--var_list /home/jarvis/Scrivania/TEST/bam-ieva/benchmark/NUOVO_BENCHMARK_DEF.tsv \
		--window 50 \
		--out ${line%.*.*}.bmk.vcf
	done
	