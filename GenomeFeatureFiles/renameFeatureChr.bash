#!/bin/bash

f_repliOri="/data/git/2017_macalpine_chr_remod/GenomeFeatureFiles/repliOri_sacCer2_Eaton_ExcluderDNAtelomere_230sites.csv"
cp "/data/genome_feature_files/yeast/replication_origins/eaton_ORC-chip-seq/eaton_acs_filter-rDNA-telomere_230-sites_sacCer2.csv" $f_repliOri

cat $f_repliOri | sed -s "s|,chr\(\w\+\),|,\1,|" | \
	sed -s "s|,I,|,1,|" | \
	sed -s "s|,II,|,2,|" | \
	sed -s "s|,III,|,3,|" | \
	sed -s "s|,IV,|,4,|" | \
	sed -s "s|,V,|,5,|" | \
	sed -s "s|,VI,|,6,|" | \
	sed -s "s|,VII,|,7,|" | \
	sed -s "s|,VIII,|,8,|" | \
	sed -s "s|,IX,|,9,|" | \
	sed -s "s|,X,|,10,|" | \
	sed -s "s|,XI,|,11,|" | \
	sed -s "s|,XII,|,12,|" | \
	sed -s "s|,XIII,|,13,|" | \
	sed -s "s|,XIV,|,14,|" | \
	sed -s "s|,XV,|,15,|" | \
	sed -s "s|,XVI,|,16,|" \
	> .tmp.csv
mv -v .tmp.csv $f_repliOri
