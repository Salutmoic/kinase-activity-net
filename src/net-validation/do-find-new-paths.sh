

for i in 0.5 0.55 0.6 0.65 0.7 0.75
do
	bsub Rscript find_unknown_path.r out/$i-T-unknown-rots-down-reg-sites.tsv
	bsub Rscript find_unknown_path.r out/$i-T-known-rots-down-reg-sites.tsv

#	bsub Rscript find_unknown_path.r out/$i-F-down-reg-sites.tsv


#	bsub Rscript find_unknown_path.r $i-vivit-T-down-reg-sites.tsv
#    bsub Rscript find_unknown_path.r $i-vivit-F-down-reg-sites.tsv

#	bsub Rscript find_unknown_path.r out/$i-T-rots-down-reg-sites.tsv

#    bsub Rscript find_unknown_path.r out/$i-F-rots-down-reg-sites.tsv

#   bsub Rscript find_unknown_path.r $i-vivit-T-rots-down-reg-sites.tsv
#    bsub Rscript find_unknown_path.r $i-vivit-F-rots-down-reg-sites.tsv


done
