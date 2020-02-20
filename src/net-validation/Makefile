inhib=MEK PIK

R_3_6=../../R-3.6.0/bin/Rscript

networks:
	bsub  Rscript in-vitro-net.r 
	bsub  Rscript in-vivo-net.r

psites:
	for i in $(inhib); do \
		echo $(inhib) ; \
		echo $$i ; \
		bsub $(R_3_6) vivo_data.r	$$i ; \
	done

paths:
	bsub  Rscript find_path.r $(cutoff) $(wtd)

new_paths:
	bsub Rscript find_unknown_path.r  out/$(cutoff)-$(wtd)-unknown-limma-reg-sites-PIK3-new.tsv
