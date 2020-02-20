inhib=MEK PIK
wtd = T
cutoff=0.5
R_3_6=../../R-3.6.0/bin/Rscript

networks:
	bsub  Rscript src/in-vitro-net.r 
	bsub  Rscript src/in-vivo-net.r

psites:
	for i in $(inhib); do \
		echo $(inhib) ; \
		echo $$i ; \
		bsub $(R_3_6) src/vivo_data.r	$$i ; \
	done

paths:
	bsub  Rscript src/find_path.r $(cutoff) $(wtd)

new_paths:
	bsub Rscript src/find_unknown_path.r  out/$(cutoff)-$(wtd)-unknown-limma-reg-sites-PIK3-new.tsv
