
for i in 0.5 0.55 0.6 0.65 0.7 0.75
do
	bsub Rscript find_path.r $i T 
	bsub Rscript find_path.r $i T 

#	bsub Rscript find_path.r $i F


#	bsub Rscript find_path.r vivit $i T
 #   bsub Rscript find_path.r vivit $i F

done
