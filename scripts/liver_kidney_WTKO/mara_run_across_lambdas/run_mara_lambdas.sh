#!/bin/sh
# Jake Yeung
# run_mara_lambdas.sh
# Run MARA for clock-driven tissue-wide module across lambdas
# to see how activities behave across global lambdas
# 2016-09-02

# default lambda: 0.0479929486250763 (-1.318823 in log10)
# use different lambdas at different log scales
# run this in R to get your lambdas and copy to shellscript (easier)
# x <- -1.318823; xvec <- c(seq(x - 5, by = 1, length.out = 5), seq(x, by = 1, length.out = 5)); paste(10^xvec, collapse = " ")
# output is
# lambdas="-6.318823 -5.318823 -4.318823 -3.318823 -2.318823 -1.318823 -0.318823 0.681177 1.681177 2.681177"
# lambdas="4.79929007479202e-07 4.79929007479202e-06 4.79929007479202e-05 0.000479929007479202 0.00479929007479202 0.0479929007479202 0.479929007479202 4.79929007479202 47.9929007479202 479.929007479202"

# in linear
# paste(seq(0.0479929 - 0.02, by = 0.005, length.out = 8), collapse = " ")
lambdas="0.0279929 0.0329929 0.0379929 0.0429929 0.0479929 0.0529929 0.0579929 0.0629929"

inscript="/home/yeung/projects/tissue-specificity/scripts/liver_kidney_WTKO/liver_kidney_new_modules.R"

n=0
maxjobs=10
for lambda in $lambdas; do
	Rscript $inscript $lambda&
	# exit 0
	# echo $lambda	
	if (( $(($((++n)) % $maxjobs)) == 0 )) ; then
		# define maxjobs and n using maxjobsn skeleton
	    wait # wait until all have finished (not optimal, but most times good enough)
	    echo $n wait
	fi
done


