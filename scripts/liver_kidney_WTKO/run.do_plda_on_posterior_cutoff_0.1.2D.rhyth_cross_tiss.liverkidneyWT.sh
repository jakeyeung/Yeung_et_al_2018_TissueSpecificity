#!/bin/sh
# Jake Yeung
# run.do_plda_on_posterior_cutoff_0.1.2D.rhyth_cross_tiss.liverkidneyWT.sh
# Run /home/yeung/projects/tissue-specificity/scripts/liver_kidney_WTKO/do_plda_on_posterior_cutoff_0.1.2D.rhyth_cross_tiss.liverkidneyWT.R
# 2016-06-28

runscript="/home/yeung/projects/tissue-specificity/scripts/liver_kidney_WTKO/do_plda_on_posterior_cutoff_0.1.2D.rhyth_cross_tiss.liverkidneyWT.R"

jcutoffhigh=3
jcutofflow=0
distfilt=10000

n=0
maxjobs=4
for distfilt in 5000 10000 20000 40000; do
	Rscript $runscript $distfilt $jcutoffhigh $jcutofflow&
	if (( $(($((++n)) % $maxjobs)) == 0 )) ; then
		# define maxjobs and n using maxjobsn skeleton
	    wait # wait until all have finished (not optimal, but most times good enough)
	    echo $n wait
	fi
done
wait
