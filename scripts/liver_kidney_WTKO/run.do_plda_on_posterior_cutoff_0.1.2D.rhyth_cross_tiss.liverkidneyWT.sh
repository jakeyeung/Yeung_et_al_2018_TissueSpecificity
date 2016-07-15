#!/bin/sh
# Jake Yeung
# run.do_plda_on_posterior_cutoff_0.1.2D.rhyth_cross_tiss.liverkidneyWT.sh
# Run /home/yeung/projects/tissue-specificity/scripts/liver_kidney_WTKO/do_plda_on_posterior_cutoff_0.1.2D.rhyth_cross_tiss.liverkidneyWT.R
# 2016-06-28

runscript="/home/yeung/projects/tissue-specificity/scripts/liver_kidney_WTKO/do_plda_on_posterior_cutoff_0.1.2D.rhyth_cross_tiss.liverkidneyWT.R"

jcutoffhigh=3
jcutofflow="0"
method="g=1001"
distfilt=20000
model="Kidney_SV129"

n=0
maxjobs=8
for method in "g=1001"; do
	for jcutoffhigh in "1.5" "2"; do
		for jcutofflow in "0" "0.5"; do
			Rscript $runscript $distfilt $jcutoffhigh $jcutofflow $method $model&
			if (( $(($((++n)) % $maxjobs)) == 0 )) ; then
				# define maxjobs and n using maxjobsn skeleton
			    wait # wait until all have finished (not optimal, but most times good enough)
			    echo $n wait
			fi
		done
	done
done
wait

# for distfilt in 40000; do
# done
