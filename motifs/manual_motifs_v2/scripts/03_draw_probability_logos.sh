set -e

input=`ls motifs/*.txt`
for x in $input
do
    tf=`basename $x | cut -d '_' -f 1`
	out=logos_probs/$tf
    Rscript scripts/_plot_affimx_logo.R $x $out
done