arglist=$1

parallel --progress -j 100 --colsep=' ' -a $arglist python /home/vladimirnoz/Projects/Codebook_Perspectives/common/scripts/annotate_with_motifs.py {1} {2} {3} {4} {5} {6}
