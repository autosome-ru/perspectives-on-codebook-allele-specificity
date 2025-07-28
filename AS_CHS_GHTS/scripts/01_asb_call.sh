selex_tfs=`ls SNPs_filtered/selex`
mixalime create --bad-maps badmaps/selex.default.badmap.bed mixalimes/selex SNPs_filtered/selex/*
mixalime fit mixalimes/selex MCNB --small-dataset-n 1000 --window-size 1000 --n-jobs 20
mixalime test mixalimes/selex --n-jobs 20
for tf in $selex_tfs
do
    mixalime combine mixalimes/selex -g SNPs_filtered/selex/$tf --subname $tf
    mixalime export pvalues mixalimes/selex as_tables/selex/ASB/$tf.tsv --subname $tf --sample-info 
done
mixalime combine mixalimes/selex -g m:SNPs_filtered/selex/*/*.vcf.gz --subname all
mixalime export pvalues mixalimes/selex as_tables_combine_everything/selex.tsv --subname all --sample-info


chipseq_tfs=`ls SNPs_filtered/chipseq`
mixalime create --bad-maps badmaps/chipseq.default.badmap.bed mixalimes/chipseq SNPs_filtered/chipseq/*
mixalime fit mixalimes/chipseq MCNB --small-dataset-n 10000 --window-size 10000 --n-jobs 20
mixalime test mixalimes/chipseq --n-jobs 20
for tf in $chipseq_tfs
do
    mixalime combine mixalimes/chipseq -g SNPs_filtered/chipseq/$tf --subname $tf
    mixalime export pvalues mixalimes/chipseq as_tables/chipseq/ASB/$tf.tsv --subname $tf --sample-info
done
mixalime combine mixalimes/chipseq -g m:SNPs_filtered/chipseq/*/*.vcf.gz --subname all
mixalime export pvalues mixalimes/chipseq as_tables_combine_everything/chipseq.tsv --subname all --sample-info
