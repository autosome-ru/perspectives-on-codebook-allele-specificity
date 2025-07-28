sarus_path=/home/vladimirnoz/Projects/Codebook_Perspectives/Chromatin/scripts/sarus.jar
fasta_path=$1
pwm_path=$2
thresholds_path=$3

java -cp "${sarus_path}" ru.autosome.SARUS "${fasta_path}" \
                            "${pwm_path}" \
                            -10000000 \
                            --pvalues-file "${thresholds_path}" \
                            --threshold-mode score \
                            --output-scoring-mode pvalue 