import shlex
import subprocess
import io
from glob import glob
from sys import argv

import pandas as pd
import numpy as np

def get_name(row):
    return f"{row['#chr']}@{row['end']}@{row['ref']}@{row['alt']}"

def tf_to_bed(df_init, motif_length):
    df = df_init.copy()
    df['start'] = df['end'] - motif_length
    df['end'] = df['end'] + motif_length - 1
    df['name'] = df_init.apply(get_name, axis=1)
    return df[['#chr', 'start', 'end', 'name']]

def iterate_fasta(string):
    for line in string.split('\n'):
        if line.startswith('>'):
            name = line.strip()[1:]
        else:
            seq = line.strip()
            yield name, seq

def get_concordance(row):
    motif_ref, motif_alt = row['ref_motif_pval'], row['alt_motif_pval']
    motif_fc = row['motif_fc']
    as_ref, as_alt = -np.log10(row['ref_comb_pval']), -np.log10(row['alt_comb_pval'])
    as_fc = (as_alt - as_ref) / np.log10(2)
    hit = min(motif_ref, motif_alt) < 0.001
    if not hit:
        return 'No Hit'
    elif as_fc * motif_fc > 0:
        return 'Concordant'
    else:
        return 'Discordant'
    
INPUT_PATH = argv[1]
PWM_PATH = argv[2]
THR_PATH = argv[3]
PERFECTOS_PATH = argv[4]
OUTPUT_PATH = argv[5]
GENOME_PATH = '/sandbox/buyanchik/genome/GRCh38.primary_assembly.genome.fa'
APE_PATH = '/home/vladimirnoz/Projects/Codebook_Perspectives/common/scripts/ape-3.0.6.jar'
PVALUE_THR = 0.001


motif_length = len(open(PWM_PATH, 'r').readlines()) - 1
table = pd.read_table(INPUT_PATH)
table['name'] = table.apply(get_name, axis=1)

bed = tf_to_bed(table, motif_length)
bed_file = bed.to_csv(index=False, sep='\t', header=False)
bedtools_command = shlex.split(f'bedtools getfasta -fi {GENOME_PATH} -bed stdin -name')
fasta = subprocess.run(bedtools_command, input=bed_file, text=True, capture_output=True).stdout
snps_file = open(PERFECTOS_PATH, 'w')
for (name, seq), (index, row) in zip(iterate_fasta(fasta), table.iterrows()):
    pos = motif_length-1
    new_str = f"{name}\t{seq[:pos]}[{row['ref']}/{row['alt']}]{seq[pos+1:]}"
    print(new_str, file=snps_file)
snps_file.close()

perfectos_command = shlex.split(f'java -cp {APE_PATH} ru.autosome.perfectosape.SNPScan {PWM_PATH} {PERFECTOS_PATH} \
            --single-motif --precalc {THR_PATH} \
            -P 1 -F 0 --log-fold-change')
perfectos_proc = subprocess.run(perfectos_command, text=True, capture_output=True)
print(perfectos_proc.stderr)
perfectos = perfectos_proc.stdout
perfectos = io.StringIO(perfectos)
colnames = 'name tf ref_motif_pos ref_motif_orient ref_seq alt_motif_pos alt_motif_orient alt_seq alleles ref_motif_pval alt_motif_pval motif_fc'.split()
perfectos_table = pd.read_table(perfectos, names=colnames, skiprows=1)
conc_table = table.merge(perfectos_table, on='name')
conc_table['motif_conc'] = conc_table.apply(get_concordance, axis=1)
conc_table.drop(columns=['name', 'tf', 'alt_seq', 'ref_seq', 'alleles'], inplace=True)
conc_table.to_csv(OUTPUT_PATH, sep='\t', index=False)
