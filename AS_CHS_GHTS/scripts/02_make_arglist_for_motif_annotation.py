from glob import glob
import os
from subprocess import run

with open('arglist.txt', 'w') as arglist:
    for file in glob('as_tables/*/ASB/*'):
        path = '/'.join(file.split('/')[:2])
        path = os.path.abspath(path)
        os.makedirs(f'{path}/ASB_motifs', exist_ok=True)
        tf = file.split('/')[-1].split('.')[0]
        pwm = f'/home/vladimirnoz/Projects/Codebook_Perspectives/motifs/manual_motifs_v2/pwms/{tf}.pwm'
        thr = f'/home/vladimirnoz/Projects/Codebook_Perspectives/motifs/manual_motifs_v2/uniform_thrs/{tf}.thr'
        snp_list = f'{path}/snp_lists'
        os.makedirs(snp_list, exist_ok=True)
        snp_list = f'{snp_list}/{tf}.txt'
        output = f'{path}/ASB_motifs/{tf}.tsv'
        print(os.path.abspath(file), pwm, thr, snp_list, output, 1, 0, file=arglist)
