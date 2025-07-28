#!/home/vladimirnoz/miniconda3/bin/python
from glob import glob
from itertools import product

tables = glob('/home/vladimirnoz/Projects/Codebook_Perspectives/UDACHA/dnase/*.tsv') + glob('/home/vladimirnoz/Projects/Codebook_Perspectives/UDACHA/atac/*.tsv')
pwms = glob('/home/vladimirnoz/Projects/Codebook_Perspectives/motifs/manual_motifs_v2/pwms/*.pwm')


arglist = open('arglist.txt', 'w')
for table, pwm in product(tables, pwms):
    thr = pwm.replace('pwms', 'uniform_thrs').replace('.pwm', '.thr')
    tf = pwm.split('/')[-1].split('.')[0]
    exp = table.split('/')[-2]
    table_name = table.split('/')[-1].split('.')[0]
    perfectos = f'/home/vladimirnoz/Projects/Codebook_Perspectives/Chromatin/SNPs_lists/{tf}@{exp}@{table_name}.txt'
    output_path = f'/home/vladimirnoz/Projects/Codebook_Perspectives/Chromatin/as_tables/{tf}@{exp}@{table_name}.tsv'
    print(table, pwm, thr, perfectos, output_path, file=arglist)
arglist.close()
