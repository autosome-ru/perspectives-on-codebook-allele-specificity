import pandas as pd
import numpy as np
from glob import glob
from tqdm.auto import tqdm
from collections import defaultdict
import os
tqdm.pandas()

def parse_phenotypes(path):
    table = pd.read_table(path, index_col='RSID',
                               dtype={x: str for x in ['ebi']})
    table = table.filter(items=['ebi'])
    return table

def parse_gtex(qtlfiles):

    result = {'cis': {}}
    print('Number of cis-eQTL files:', len(qtlfiles))

    for qtlfile in tqdm(list(qtlfiles)):
        with open(qtlfile) as qfile:

            tis = qtlfile[qtlfile.rfind('/') + 1:qtlfile.find('.v8.')]
            for line in qfile:

                if line.startswith('variant_id'):
                    tit = line.strip('\n').split('\t')
                    titlen = len(tit)
                    continue

                a = line.strip('\n').split('\t')
                a = {tit[x]: a[x] for x in range(titlen)}

                chrpos = '_'.join(a['variant_id'].split('_')[:2])
                result['cis'].setdefault(chrpos, (set(), set()))[0].add(tis)
                result['cis'][chrpos][1].add(a['gene_id'])

    return result

def annotate_with_gtex(row, qtl_info: dict, cis_trans: str, gene_or_pheno: int) -> str:
    if not row['posID'] in qtl_info[cis_trans]:
        return '-'
    return ';'.join(sorted(qtl_info[cis_trans][row['posID']][gene_or_pheno]))

def parse_cl_tf(path):
    tables = glob(f'{path}/*.tsv')
    to_concat = []
    for t in tqdm(tables):
        which = t.split('/')[-1][:-4]
        df = pd.read_table(t)
        df.query('(fdrp_bh_ref < 0.05) | (fdrp_bh_alt < 0.05)')
        df['which'] = which
        to_concat.append(df[['ID', 'which']])
    df = pd.concat(to_concat).sort_values(by='ID')
    return df.groupby('ID').agg(lambda x: ';'.join(x))
                                

def parse_adastra(path):
    by_tf = parse_cl_tf(f'{path}/TF').rename({'which': 'ADASTRA_TF'}, axis=1)
    by_cl = parse_cl_tf(f'{path}/CL').rename({'which': 'ADASTRA_CL'}, axis=1)
    return by_tf, by_cl
 


def annotate(snps_path, qtl_info, phenotypes_info, adastra_info):
    snps_positions = pd.read_table(snps_path)
    snps_positions['posID'] = snps_positions['#chr'] + '_' + snps_positions['end'].astype(str)
    annotation = snps_positions.join(phenotypes_info, how='left', on='id')
    
    annotation['eQTL_cis'] = annotation.apply(lambda x: annotate_with_gtex(x, qtl_info, 'cis', 0), axis=1)
    annotation['eQTL_cis_gene'] = annotation.apply(lambda x: annotate_with_gtex(x, qtl_info, 'cis', 1), axis=1)
    #annotation['ADASTRA_TFs'] = annotation.apply(lambda x: annotate_with_adastra(x, adastra_info), axis=1)
    annotation = annotation.merge(adastra_info[0], how='left', left_on='id', right_on='ID')
    annotation = annotation.merge(adastra_info[1], how='left', left_on='id', right_on='ID')
    
    return annotation.drop('posID', axis=1)

def annotate_proj(proj_path, out_path, qtl_info, phenotypes_info, adastra_info):
    #out_path = f'annotated'
    os.makedirs(out_path, exist_ok=True)
    snps_paths = glob(f'{proj_path}/*.tsv')
    for file in tqdm(snps_paths):
        annotation = annotate(file, qtl_info, phenotypes_info, adastra_info)
        filepath = file.split('/')[-1]
        annotation.to_csv(f'{out_path}/{filepath}', index=False, sep='\t', na_rep='-')
    
    
    
qtl_files = glob('/home/vladimirnoz/Projects/Codebook_Perspectives/common/Databases_2025/GTEx_Analysis_v8_eQTL/*signif*.txt')
qtl_info = parse_gtex(qtl_files)
phenotypes_info = parse_phenotypes('/home/vladimirnoz/Projects/Codebook_Perspectives/common/Databases_2025/ebi.tsv')
adastra_info = parse_adastra('/home/vladimirnoz/Projects/Codebook_Perspectives/common/Databases_2025/adastra_mabel')

#annotate_proj('as_tables/chipseq/ASB_motifs', 'as_tables/chipseq/ASB_phenotypes', qtl_info, phenotypes_info, adastra_info)
#annotate_proj('as_tables/selex/ASB_motifs', 'as_tables/selex/ASB_phenotypes', qtl_info, phenotypes_info, adastra_info)
annotate_proj('joint_tables/ASB', 'joint_tables/ASB_Phenotypes', qtl_info, phenotypes_info, adastra_info)


