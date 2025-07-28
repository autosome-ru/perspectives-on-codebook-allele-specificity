import os
from pathlib import Path
from collections import Counter, defaultdict
from glob import glob
from sys import argv
from concurrent.futures import ProcessPoolExecutor as Executor

import pandas as pd
import numpy as np
from tqdm.auto import tqdm
from pqdm.processes import pqdm
from scipy.stats import fisher_exact, false_discovery_control
from scipy.stats import combine_pvalues

def chr_effect(row):
    FDR_less = row['fdr_comb_less']
    FDR_greater = row['fdr_comb_greater']
    if FDR_less < FDR_greater:
        nr = ['less', FDR_less]
    else:
        nr = ['greater', FDR_greater]
    return nr

def series_counter(series):
    counter = series.value_counts()
    dd = defaultdict(int)
    dd.update(counter.to_dict())
    return dd

def process_table(path):
    global HIT_THR
    df = pd.read_table(path)
    df['hit'] = df[['ref_motif_pval', 'alt_motif_pval']].min(axis=1) < HIT_THR
    df['x'] = 'No Hit'
    df['motif_conc'] = ~(df['hit']) * df['x'] + df['hit'] * df['motif_conc']
    df['signif'] = df['fdr_comb_pval'] < 0.05
    counter_signif = series_counter(df.query('signif')['motif_conc'])
    if counter_signif['Concordant'] + counter_signif['Discordant'] < 10:
        return
    df = df.query('fdr_comb_pval < 0.05 | fdr_comb_pval > 0.5')
    cont_table_hits = pd.crosstab(df['hit'], df['signif'])
    cont_table_hits.index = pd.Series(cont_table_hits.index).apply(lambda x: 1 if x else 0)
    cont_table_hits.columns = pd.Series(cont_table_hits.columns).apply(lambda x: 1 if x else 0)
    cont_table_hits = cont_table_hits.loc[[1, 0], [1, 0]]
    _, fisher_hits = fisher_exact(cont_table_hits, alternative='greater')
    if fisher_hits > 0.05:
        return
    
    counter_signif = series_counter(df.query('fdr_comb_pval < 0.05 & abs(motif_fc) > 1')['motif_conc'])
    counter_nonsignif = series_counter(df.query('fdr_comb_pval > 0.5 & abs(motif_fc) > 1')['motif_conc'])
    cont_table = [[counter_signif['Concordant'], counter_signif['Discordant']],
                  [counter_nonsignif['Concordant'], counter_nonsignif['Discordant']]]
    fisher_OR_less, fisher_PV_less = fisher_exact(cont_table, alternative='less')
    fisher_OR_greater, fisher_PV_greater = fisher_exact(cont_table, alternative='greater')
    tf = path.split('/')[-1].split('.')[0].split('@')[0]
    sample = '@'.join(path.split('/')[-1].split('.')[0].split('@')[1:])
    dct = {'TF': tf, 
           'sample': sample,
           'signif_NoHit': counter_signif['No Hit'],
           'signif_Concordant': counter_signif['Concordant'],
           'signif_Discordant': counter_signif['Discordant'],
           'nonsignif_NoHit': counter_nonsignif['No Hit'],
           'nonsignif_Concordant': counter_nonsignif['Concordant'],
           'nonsignif_Discordant': counter_nonsignif['Discordant'],
           'fisher_PV_less': fisher_PV_less,
           'fisher_PV_greater': fisher_PV_greater,
           'fisher_OR_less': fisher_OR_less,
           'fisher_OR_greater': fisher_OR_greater,
          }
    return dct


EXP_NAME = argv[1]
HIT_THR = 0.0005
DIR_TABLES = 'as_tables'
OUTPUT_PATH_DIR = Path(f'agg_{EXP_NAME}')
os.makedirs(OUTPUT_PATH_DIR, exist_ok=True)

exp = EXP_NAME
paths = glob(f'{DIR_TABLES}/*@{exp}@*.tsv')
cells = set(x.split('/')[-1].split('@')[2][:-4] for x in paths)
tables = []
for i, cell in enumerate(sorted(cells)):
    print(f'Processing {cell} ({i+1}/{len(cells)})!!')
    paths = glob(f'{DIR_TABLES}/*@{exp}@{cell}.tsv')
    table = pqdm(paths, process_table, n_jobs=250)
    table = filter(lambda x: x is not None, table)
    df_init = pd.DataFrame(table)
    df_init = pd.DataFrame(df_init).fillna(0)
    tables.append(df_init)
    print('-'*50)
df_init = pd.concat(tables)
df_init.to_csv(OUTPUT_PATH_DIR / f'pvalues_init.tsv', sep='\t', index=False)

for condition in ('less', 'greater'):
    df = df_init.copy()
    df['fisher_PV'] = df[f'fisher_PV_{condition}']
    df['fisher_OR'] = df[f'fisher_OR_{condition}']
    aggregated = df.groupby('TF').agg(
        median_OR=pd.NamedAgg(column="fisher_OR", aggfunc=lambda x: x.replace(np.inf, np.nan).median(skipna=True)),
        mean_OR=pd.NamedAgg(column="fisher_OR", aggfunc=lambda x: x.replace(np.inf, np.nan).mean(skipna=True)),
        agg_num=pd.NamedAgg(column="fisher_OR", aggfunc='count'),
        total_conc=pd.NamedAgg(column="signif_Concordant", aggfunc='sum'),
        total_disc=pd.NamedAgg(column="signif_Discordant", aggfunc='sum'),
        comb_pval=pd.NamedAgg(column="fisher_PV", aggfunc=lambda p: combine_pvalues(p).pvalue),
    )
    aggregated['fdr_comb'] = false_discovery_control(aggregated['comb_pval'])
    aggregated.sort_values(by='comb_pval', inplace=True)
    aggregated.to_csv(OUTPUT_PATH_DIR / f'pvalues_agg_{condition}.tsv', sep='\t')