#!/usr/bin/env python3

import os
import pandas as pd
import numpy as np
import h5py
from random import shuffle
from multiprocessing import Pool

AA_LIST = ['G','P','A','V','L','I','M','C','F','Y','W','H','K','R','Q','N','E','D','S','T']
POSITIONS_LIST = list(range(-5,5))

def scale_neg1_pos1(
    
    series:pd.Series

):

    series = 2 * series / (series.max() - series.min())
    series = series - series.min() - 1

    return series

def scale_01(
    
    series:pd.Series

):

    series = series / (series.max() - series.min())
    series = series - series.min()

    return series

def score_sequence(
    
    seq_str:str,
    pssm_df:pd.DataFrame

):

    if len(seq_str) != len(pssm_df):

        raise Exception('Sequence length cannot be different from pssm length')

    """ multiply probabilities instead of adding scores
    p = 1

    for i in range(len(seq_str)):

        pos = list(pssm_df.index)[i]
        p = p * pssm_df.loc[pos, seq_str[i]]
    """
    p = 0

    if pssm_df.loc[0, seq_str[5]] != 0:

        for i in range(len(seq_str)):

            if seq_str[i] != '_':

                pos = list(pssm_df.index)[i]
                p = p + pssm_df.loc[pos, seq_str[i]]

    return p

def pssm_scoring(
    
    seq: str,
    pssm_df_dict: dict

):
    
    record = {}
    
    for kinase in pssm_df_dict.keys():
        
        p = score_sequence(seq, pssm_df_dict[kinase])
        record[kinase] = p
    
    out_series = pd.Series(record)
    
    return out_series

# take the largest deviation between test and null curves
def ks_statistic(
    
    test_var_series:pd.Series,
    null_var_df:pd.DataFrame,
    plot_bool:bool=False,
    out_plot_dir_str:str=''

):

    kinase = test_var_series.name
    
    deltas = []
    running_sum_test = 0
    running_sum_null = 0
    
    for i in range(len(test_var_series)):
    
        running_sum_test = running_sum_test + test_var_series[i]
        running_sum_null = running_sum_null + null_var_df[kinase][i]
        deltas.append(running_sum_test - running_sum_null)
    
    max_delta = max(deltas)
    min_delta = min(deltas)
    
    if plot_bool:
    
        data = pd.Series(deltas)
    
        import matplotlib.pyplot as plt
        import seaborn as sns
        
        if not os.path.exists(out_plot_dir_str):
            os.makedirs(out_plot_dir_str)
    
        plt.clf()
        fig, ax = plt.subplots(figsize=(4, 3))
    
        ax = sns.lineplot(
            data=data,
            linewidth=2.5
        )
    
        ax.set(
            title=f'{kinase}',
            xlabel='Phosphosite rank',
            ylabel='Running Sum Statistic'
        )
    
        sns.despine()
        plt.tight_layout()
    
        plt.savefig(os.path.join(out_plot_dir_str, f'{kinase}.svg'))
        plt.close()
    
    if abs(max_delta) > abs(min_delta):
    
        return max_delta
    
    else:
    
        return min_delta
""" # take the deviation at the end of the test and null curves
def ks_statistic(
    
    test_var_series:pd.Series,
    null_var_df:pd.DataFrame,
    plot_bool:bool=False,
    out_plot_dir_str:str=''

):
    
    kinase = test_var_series.name
    
    running_sum_test = 0
    running_sum_null = 0
    
    running_sum_test_list = [0]
    running_sum_null_list = [0]
    
    for i in range(len(test_var_series)):
        
        running_sum_test = running_sum_test + test_var_series[i]
        running_sum_null = running_sum_null + null_var_df[kinase][i]
        
        running_sum_test_list.append(running_sum_test)
        running_sum_null_list.append(running_sum_null)
   
    if plot_bool:
        
        import matplotlib.pyplot as plt
        import seaborn as sns
        
        if not os.path.exists(out_plot_dir_str):
            os.makedirs(out_plot_dir_str)
        
        data = pd.DataFrame({
            f'{kinase}': running_sum_test_list,
            'Null': running_sum_null_list
        })
        
        plt.clf()
        fig, ax = plt.subplots(figsize=(4, 3))
        
        ax = sns.lineplot(
            data=pd.DataFrame(data),
            linewidth=2.5
        )
        
        ax.set(
            title=f'{kinase}',
            xlabel='Phosphosite rank',
            ylabel='Statistic'
        )
        
        sns.despine()
        plt.tight_layout()
        
        plt.savefig(os.path.join(out_plot_dir_str, f'{kinase}.svg'))
        plt.close()
    
    delta = running_sum_test - running_sum_null
    
    return delta
"""

def compute_null_ks(
    
    pssm_scoring_df:pd.DataFrame,
    seqrnk_null_scores_df:pd.DataFrame,
    seqrnk_series:pd.Series,

):

    idx_list = list(range(len(seqrnk_series)))

    # permute phosphosites
    shuffle(idx_list)
    shuffled_pssm_scoring_df = pssm_scoring_df.reindex(idx_list)
    shuffled_pssm_scoring_df.index = list(range(len(seqrnk_series)))

    # multiply the kinase pssm scores by the phosphosites seqrnk
    shuffled_seqrnk_pssm_scores_df = shuffled_pssm_scoring_df.apply(
        lambda x: x * seqrnk_series,
        axis=0
    )

    # compute ks between shuffled pssm scores and null scores
    ks_series = shuffled_seqrnk_pssm_scores_df.apply(
        ks_statistic,
        args=[seqrnk_null_scores_df, False],
        axis=0
    )

    ks_df = pd.DataFrame([ks_series])

    return ks_df

def compute_ks_empirical_distrib(
    
    pssm_scoring_df:pd.DataFrame,
    null_scoring_df:pd.DataFrame,
    seqrnk_series:pd.Series,
    n:int=1000,
    n_proc:int=1

):
    # multiply the null pssm scores by the phosphosites seqrnk
    seqrnk_null_scores_df = null_scoring_df.apply(
        lambda x: x * seqrnk_series,
        axis=0
    )
    
    arg1 = [pssm_scoring_df for i in range(n)] 
    arg2 = [seqrnk_null_scores_df for i in range(n)]
    arg3 = [seqrnk_series for i in range(n)]
    
    with Pool(processes=n_proc) as pool:
        dfs_list = pool.starmap(compute_null_ks, zip(arg1, arg2, arg3))
    
    out_df = pd.concat(dfs_list, axis=0)
    out_df.index = list(range(len(out_df)))
    
    return out_df

def compute_ks_pvalues(
    
    ks_empirical_distrib_df: pd.DataFrame,
    ks_series: pd.Series,
    plot_bool:bool=False,
    out_plot_dir_str:str=''

):

    if plot_bool:

        import matplotlib.pyplot as plt
        import seaborn as sns
        
        if not os.path.exists(out_plot_dir_str):
            os.makedirs(out_plot_dir_str)
        
    ks_pvalue_series = pd.Series()
    ks_pvalue_series.name = 'p value'
    
    for kinase_str in ks_empirical_distrib_df.columns:
    
        empirical_list = ks_empirical_distrib_df[kinase_str]
    
        value_float = ks_series[kinase_str]
    
        num_int = sum([x>value_float for x in empirical_list])
        den_int = len(empirical_list)
        pvalue = num_int / den_int
    
        if pvalue > 0.5:
            num_int = den_int - num_int
            pvalue = num_int / den_int
        ks_pvalue_series[kinase_str] = pvalue
    
        if plot_bool:
            
            plt.clf()
            fig, ax = plt.subplots(figsize=(5, 3))
            
            ax = sns.displot(
                data=empirical_list,
                height=3,
                aspect=5/3
            )
            
            ax.set(
                title=f'Empirical ks distribution for {kinase_str} ({len(empirical_list)} perm.)',
                xlabel='KS',
                ylabel='Count'
            )
            
            plt.axvline(x=value_float)
            
            sns.despine()
            plt.tight_layout()
            
            plt.savefig(os.path.join(out_plot_dir_str, f'{kinase_str}_ks_dist.svg'))
            plt.close()
    
    return(ks_pvalue_series)

def read_pssms(
    
    pssms_h5_file:str

):
    
    pssms_h5 = h5py.File(pssms_h5_file, 'r')
    pssm_df_dict = {}
    
    for kinase in pssms_h5.keys():
        
        pssm_df_dict[kinase] = pd.DataFrame(pssms_h5[kinase])
        pssm_df_dict[kinase].columns = AA_LIST
        pssm_df_dict[kinase].index = POSITIONS_LIST
    
    return pssm_df_dict

def read_seqrnk(
    
    seqrnk_file:str

):

    seqrnk = pd.read_csv(seqrnk_file, sep='\t', header=None)
    seqrnk.columns = ['Sequence', 'Score']

    return seqrnk

def kinase_activities(
    
    seqrnk_file:str,
    pssms_h5_file:str,
    n_perm:int=10000,
    n_proc:int=1,
    plot_figures:bool=False,
    out_plot_dir:str='',

):
    
    pssm_df_dict = read_pssms(pssms_h5_file)
    seqrnk = read_seqrnk(seqrnk_file)
    
    # tmp: filter out Y phosphosites
    seqrnk = seqrnk.loc[
        [seqrnk['Sequence'][i][5] == 'S' or seqrnk['Sequence'][i][5] == 'T' for i in range(len(seqrnk))]
    ]
    seqrnk.index = range(len(seqrnk))

    # score phosphosite sequences with each PSSM
    seq_series = seqrnk['Sequence']
    arg1 = [seq_series[i] for i in range(len(seq_series))]
    arg2 = [pssm_df_dict for i in range(len(seq_series))]
    with Pool(processes=n_proc) as pool:
        dfs_list = pool.starmap(pssm_scoring, zip(arg1, arg2))
    pssm_scoring_df = pd.concat(dfs_list, axis=1).T
    pssm_scoring_df.index = list(range(len(seq_series)))

    # scale the PSSM scores between 0 and 1 for each phosphosite
    pssm_scoring_scaled01_df = pssm_scoring_df.apply(
        scale_01,
        axis=1
    )

    # obtain mean (null) pssm scores for each psssm
    # TODO can be made just a Series (following code to be optimised)
    pssm_null_scoring_df = pssm_scoring_scaled01_df.apply(
        lambda x: pd.Series([x.mean()] * len(x), index=x.index),
        axis=0
    )

    # compute empirical distribution of ks statistics
    ks_empirical_distrib_df = compute_ks_empirical_distrib(
        pssm_scoring_df=pssm_scoring_scaled01_df,
        null_scoring_df=pssm_null_scoring_df,
        seqrnk_series=seqrnk['Score'],
        n=n_perm,
        n_proc=n_proc
    )

    # multiply the pssm scores by the phosphosites seqrnk
    seqrnk_pssm_scoring_row_scaled_df = pssm_scoring_scaled01_df.apply(
        lambda x: x * seqrnk['Score'],
        axis=0
    )

    # multiply the null scores by the phosphosites seqrnk
    seqrnk_pssm_null_scoring_df = pssm_null_scoring_df.apply(
        lambda x: x * seqrnk['Score'],
        axis=0
    )

    # compute ks statistics
    ks_series = seqrnk_pssm_scoring_row_scaled_df.apply(
        ks_statistic,
        args=[seqrnk_pssm_null_scoring_df, plot_figures, out_plot_dir],
        axis=0
    )
    ks_series.name = 'KS'

    # compute ks pvalues
    ks_pvalue_series = compute_ks_pvalues(
        ks_empirical_distrib_df=ks_empirical_distrib_df,
        ks_series=ks_series,
        plot_bool=plot_figures,
        out_plot_dir_str=out_plot_dir
    )

    # compute FPR
    results_df = pd.concat([ks_series, ks_pvalue_series], axis=1)
    results_df['FPR p value'] = results_df['p value'] * len(results_df)
    results_df['FPR p value'].loc[results_df['FPR p value'] > 1] = 1

    # compute activity score as -log10(p value), signed with the sign of ks,
    # capped based on smallest computable p value
    activity_score_series = results_df['p value'].copy()
    min_non0_fpr_float = min(activity_score_series.loc[activity_score_series != 0])
    max_abs_activity_score_float = -np.log10(min_non0_fpr_float)
    activity_score_series.loc[activity_score_series == 0] = np.nextafter(np.float32(0), np.float32(1))
    activity_score_series = -np.log10(activity_score_series)
    activity_score_series.loc[activity_score_series > max_abs_activity_score_float] = max_abs_activity_score_float
    results_df['Activity Score'] = activity_score_series * results_df['KS'].apply(lambda x: 1 if x > 0 else -1)
    
    # export results
    print(results_df.to_csv(sep='\t', header=True, index=True))
