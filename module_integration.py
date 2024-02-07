'''
Run module_inference_script.py first. This script accepts the inferred module dynamics.
'''

import matplotlib.pyplot as plt
from numpy import arange
from pandas import DataFrame, read_csv
from scipy.stats import pearsonr, rankdata
from seaborn import heatmap

######################################################################

def main():
    # read in phospho module data
    print('importing phospho modules ... ', end='')
    df_phos_modules = read_csv('phospho_module_dynamics.csv')
    df_phos_modules = df_phos_modules.set_index(df_phos_modules.columns[0])
    n_phos_modules = len(df_phos_modules.index)
    df_phos_error = read_csv('phospho_module_error.csv')
    df_phos_error = df_phos_error.set_index(df_phos_error.columns[0])
    print('done')
    
    # plot each module's dynamics
    print('plotting phospho module dynamics ... ', end='')
    for i,row in df_phos_modules.iterrows():
        x = arange(len(row))
        y = row
        fig = plt.figure(figsize=(7,2))
        plt.plot(x, y, c='firebrick', zorder=2)
        error = df_phos_error.loc[i]
        plt.fill_between(x, y-error, y+error, alpha=0.1, facecolor='firebrick', edgecolor='firebrick', zorder=2)
        plt.xticks(ticks=x, labels=[k.split(' Abundance')[0] for k in row.index], rotation=30, ha='right', va='top')
        x_min, x_max = -0.2, len(x)-0.8
        plt.xlim(x_min,x_max)
        plt.xlabel('Time')
        plt.ylabel('Mean phosphosite\nabundance (z-score)')
        if n_phos_modules > 26: m = i[0] + r'$_{' + i[1:] + '}$'
        else: m = i
        plt.title('Signaling module ' + m + '\n')
        plt.axvline(x=7, color='k', linestyle='--', linewidth=1, zorder=1)
        xaxis_length = x_max-x_min
        plt.annotate('\nVEM on', xy=(0.2/xaxis_length, 1.05), xycoords='axes fraction', ha='left', va='center', fontsize=10, color='k', weight='bold')
        plt.annotate('\nVEM off', xy=((7+0.4)/xaxis_length, 1.05), xycoords='axes fraction', ha='left', va='center', fontsize=10, color='k', weight='bold')
        plt.gca().spines['right'].set_visible(False)
        plt.gca().spines['top'].set_visible(False)
        plt.savefig('dynamics_phos_' + i + '.png', dpi=400, bbox_inches='tight')
        plt.close()
    print('done\n')
    
    ####################################################################################################
    
    # read in RNA module data
    print('importing RNA modules ... ', end='')
    df_RNA_modules = read_csv('RNA_module_dynamics.csv')
    df_RNA_modules = df_RNA_modules.set_index(df_RNA_modules.columns[0])
    n_RNA_modules = len(df_RNA_modules.index)
    df_RNA_error = read_csv('RNA_module_error.csv')
    df_RNA_error = df_RNA_error.set_index(df_RNA_error.columns[0])
    print('done')
    
    # plot each module's dynamics
    print('plotting RNA module dynamics ... ', end='')
    for i,row in df_RNA_modules.iterrows():
        x = arange(len(row))
        y = row
        fig = plt.figure(figsize=(7,2))
        plt.plot(x, y, c='mediumblue', zorder=2)
        error = df_RNA_error.loc[i]
        plt.fill_between(x, y-error, y+error, alpha=0.1, facecolor='mediumblue', edgecolor='mediumblue', zorder=2)
        plt.xticks(ticks=x, labels=[k.split(' Abundance')[0] for k in row.index], rotation=30, ha='right', va='top')
        x_min, x_max = -0.2, len(x)-0.8
        plt.xlim(x_min,x_max)
        plt.xlabel('Time')
        plt.ylabel('Mean transcript\nabundance (z-score)')
        if n_RNA_modules > 26: m = i[0] + r'$_{' + i[1:] + '}$'
        else: m = i
        plt.title('Transcriptional module ' + m + '\n')
        plt.axvline(x=5, color='k', linestyle='--', linewidth=1, zorder=1)
        xaxis_length = x_max-x_min
        plt.annotate('\nVEM on', xy=(0.2/xaxis_length, 1.05), xycoords='axes fraction', ha='left', va='center', fontsize=10, color='k', weight='bold')
        plt.annotate('\nVEM off', xy=((5+0.4)/xaxis_length, 1.05), xycoords='axes fraction', ha='left', va='center', fontsize=10, color='k', weight='bold')
        plt.gca().spines['right'].set_visible(False)
        plt.gca().spines['top'].set_visible(False)
        plt.savefig('dynamics_RNA_' + i + '.png', dpi=400, bbox_inches='tight')
        plt.close()
    print('done\n')
    
    ###################################################################################################################
    
    # since the phospho and RNA-seq data don't have all the same timepoints, get
    # the ones in common between the two for integration
    quan_fields_shared = set(df_phos_modules.columns).intersection(set(df_RNA_modules.columns))
    quan_fields_shared = [c for c in df_phos_modules if c in quan_fields_shared]
    
    # integrate phospho and RNA modules by calculating Pearson correlation between every pair
    print('performing module integration ... ', end='')
    phospho_col, RNA_col, r_col, pval_col = [], [], [], []
    df_sync_hm = DataFrame(index=df_phos_modules.index, columns=df_RNA_modules.index)
    for i,row_RNA in df_RNA_modules.iterrows():
        # make columns for the heatmap (hm) dataframes
        sync_col_hm = []
        for j,row_phos in df_phos_modules.iterrows():
            r, pval = pearsonr(row_phos[quan_fields_shared], row_RNA[quan_fields_shared])
            phospho_col.append(j)
            RNA_col.append(i)
            r_col.append(r)
            pval_col.append(pval)
            sync_col_hm.append(r)
        df_sync_hm[i] = sync_col_hm
    print('done')
    
    # B-H correct p-values
    pval_rank = rankdata(pval_col)
    n_pvals = len(pval_col)
    FDR_col = [0]*n_pvals
    for i in range(n_pvals):
        FDR = pval_col[i]*n_pvals/pval_rank[i]
        if FDR > 1: FDR = 1
        FDR_col[i] = FDR
    
    # export dataframe of correlations and pvalues
    print('calculating statistics, plotting heatmap ... ', end='')
    df_stats = DataFrame()
    df_stats['Phospho Module'] = phospho_col
    df_stats['RNA Module'] = RNA_col
    df_stats['Pearson Correlation'] = r_col
    df_stats['P-value'] = pval_col
    df_stats['Corrected p-value (FDR)'] = FDR_col
    df_stats = df_stats.sort_values('Pearson Correlation', ascending=False)
    try: df_stats.to_csv('integration_stats.csv', index=False)
    except PermissionError: print('integration_stats.csv is open and therefore was not updated')
    
    # plot heatmap of scores (correlations)
    max_sync_score = df_sync_hm.max().max()
    fig = plt.figure(figsize=(6,4))
    g = heatmap(df_sync_hm, cmap='seismic', vmin=-1, vmax=1, xticklabels=True, yticklabels=True,
                linewidths=1, linecolor='k', cbar_kws={'label': r"Pearson's $r$", 'ticks': [-1,0,1], 'shrink': 0.3})
    g.tick_params(left=False, bottom=False)
    plt.xlabel('Transcriptional module')
    plt.ylabel('Signaling module')
    if n_RNA_modules > 26: xlabels = [c[0] + r'$_{' + c[1:] + '}$' for c in df_sync_hm.columns]
    else: xlabels = df_sync_hm.columns
    g.set_xticklabels(xlabels, rotation=0, fontsize=8, ha='center', va='center')
    if n_phos_modules > 26: ylabels = [c[0] + r'$_{' + c[1:] + '}$' for c in df_sync_hm.index]
    else: ylabels = df_sync_hm.index
    g.set_yticklabels(ylabels, rotation=0, fontsize=8, ha='center', va='center')
    # annotate pairs that are significantly correlated
    df_stats['Module Pair'] = df_stats['Phospho Module'] + '+' + df_stats['RNA Module']
    sig_pairs = df_stats[df_stats['Corrected p-value (FDR)'] < 0.05]
    for i,row in sig_pairs.iterrows():
        phos_mod = row['Phospho Module']
        RNA_mod = row['RNA Module']
        x = list(df_sync_hm.columns).index(RNA_mod)+0.5
        y = list(df_sync_hm.index).index(phos_mod)+0.6
        pair = phos_mod + '+' + RNA_mod
        plt.annotate('*', xy=(x, y), fontsize=8, ha='center', va='center', color='w', weight='bold')
    plt.savefig('scores_sync.png', dpi=400, bbox_inches='tight')
    plt.close()
    print('done')
    
    ####################################################################################################3
    
    # plot pairs of modules together ("coplots")
    print('plotting coplots ... ', end='')
    x_ticks = [
        '0m', '15m+', '45m+', '90m+', '6h+', '12h+',
        '24h+', '36h+', '48h+', '72h+', '72h+/12h-',
        '72h+/24h-', '72h+/48h-', '72h+/72h-', '72h+/6d-'
    ]
    x_phos = [x_ticks.index(label.split(' Abundance')[0]) for label in df_phos_modules.columns]
    x_RNA = [x_ticks.index(label.split(' Abundance')[0]) for label in df_RNA_modules.columns]
    x_all = arange(len(x_ticks))
    module_coplots = []
    for i,row in df_stats.iterrows():
        if row['Corrected p-value (FDR)'] < 0.05: module_coplots.append((row['Phospho Module'], row['RNA Module']))
    for m_phos, m_RNA in module_coplots:
        fig = plt.figure(figsize=(6,2))
        y_phos = df_phos_modules.loc[m_phos]
        if n_phos_modules > 26: m = m_phos[0] + r'$_{' + m_phos[1:] + '}$'
        else: m = m_phos
        plt.plot(x_phos, y_phos, c='firebrick', label='Signaling\nmodule ' + m, zorder=2)
        error_phos = df_phos_error.loc[m_phos]
        plt.fill_between(x_phos, y_phos-error_phos, y_phos+error_phos, alpha=0.1, facecolor='firebrick', edgecolor='firebrick', zorder=2)
        y_RNA = df_RNA_modules.loc[m_RNA]
        if n_RNA_modules > 26: m = m_RNA[0] + r'$_{' + m_RNA[1:] + '}$'
        else: m = m_RNA
        plt.plot(x_RNA, y_RNA, c='mediumblue', label='Transcriptional\nmodule ' + m, zorder=2)
        error_RNA = df_RNA_error.loc[m_RNA]
        plt.fill_between(x_RNA, y_RNA-error_RNA, y_RNA+error_RNA, alpha=0.1, facecolor='mediumblue', edgecolor='mediumblue', zorder=2)
        plt.xticks(ticks=x_all, labels=x_ticks, rotation=30, ha='right', va='top')
        x_min, x_max = -0.2, len(x_all)-0.8
        plt.xlim(x_min,x_max)
        plt.xlabel('Time')
        plt.ylabel('Mean abundance\n(z-score)')
        plt.gca().spines['right'].set_visible(False)
        plt.gca().spines['top'].set_visible(False)
        plt.legend(frameon=False, bbox_to_anchor=(1.0,0.5), loc='center left')
        plt.axvline(x=9, color='k', linestyle='--', linewidth=1, zorder=1)
        xaxis_length = x_max-x_min
        plt.annotate('\nVEM on', xy=(0.2/xaxis_length, 1.05), xycoords='axes fraction', ha='left', va='center', fontsize=10, color='k', weight='bold')
        plt.annotate('\nVEM off', xy=((9+0.4)/xaxis_length, 1.05), xycoords='axes fraction', ha='left', va='center', fontsize=10, color='k', weight='bold')
        r = df_stats[(df_stats['Phospho Module'] == m_phos) & (df_stats['RNA Module'] == m_RNA)]['Pearson Correlation'].values[0]
        p = df_stats[(df_stats['Phospho Module'] == m_phos) & (df_stats['RNA Module'] == m_RNA)]['Corrected p-value (FDR)'].values[0]
        plt.annotate(r'$r^{2}='+str(round(r**2,2)) + '$', xy=(1.03,0.15), xycoords='axes fraction', ha='left', va='center', fontsize=8, color='k')
        plt.annotate(r'$P = {:0.3f}$'.format(p), xy=(1.03,0.05), xycoords='axes fraction', ha='left', va='center', fontsize=8, color='k')
        plt.savefig('dynamics_coplot_' + m_phos + m_RNA + '.png', dpi=400, bbox_inches='tight')
        plt.close()
    print('done')

if __name__ == '__main__':
    main()
