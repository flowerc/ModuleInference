from math import floor
import matplotlib as mpl
import matplotlib.pyplot as plt
from numpy import array
from pandas import concat, DataFrame, read_excel
from scipy.cluster.hierarchy import distance, fcluster, linkage
from scipy.stats import fisher_exact, rankdata
import seaborn as sns
from string import ascii_uppercase
# plot formatting
sns.set(font_scale=1.5)
mpl.style.use('default')
alphabet = list(ascii_uppercase)

def main():
    # read in phospho data
    print('importing and formatting phosphoproteomics data ... ', end='')
    df_phos = read_excel('Supplementary Table S1.xlsx', sheet_name='Time-series phospho - overlap')
    df_phos = df_phos.set_index('Phosphosites').sort_index()
    df_phos = df_phos.drop('Protein Descriptions', axis=1)
    n_sites = len(df_phos.index)
    # average over all three replicates
    done = []
    for c in df_phos.columns:
        condition = c.split(' BR')[0]
        if condition in done: continue
        cols = [condition + ' BR' + str(BR) + '\n(TMT intensity)' for BR in [1,2,3]]
        df_phos[condition] = df_phos[cols].mean(axis=1)
        df_phos = df_phos.drop(cols, axis=1)
        if condition not in done: done.append(condition)
    # drop DMSO samples
    df_phos = df_phos.drop(['24h-', '72h-', '6d-', '9d-'], axis=1)
    # z-score normalize (doesn't change clustering/correlation but will be necessary for plotting module dynamics later)
    df_phos = df_phos.sub(df_phos.mean(axis=1), axis=0).div(df_phos.std(axis=1), axis=0)
    print('done')
    
    # calculate pairwise correlation matrix
    print('calculating correlation matrix ... ', end='')
    df_phos_corr = df_phos.T.corr(method='pearson')
    print('done')
    
    # calculate distances and linkages for clustering
    print('calculating phospho linkages ... ', end='')
    d = distance.pdist(df_phos_corr, metric='euclidean')
    L = linkage(d, method='ward')   # change method to average to reproduce seaborn default clustering
    print('done')
    
    # plot clustered heatmap
    print('making clustered heatmap ... ', end='')
    g = sns.clustermap(df_phos_corr, row_linkage=L, col_linkage=L, cmap='seismic', vmin=-1, vmax=1, center=0,
                       xticklabels=False, yticklabels=False, cbar_kws={'ticks': [-1,0,1], 'orientation': 'horizontal'})
    plt.title("Pearson's " + r'$r$', size=16)
    print('done')
    
    # extract modules
    print('extracting modules ... ', end='')
    dist_threshold_factor = 0.037        # increasing the value of this hyerparameter results in fewer, larger modules. decreasing has the opposite effect.
    modules = fcluster(L, t=dist_threshold_factor*n_sites, criterion='distance')
    df_modules = DataFrame()
    df_modules['Site'] = df_phos_corr.index
    n_modules = len(set(modules))
    if n_modules > 26: module_labels = [alphabet[m-1-floor((m-1)/26)*26] + str(floor((m-1)/26)+1) for m in modules]
    else: module_labels = [alphabet[m-1] for m in modules]
    df_modules['Module'] = module_labels
    reordered_row_incides = g.dendrogram_row.reordered_ind
    df_phos_corr = df_phos_corr.reset_index()
    df_phos_corr = df_phos_corr.reset_index().reindex(reordered_row_incides)
    df_phos_corr = df_phos_corr.reset_index().sort_values('Phosphosites')
    df_modules.index = df_phos_corr.index
    df_modules = df_modules.sort_index()
    print('done')
    
    # format heatmap
    print('exporting unannotated heatmap ... ', end='')
    g.ax_col_dendrogram.set_visible(False)
    g.ax_cbar.set_position([0.05, 0.78, g.ax_row_dendrogram.get_position().width*0.65, 0.01])
    row_dend = g.ax_row_dendrogram.get_position()
    g.ax_row_dendrogram.set_position([row_dend.x0+row_dend.width*(1-0.6), row_dend.y0, row_dend.width*0.6, row_dend.height])
    g = g.ax_heatmap
    g.set_xlabel('Phosphosites', size=24)
    g.set_ylabel('Phosphosites', size=24)
    plt.savefig('corr_phos_unannotated.png', bbox_inches='tight', dpi=400)
    print('done')
    
    # annotate modules
    print('annotating modules ... ', end='')
    margin = n_sites/100
    df_module_dynamics = DataFrame(columns=df_phos.columns)
    df_module_error = DataFrame(columns=df_phos.columns)
    for i,m in enumerate(df_modules['Module'].unique()):
        module_members = df_modules[df_modules['Module'] == m]
        module_start = min(module_members.index)
        module_end = max(module_members.index)
        
        # plot boxes and labels around each module
        if module_start == 0: module_start = 0.5
        if module_end == n_sites-1: module_end -= 0.5
        g.hlines(module_start, module_start, module_end+1, linewidth=3, color='k')
        g.hlines(module_end+1, module_start, module_end+1, linewidth=3, color='k')
        g.vlines(module_start, module_start-margin/5, module_end+1+margin/5, linewidth=3, color='k')
        g.vlines(module_end+1, module_start-margin/5, module_end+1+margin/5, linewidth=3, color='k') 
        if n_modules > 26: label = m[0] + r'$_{' + m[1:] + '}$'
        else: label = m
        if i % 2 == 0:
            # position label in bottom left
            g.text(module_start-margin, module_end+margin, label, ha='right', va='top', fontsize=18,
                   color='w', weight='bold', bbox=dict(facecolor='k', edgecolor='k', boxstyle='round,pad=0.1'))
        else:
            # position label in top right
            g.text(module_end+margin, module_start-margin, label, ha='left', va='bottom', fontsize=18,
                   color='w', weight='bold', bbox=dict(facecolor='k', edgecolor='k', boxstyle='round,pad=0.1'))
        
        # add consensus module dynamics and error bars to dataframes to export after this loop
        m_timeseries = df_phos.loc[module_members['Site']]
        df_module_dynamics.loc[m] = m_timeseries.mean()
        df_module_error.loc[m] = m_timeseries.std()
        
    g.set_title('Inferred signaling modules\n', size=28)
    plt.savefig('corr_phos_annotated.png', bbox_inches='tight', dpi=400)
    plt.close()
    print('done')
    
    # perform representation analysis for a defined set of libraries
    libraries = [
        'library_BioCarta_v2023_1',
        'library_GObp_v2023_1',
        'library_GOcc_v2023_1',
        'library_GOmf_v2023_1',
        'library_Hallmarks_v2023_1',
        'library_KEGG_v2023_1',
        'library_PID_v2023_1',
        'library_Reactome_v2023_1',
        'library_TFT_v2023_1',
        'library_WikiPathways_v2023_1'
    ]
    lib_title_map = {
        'library_BioCarta_v2023_1': 'BioCarta Pathways',
        'library_GObp_v2023_1': 'GO Biological Processes',
        'library_GOcc_v2023_1': 'GO Cellular Compartments',
        'library_GOmf_v2023_1': 'GO Molecular Functions',
        'library_Hallmarks_v2023_1': 'MSigDB Hallmark Gene Sets',
        'library_KEGG_v2023_1': 'KEGG Pathways',
        'library_PID_v2023_1': 'PID Pathways',
        'library_Reactome_v2023_1': 'Reactome Pathways',
        'library_TFT_v2023_1': 'Transcription Factor Targets',
        'library_WikiPathways_v2023_1': 'WikiPathways'
    }
    
    # before performing gene set overrepresentation analysis in the following loop,
    # get gene background from the data (all sites in the dataset that map uniquely to one gene)
    background = set([s for s in df_modules['Site'] if ';' not in s])
    
    # perform phosphosite-set representation analysis
    print('phosphosite-set representation analysis:')
    df_rep_list = []
    for lib in libraries:
        lib_name_reformat = lib_title_map[lib]
        lib_filename = 'libraries\\' + lib + '.gmt'
        for m in df_modules['Module'].unique():
            print('library: ' + lib_name_reformat + '\tphos module: ' + m + ' ... ', end='')
            # get foreground (all sites in the module that map uniquely to one gene)
            foreground = list(df_modules[df_modules['Module'] == m]['Site'])
            foreground = set([s for s in foreground if ';' not in s])
            
            module_list, lib_list, gene_set_name_list, OR_list, p_list, overlap_list, overlap_IDs_list = [], [], [], [], [], [], []
            with open(lib_filename) as file:
                for line in file:
                    line_list = line.rstrip().split('\t')
                    query_set_name = line_list[0]
                    query_set = line_list[2:]
                    query_set = [s for s in background if s.split('-p')[0] in query_set]
                    query_set = set(query_set).intersection(background)
                    if len(query_set) == 0: continue
                    # perform one-sided Fisher's exact test
                    M = len(background)
                    n = len(query_set)
                    N = len(foreground)
                    overlap = sorted(list(foreground.intersection(query_set)))
                    k = len(overlap)
                    table = array([[k, n-k], [N-k, M-(n+N)+k]])
                    odds_ratio, p = fisher_exact(table, alternative='greater')
                    gene_set_name_list.append(query_set_name)
                    OR_list.append(odds_ratio)
                    p_list.append(p)
                    module_list.append(m)
                    lib_list.append(lib_name_reformat)
                    overlap_list.append("'" + str(k) + '/' + str(n) + "'")
                    overlap_IDs_list.append('; '.join(overlap))
            
            # adjusted p-values by BH procedure.
            pval_rank = rankdata(p_list)
            n_pvals = len(p_list)
            FDR_list = [0]*n_pvals
            for i in range(n_pvals):
                FDR = p_list[i]*n_pvals/pval_rank[i]
                if FDR > 1: FDR = 1
                FDR_list[i] = FDR
            
            df_rep = DataFrame()
            df_rep['Module'] = module_list
            df_rep['Library'] = lib_list
            df_rep['Gene Set'] = gene_set_name_list
            df_rep['Overlap'] = overlap_list
            df_rep['Sites'] = overlap_genes_list
            df_rep['Odds Ratio'] = OR_list
            df_rep['P-value'] = p_list
            df_rep['Corrected p-value'] = FDR_list
            df_rep = df_rep[df_rep['Corrected p-value'] < 0.05]
            df_rep_list.append(df_rep)
            print('done')
    
    # concatenate dataframes to create a master df for representation analysis,
    # re-sort by FDR, export
    df_rep = concat(df_rep_list).reset_index(drop=True)
    df_rep = df_rep.sort_values('Corrected p-value')
    print('\n')
    
    # export module members, module dynamics, and module representation stats
    df_modules.to_csv('phospho_module_members.csv', index=False)
    df_module_dynamics.to_csv('phospho_module_dynamics.csv')
    df_module_error.to_csv('phospho_module_error.csv')
    df_rep.to_csv('geneset_representation_stats_phos.csv', index=False)
    
    ############################################################################################################
    
    # read in transcriptome data
    print('importing and formatting RNA-seq data ... ', end='')
    df_RNA = read_excel('Supplementary Table S2.xlsx', sheet_name='Time-series RNA-seq - FPKM')
    df_RNA = df_RNA.set_index('Gene ID').sort_index()
    df_RNA = df_RNA.drop('Ensembl ID', axis=1)
    col_names = [c.split('\n(FPKM)')[0] for c in df_RNA.columns]
    df_RNA = df_RNA.rename(columns={k:v for k,v in zip(df_RNA.columns, col_names)})
    # filter out transcripts with 0s
    df_RNA['Min FPKM'] = df_RNA.min(axis=1)
    df_RNA = df_RNA[df_RNA['Min FPKM'] > 0]
    df_RNA = df_RNA.drop('Min FPKM', axis=1)
    n_genes = len(df_RNA.index)
    # zscore normalize (doesn't change clustering/correlation but will be necessary for plotting module dynamics later)
    df_RNA = df_RNA.sub(df_RNA.mean(axis=1), axis=0).div(df_RNA.std(axis=1), axis=0)
    print('done')
    
    # calculate pairwise correlation matrix
    print('calculating correlation matrix ... ', end='')
    df_RNA_corr = df_RNA.T.corr(method='pearson')
    print('done')
    
    # calculate distances and linkages for clustering
    # (this can take some time)
    print('calculating RNA-seq linkages ... ', end='')
    d = distance.pdist(df_RNA_corr, metric='euclidean')
    L = linkage(d, method='ward')
    print('done')
    
    # plot clustered heatmap
    print('making clustered heatmap ... ', end='')
    g = sns.clustermap(df_RNA_corr, row_linkage=L, col_linkage=L, cmap='seismic', vmin=-1, vmax=1, center=0,
                       xticklabels=False, yticklabels=False, cbar_kws={'ticks': [-1,0,1], 'orientation': 'horizontal'})
    plt.title("Pearson's " + r'$r$', size=16)
    print('done')
    
    # extract modules
    print('extracting modules ... ', end='')
    modules = fcluster(L, t=dist_threshold_factor*n_genes, criterion='distance')
    df_modules = DataFrame()
    df_modules['Gene ID'] = df_RNA_corr.index
    n_modules = len(set(modules))
    if n_modules > 26: module_labels = [alphabet[m-1-floor((m-1)/26)*26] + str(floor((m-1)/26)+1) for m in modules]
    else: module_labels = [alphabet[m-1] for m in modules]
    df_modules['Module'] = module_labels
    reordered_row_incides = g.dendrogram_row.reordered_ind
    df_RNA_corr = df_RNA_corr.reset_index()
    df_RNA_corr = df_RNA_corr.reset_index().reindex(reordered_row_incides)
    df_RNA_corr = df_RNA_corr.reset_index().sort_values('Gene ID')
    df_modules.index = df_RNA_corr.index
    df_modules = df_modules.sort_index()
    print('done')
    
    # format heatmap
    print('exporting unannotated heatmap ... ', end='')
    g.ax_col_dendrogram.set_visible(False)
    g.ax_cbar.set_position([0.05, 0.78, g.ax_row_dendrogram.get_position().width*0.65, 0.01])
    row_dend = g.ax_row_dendrogram.get_position()
    g.ax_row_dendrogram.set_position([row_dend.x0+row_dend.width*(1-0.6), row_dend.y0, row_dend.width*0.6, row_dend.height])
    g = g.ax_heatmap
    g.set_xlabel('Genes', size=24)
    g.set_ylabel('Genes', size=24)
    plt.savefig('corr_RNA_unannotated.png', bbox_inches='tight', dpi=400)
    print('done')
    
    # annotate modules
    print('annotating heatmap ... ', end='')
    margin = n_genes/100
    df_module_dynamics = DataFrame(columns=df_RNA.columns)
    df_module_error = DataFrame(columns=df_RNA.columns)
    for i,m in enumerate(df_modules['Module'].unique()):
        module_members = df_modules[df_modules['Module'] == m]
        module_start = min(module_members.index)
        module_end = max(module_members.index)
        # plot boxes and labels around each module
        if module_start == 0: module_start = 0.5
        if module_end == n_genes-1: module_end -= 0.5
        g.hlines(module_start, module_start, module_end+1, linewidth=3, color='k')
        g.hlines(module_end+1, module_start, module_end+1, linewidth=3, color='k')
        g.vlines(module_start, module_start-margin/5, module_end+1+margin/5, linewidth=3, color='k')
        g.vlines(module_end+1, module_start-margin/5, module_end+1+margin/5, linewidth=3, color='k')
        if n_modules > 26: label = m[0] + r'$_{' + m[1:] + '}$'
        else: label = m
        if i % 2 == 0:
            # position label in bottom left
            g.text(module_start-margin, module_end+margin, label, ha='right', va='top', fontsize=18,
                   color='w', weight='bold', bbox=dict(facecolor='k', edgecolor='k', boxstyle='round,pad=0.1'))
        else:
            # position label in top right
            g.text(module_end+margin, module_start-margin, label, ha='left', va='bottom', fontsize=18,
                   color='w', weight='bold', bbox=dict(facecolor='k', edgecolor='k', boxstyle='round,pad=0.1'))
        # add consensus module dynamics and error bars to dataframes to export after this loop
        m_timeseries = df_RNA.loc[module_members['Gene ID']]
        df_module_dynamics.loc[m] = m_timeseries.mean()
        df_module_error.loc[m] = m_timeseries.std()
    
    g.set_title('Inferred transcriptional modules\n', size=28)
    plt.savefig('corr_RNA_annotated.png', bbox_inches='tight', dpi=400)
    plt.close()
    print('done')
    
    # before performing gene set representation analysis in the following loop,
    # get gene background from the data (all genes in the dataset)
    background = set(df_modules['Gene ID'])
    
    # perform gene set representation analysis
    print('gene set representation analysis:')
    df_rep_list = []
    for lib in libraries:
        lib_name_reformat = lib_title_map[lib]
        lib_filename = 'libraries\\' + lib + '.gmt'
        for m in df_modules['Module'].unique():
            print('library: ' + lib_name_reformat + '\tRNA module: ' + m + ' ... ', end='')
            # get query gene list from the module
            foreground = set(df_modules[df_modules['Module'] == m]['Gene ID'])
            
            module_list, lib_list, gene_set_name_list, OR_list, p_list, overlap_list, overlap_IDs_list = [], [], [], [], [], [], []
            with open(lib_filename) as file:
                for line in file:
                    line_list = line.rstrip().split('\t')
                    query_set_name = line_list[0]
                    query_set = line_list[2:]
                    query_set = set(query_set).intersection(background)
                    if len(query_set) == 0: continue
                    # perform one-sided Fisher's exact test
                    M = len(background)
                    n = len(query_set)
                    N = len(foreground)
                    overlap = sorted(list(foreground.intersection(query_set)))
                    k = len(overlap)
                    table = array([[k, n-k], [N-k, M-(n+N)+k]])
                    odds_ratio, p = fisher_exact(table, alternative='greater')
                    gene_set_name_list.append(query_set_name)
                    OR_list.append(odds_ratio)
                    p_list.append(p)
                    module_list.append(m)
                    lib_list.append(lib_name_reformat)
                    overlap_list.append("'" + str(k) + '/' + str(n) + "'")
                    overlap_IDs_list.append('; '.join(overlap))
            
            # adjusted p-values by BH procedure.
            pval_rank = rankdata(p_list)
            n_pvals = len(p_list)
            FDR_list = [0]*n_pvals
            for i in range(n_pvals):
                FDR = p_list[i]*n_pvals/pval_rank[i]
                if FDR > 1: FDR = 1
                FDR_list[i] = FDR
            
            df_rep = DataFrame()
            df_rep['Module'] = module_list
            df_rep['Library'] = lib_list
            df_rep['Gene Set'] = gene_set_name_list
            df_rep['Overlap'] = overlap_list
            df_rep['Genes'] = overlap_genes_list
            df_rep['Odds Ratio'] = OR_list
            df_rep['P-value'] = p_list
            df_rep['Corrected p-value'] = FDR_list
            df_rep = df_rep[df_rep['Corrected p-value'] < 0.05]
            df_rep_list.append(df_rep)
            print('done')
    
    # concatenate dataframes to create a master df for representation analysis,
    # re-sort by FDR, export
    df_rep = concat(df_rep_list).reset_index(drop=True)
    df_rep = df_rep.sort_values('Corrected p-value')
    print('\n')
    
    # export module members and module dynamics
    df_modules.to_csv('RNA_module_members.csv', index=False)
    df_module_dynamics.to_csv('RNA_module_dynamics.csv')
    df_module_error.to_csv('RNA_module_error.csv')
    df_rep.to_csv('geneset_representation_stats_RNA.csv', index=False)

# main code
if __name__ == '__main__':
    main()
