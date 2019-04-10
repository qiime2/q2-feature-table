# ----------------------------------------------------------------------------
# Copyright (c) 2016-2019, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os.path
import pkg_resources

import q2templates
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import qiime2


TEMPLATES = pkg_resources.resource_filename('q2_feature_table._heatmap',
                                            'assets')


heatmap_choices = {
    'metric': {'braycurtis', 'canberra', 'chebyshev', 'cityblock',
               'correlation', 'cosine', 'dice', 'euclidean', 'hamming',
               'jaccard', 'kulsinski', 'mahalanobis', 'matching', 'minkowski',
               'rogerstanimoto', 'russellrao', 'seuclidean', 'sokalmichener',
               'sokalsneath', 'sqeuclidean', 'yule'},
    'method': {'single', 'complete', 'average', 'weighted', 'centroid',
               'median', 'ward'},
    'cluster': {'samples', 'features', 'both', 'none'},
    'color_scheme': {'Accent', 'Accent_r', 'Blues', 'Blues_r', 'BrBG',
                     'BrBG_r', 'BuGn', 'BuGn_r', 'BuPu', 'BuPu_r', 'CMRmap',
                     'CMRmap_r', 'Dark2', 'Dark2_r', 'GnBu', 'GnBu_r',
                     'Greens', 'Greens_r', 'Greys', 'Greys_r', 'OrRd',
                     'OrRd_r', 'Oranges', 'Oranges_r', 'PRGn', 'PRGn_r',
                     'Paired', 'Paired_r', 'Pastel1', 'Pastel1_r', 'Pastel2',
                     'Pastel2_r', 'PiYG', 'PiYG_r', 'PuBu', 'PuBuGn',
                     'PuBuGn_r', 'PuBu_r', 'PuOr', 'PuOr_r', 'PuRd', 'PuRd_r',
                     'Purples', 'Purples_r', 'RdBu', 'RdBu_r', 'RdGy',
                     'RdGy_r', 'RdPu', 'RdPu_r', 'RdYlBu', 'RdYlBu_r',
                     'RdYlGn', 'RdYlGn_r', 'Reds', 'Reds_r', 'Set1', 'Set1_r',
                     'Set2', 'Set2_r', 'Set3', 'Set3_r', 'Spectral',
                     'Spectral_r', 'Vega10', 'Vega10_r', 'Vega20', 'Vega20_r',
                     'Vega20b', 'Vega20b_r', 'Vega20c', 'Vega20c_r', 'Wistia',
                     'Wistia_r', 'YlGn', 'YlGnBu', 'YlGnBu_r', 'YlGn_r',
                     'YlOrBr', 'YlOrBr_r', 'YlOrRd', 'YlOrRd_r', 'afmhot',
                     'afmhot_r', 'autumn', 'autumn_r', 'binary', 'binary_r',
                     'bone', 'bone_r', 'brg', 'brg_r', 'bwr', 'bwr_r', 'cool',
                     'cool_r', 'coolwarm', 'coolwarm_r', 'copper', 'copper_r',
                     'cubehelix', 'cubehelix_r', 'flag', 'flag_r',
                     'gist_earth', 'gist_earth_r', 'gist_gray', 'gist_gray_r',
                     'gist_heat', 'gist_heat_r', 'gist_ncar', 'gist_ncar_r',
                     'gist_rainbow', 'gist_rainbow_r', 'gist_stern',
                     'gist_stern_r', 'gist_yarg', 'gist_yarg_r', 'gnuplot',
                     'gnuplot2', 'gnuplot2_r', 'gnuplot_r', 'gray', 'gray_r',
                     'hot', 'hot_r', 'hsv', 'hsv_r', 'icefire', 'icefire_r',
                     'inferno', 'inferno_r', 'jet', 'jet_r', 'magma',
                     'magma_r', 'mako', 'mako_r', 'nipy_spectral',
                     'nipy_spectral_r', 'ocean', 'ocean_r', 'pink', 'pink_r',
                     'plasma', 'plasma_r', 'prism', 'prism_r', 'rainbow',
                     'rainbow_r', 'rocket', 'rocket_r', 'seismic', 'seismic_r',
                     'spectral', 'spectral_r', 'spring', 'spring_r', 'summer',
                     'summer_r', 'tab10', 'tab10_r', 'tab20', 'tab20_r',
                     'tab20b', 'tab20b_r', 'tab20c', 'tab20c_r', 'terrain',
                     'terrain_r', 'viridis', 'viridis_r', 'vlag', 'vlag_r',
                     'winter', 'winter_r'}
}


_clustering_map = {'both': {'col_cluster': True, 'row_cluster': True},
                   'samples': {'col_cluster': False, 'row_cluster': True},
                   'features': {'col_cluster': True, 'row_cluster': False},
                   'none': {'col_cluster': False, 'row_cluster': False}}


def _munge_metadata(metadata, table, cluster):
    metadata = metadata.filter_ids(table.index)
    column_name = metadata.name

    metadata_df = metadata.to_dataframe()
    metadata_df = metadata_df.fillna('[No Value]')
    metadata_df['merged-id'] = metadata_df[column_name].str.cat(
        metadata_df.index, sep=' | ')
    # Inner join here because we have already validated that all sample IDs in
    # the table are present in the metadata, and the metadata has been filtered
    # to only include the table's IDs.
    table = table.join(metadata_df, how='inner')
    # It doesn't make sense to sort the samples if clustering is enabled on
    # the sample axis (e.g. `both` or `samples`).
    if cluster == 'features':
        table.sort_values(column_name, axis=0, ascending=True, inplace=True)
    table.set_index('merged-id', inplace=True)
    table.index.name = '%s | %s' % (column_name, metadata_df.index.name)
    table.drop([column_name], axis=1, inplace=True)
    return table


def heatmap(output_dir, table: pd.DataFrame,
            metadata: qiime2.CategoricalMetadataColumn = None,
            normalize: bool = True, title: str = None,
            metric: str = 'euclidean', method: str = 'average',
            cluster: str = 'both', color_scheme: str = 'rocket') -> None:
    if table.empty:
        raise ValueError('Cannot visualize an empty table.')

    if metadata is not None:
        table = _munge_metadata(metadata, table, cluster)

    cbar_label = 'frequency'
    if normalize:
        table = table.apply(lambda x: np.log10(x + 1))
        cbar_label = 'log10 frequency'

    # Hard-coded values for reasonable plots
    scaletron, labelsize, dpi = 50, 8, 100
    sns.set(rc={'xtick.labelsize': labelsize, 'ytick.labelsize': labelsize,
                'figure.dpi': dpi})
    width, height = table.shape[1] / scaletron, table.shape[0] / scaletron

    heatmap_plot = sns.clustermap(table, method=method, metric=metric,
                                  **_clustering_map[cluster],
                                  cmap=color_scheme,
                                  xticklabels=True, yticklabels=True,
                                  cbar_kws={'label': cbar_label})
    if title is not None:
        heatmap_plot.fig.suptitle(title)

    hm = heatmap_plot.ax_heatmap.get_position()
    cbar = heatmap_plot.cax.get_position()
    row = heatmap_plot.ax_row_dendrogram.get_position()
    col = heatmap_plot.ax_col_dendrogram.get_position()

    # Resize the plot to set cell aspect-ratio to square
    heatmap_plot.ax_heatmap.set_position([hm.x0, hm.y0, width, height])
    heatmap_plot.cax.set_position([cbar.x0, hm.y0 + height, cbar.width,
                                   cbar.height])
    heatmap_plot.ax_row_dendrogram.set_position([row.x0, row.y0, row.width,
                                                 height])
    heatmap_plot.ax_col_dendrogram.set_position([col.x0, hm.y0 + height, width,
                                                 col.height])

    # https://stackoverflow.com/a/34697479/3776794
    plt.setp(heatmap_plot.ax_heatmap.xaxis.get_majorticklabels(), rotation=90)
    plt.setp(heatmap_plot.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)

    for ext in ['png', 'svg']:
        img_fp = os.path.join(output_dir, 'feature-table-heatmap.%s' % ext)
        heatmap_plot.savefig(img_fp)

    index_fp = os.path.join(TEMPLATES, 'index.html')
    q2templates.render(index_fp, output_dir, context={'normalize': normalize})
