# ----------------------------------------------------------------------------
# Copyright (c) 2016-2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os.path
import pkg_resources

import biom
import q2templates
import numpy as np
import pandas as pd
import seaborn as sns


TEMPLATES = pkg_resources.resource_filename(
    'q2_feature_table._core_features', 'core_features_assets')


def core_features(output_dir, table: biom.Table, min_fraction: float = 0.5,
                  max_fraction: float = 1.0, steps: int = 11) -> None:
    if max_fraction < min_fraction:
        raise ValueError('min_fraction (%r) parameter must be less than '
                         'max_fraction (%r) parameter.' %
                         (min_fraction, max_fraction))

    index_fp = os.path.join(TEMPLATES, 'index.html')
    context = {
        'num_samples': table.shape[1],
        'num_features': table.shape[0]
    }

    if min_fraction == max_fraction:
        fractions = [min_fraction]
    else:
        fractions = np.linspace(min_fraction, max_fraction, steps)

    rounded_fractions = _round_fractions(fractions)

    data = []
    file_links = []
    for fraction, rounded_fraction in zip(fractions, rounded_fractions):
        core_features = _get_core_features(table, fraction)
        core_feature_count = len(core_features)
        data.append([fraction, core_feature_count])

        if core_feature_count > 0:
            core_feature_fn = 'core-features-%s.tsv' % rounded_fraction
            core_feature_fp = os.path.join(output_dir, core_feature_fn)

            file_links.append("<a href='./%s'>TSV</a>" % core_feature_fn)

            core_features.to_csv(core_feature_fp, sep='\t',
                                 index_label='Feature ID')
        else:
            file_links.append('No core features')

    df = pd.DataFrame(data, columns=['Fraction of samples', 'Feature count'])
    df['Fraction of features'] = df['Feature count'] / table.shape[0]
    df['Feature list'] = file_links

    # newer versions of seaborn don't like dataframes with fewer than two rows
    if len(fractions) > 1:
        ax = sns.regplot(data=df, x='Fraction of samples', y='Feature count',
                         fit_reg=False)

        # matplotlib will issue a UserWarning if attempting to set left and
        # right bounds to the same value.
        ax.set_xbound(min(fractions), max(fractions))
        ax.set_ybound(0, max(df['Feature count']) + 1)

        ax.get_figure().savefig(
            os.path.join(output_dir, 'core-feature-counts.svg'))
        context['show_plot'] = True

    context['table_html'] = q2templates.df_to_html(df, index=False,
                                                   escape=False)

    q2templates.render(index_fp, output_dir, context=context)


def _get_core_features(table, fraction):
    filter_f = _get_filter_to_core_f(table, fraction)
    feature_filtered_table = table.filter(
        filter_f, axis='observation', inplace=False)
    index = []
    data = []
    for values, id_, _ in feature_filtered_table.iter(axis='observation'):
        index.append(id_)
        data.append(_seven_number_summary(values))
    if len(data) > 0:
        return pd.DataFrame(data, index=index).sort_values(
            by='50%', ascending=False)
    else:
        return pd.DataFrame()


def _get_filter_to_core_f(table, fraction):
    # determine the number of samples that must have a non-zero value for a
    # feature to be considered part of the core
    min_count = fraction * table.shape[1]

    def f(values, *_):
        return np.count_nonzero(values) >= min_count
    return f


def _seven_number_summary(a):
    # TODO this should probably be publicly accessible throughout QIIME 2. it's
    # also currently implemented in q2-demux `summarize` and q2-diversity
    # `alpha-rarefaction`
    stats = pd.Series(a, dtype=float).describe(
        percentiles=[0.02, 0.09, 0.25, 0.5, 0.75, 0.91, 0.98])
    drop_cols = stats.index.isin(['std', 'mean', 'min', 'max', 'count'])
    stats = stats[~drop_cols]
    return stats


def _round_fractions(fractions):
    """Format float fractions as rounded fraction strings.

    Fractions are rounded as short as possible without losing the precision
    necessary to distinguish one fraction from another. The rounding starts at
    3 decimal places and continues up to 10 decimal places, if a shorter
    rounded representation isn't converged upon earlier. If rounding at 10
    decimal places doesn't provide enough precision, each fraction's full
    float precision will be used (i.e. as returned by `repr(float)`).

    For example, [0.12345, 0.12351] would be rounded to four decimal places:
    ['0.1234', '0.1235']. 3 decimal places would be tried first but there isn't
    enough precision to keep both fractions distinguishable.

    """
    for decimals in range(3, 11):
        rounded = ['{fraction:1.{decimals}f}'.format(fraction=fraction,
                                                     decimals=decimals)
                   for fraction in fractions]
        if len(rounded) == len(set(rounded)):
            return rounded
    return [repr(fraction) for fraction in fractions]
