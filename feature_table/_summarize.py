# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------


import ipywidgets
from ipywidgets import interactive, fixed, IntSlider
from IPython.display import display
from scipy.optimize import minimize_scalar

def _summarize_even_sampling_depth(even_sampling_depth, counts):
    samples_retained = (counts >= even_sampling_depth)
    num_samples_retained = samples_retained.sum()
    num_sequences_retained = num_samples_retained * even_sampling_depth
    return samples_retained, num_samples_retained, num_sequences_retained

def _get_depth_for_max_sequence_count(counts):
    """Find the even sampling depth that retains the most sequences."""
    count_summary = counts.describe()
    def f(d):
        return -1 * _summarize_even_sampling_depth(d, counts)[2]

    res = minimize_scalar(f,
                          bounds=(count_summary['min'], count_summary['max']),
                          method='bounded')
    return int(np.floor(res.x))

def get_default_even_sampling_depth(biom):
    counts = biom.sum()
    return _get_depth_for_max_sequence_count(counts)

def explore_sampling_depth(biom):
    import seaborn as sns
    counts = biom.sum()
    count_summary = counts.describe()
    total_num_samples = len(counts)
    total_num_sequences = counts.sum()
    depth_for_max_sequence_count = _get_depth_for_max_sequence_count(counts)
    sampling_depth_slider = IntSlider(min=count_summary['min'],
                                      max=count_summary['max'],
                                      step=10 ** (math.log(count_summary['max'], 10) - 2),
                                      value=depth_for_max_sequence_count)
    default_samples_retained, default_num_samples_retained, default_num_sequences_retained = \
            _summarize_even_sampling_depth(depth_for_max_sequence_count, counts)

    default_percent_samples_retained = default_num_samples_retained * 100 / total_num_samples
    default_percent_sequences_retained = default_num_sequences_retained * 100 / total_num_sequences

    label_s = "Depth {0}: {1:.2f}% of sequences and {2:.2f}% of samples retained."

    def f(even_sampling_depth):
        samples_retained, num_samples_retained, num_sequences_retained = \
            _summarize_even_sampling_depth(even_sampling_depth, counts)
        percent_samples_retained = num_samples_retained * 100 / total_num_samples
        percent_sequences_retained = num_sequences_retained * 100 / total_num_sequences
        ax = sns.distplot(counts)
        ax.set_xlabel("Number of sequences per sample")
        ax.set_ylabel("Frequency")
        line_label = label_s.format(depth_for_max_sequence_count,
                                    default_percent_sequences_retained,
                                    default_percent_samples_retained)
        ax.plot([depth_for_max_sequence_count, depth_for_max_sequence_count], ax.get_ylim(),
                'k--', label=line_label)

        line_label = label_s.format(even_sampling_depth,
                                    percent_sequences_retained,
                                    percent_samples_retained)
        ax.plot([even_sampling_depth, even_sampling_depth], ax.get_ylim(),
                'k-', label=line_label)
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    def reset_depth(_):
        sampling_depth_slider.value = depth_for_max_sequence_count

    reset = ipywidgets.Button(icon='fa-refresh')
    reset.on_click(reset_depth)

    w = interactive(f, even_sampling_depth=sampling_depth_slider)
    display(ipywidgets.HBox(children=[w, reset]))
