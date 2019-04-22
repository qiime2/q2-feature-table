# ----------------------------------------------------------------------------
# Copyright (c) 2016-2019, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import json


def vega_spec(sample_metadata, sample_frequencies):

    values = []
    if sample_metadata:
        sample_metadata = sample_metadata.filter_ids(
            sample_frequencies.index)
        df = sample_metadata.to_dataframe()

        for i, row in df.iterrows():
            values.append({
                'id': i,
                'metadata': {j: row[j] for j in df.columns},
                'frequency': sample_frequencies[i]
            })

        metadata_categories = list(df.columns.values)
        max_frequency = int(max(sample_frequencies.values.tolist()))

    else:
        for freq in sample_frequencies:
            # add ids
            values.append({
                'frequency': i
            })
        metadata_categories = []
        max_frequency = 0

    spec = {
  "$schema": "https://vega.github.io/schema/vega/v4.json",
  "autosize": {"contains": "content", "type": "fit-x", "resize":True},
  "width": 800,
  "data": [
    {
      "transform": [
        {
          "as": "groupbyCat",
          "expr": "datum.metadata[category]",
          "type": "formula"
        }
      ],
      "name": "table",
      "values": values
    },
    {
      "transform": [
        {
          "type": "aggregate",
          "groupby": [
            "groupbyCat"
          ]
        }
      ],
      "name": "grouped",
      "source": "table"
    },
    {
      "transform": [
        {
          "expr": "datum.frequency >= SamplingDepth",
          "type": "filter"
        },
        {
          "type": "aggregate",
          "groupby": [
            "groupbyCat"
          ]
        }
      ],
      "name": "groupedAndRetained",
      "source": "table"
    }
  ],
  "signals": [
    {"name": "chartHeight", "value": 400},
    {"name": "chartOffset", "value": 20},
    {"name": "height", "update": "chartHeight + chartOffset"},
    {
      "name": "tooltip",
      "value": {},
      "on": [
        {
          "events": "rect:mouseover",
          "update": "datum"
        },
        {
          "events": "rect:mouseout",
          "update": "{}"
        }
      ]
    },
    {
      "bind": {
        "input": "select",
        "options": metadata_categories,
        "element": "#metadata-category"
      },
      "name": "category",
      "value": metadata_categories[0]
    },
    {
      "bind": {
        "input": "range",
        "element": "#sampling-depth-slider",
        "max": max_frequency,
        "min": 0,
        "step": 1
      },
      "name": "SamplingDepth",
      "value": 0
    },
    {
      "name": "rotateLabelsCheckbox",
      "value": False,
      "bind": {
        "input": "checkbox",
        "element": "#rotate-labels"
      }
    },
    {
      "name": "rotateLabels",
      "update": "if(rotateLabelsCheckbox, 90, 0)"
    }
  ],
  "marks": [
    {
      "description": "Interactive rarefaction summary plot",
      "name": "rarefactionplot",
      "type": "group",
      "encode": {
        "enter": {
          "x": {
            "value": 0
          },
          "y": {
            "signal": "chartOffset"
          },
          "width": {
            "value": 200
          },
          "height": {
            "signal": "chartHeight"
          }
        }
      },
      "marks": [
        {
          "type": "rect",
          "from": {
            "data": "grouped"
          },
          "encode": {
            "enter": {
              "x": {
                "scale": "xscale",
                "field": "groupbyCat"
              },
              "width": {
                "scale": "xscale",
                "band": 1
              },
              "y": {
                "scale": "yscale",
                "field": "count"
              },
              "y2": {
                "scale": "yscale",
                "value": 0
              }
            },
            "update": {
              "fill": {
                "value": "#D3D3D3"
              }
            }
          }
        },
        {
          "type": "rect",
          "from": {
            "data": "groupedAndRetained"
          },
          "encode": {
            "enter": {
              "x": {
                "scale": "xscale",
                "field": "groupbyCat"
              },
              "width": {
                "scale": "xscale",
                "band": 1
              },
              "y": {
                "scale": "yscale",
                "field": "count"
              },
              "y2": {
                "scale": "yscale",
                "value": 0
              }
            },
            "update": {
              "y": {
                "scale": "yscale",
                "field": "count"
              },
              "y2": {
                "scale": "yscale",
                "value": 0
              },
              "fill": {
                "value": "steelBlue"
              }
            },
            "hover": {
              "fill": {
                "value": "red"
              }
            }
          }
        },
        {
          "type": "text",
          "encode": {
            "enter": {
              "align": {
                "value": "center"
              },
              "baseline": {
                "value": "top"
              },
              "fill": {
                "value": "#333"
              }
            },
            "update": {
              "x": {
                "scale": "xscale",
                "signal": "tooltip.groupbyCat",
                "band": 0.5
              },
              "y": {
                "scale": "yscale",
                "signal": "tooltip.count",
                "offset": -12
              },
              "text": {
                "signal": "tooltip.count"
              },
              "fillOpacity": [
                {
                  "test": "datum === tooltip",
                  "value": 0
                },
                {
                  "value": 1
                }
              ]
            }
          }
        }
      ],
      "scales": [
        {
          "name": "xscale",
          "type": "band",
          "domain": {
            "data": "grouped",
            "field": "groupbyCat",
            "sort": True
          },
          "range": "width",
          "padding": 0.1,
          "round": True
        },
        {
          "name": "yscale",
          "domain": {
            "fields": [
              {
                "data": "grouped",
                "field": "count"
              }
            ]
          },
          "nice": True,
          "range": [{"signal": "chartHeight"}, 0]
        }
      ],
      "axes": [
        {
          "orient": "bottom",
          "scale": "xscale",
          "title": {
            "signal": "category"
          },
          "labelAngle": {"signal": "rotateLabels"},
          "labelAlign": "left"
        },
        {
          "orient": "left",
          "scale": "yscale",
          "title": "Number of Samples"
        }
      ]
    }
  ]
}
    vega_spec = json.dumps(spec)
    return vega_spec
