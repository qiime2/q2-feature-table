import json
import os
import pkg_resources
from distutils.dir_util import copy_tree

import qiime2
import q2templates

# data is being passed dynamically

def vega_spec(table, sample_metadata, sample_summary, sample_frequencies):
    values = []

    if sample_metadata:
        sample_metadata = sample_metadata.filter_ids(
            sample_frequencies.index)
        pandadataframe = sample_metadata.to_dataframe()

        for i, row in pandadataframe.iterrows():
            values.append({
            'id': i,
            'metadata': {j: row[j] for j in pandadataframe.columns},
            'frequency': sample_frequencies[i]
            })
        json.dumps(values)

    else:
        for i in sample_frequencies:
            # add ids
            values.append({
            'frequency': i
            })
        json.dumps(values)
    # pass data as json format to the vega_spec function
    context = dict()
    spec = {
 "$schema": "https://vega.github.io/schema/vega/v4.json",
  "width": 600,
  "height": 400,
  "padding": 5,
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
      "name": "grouped2",
      "source": "table"
    }
  ],
  "signals": [
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
        "options": [
          # input all metadata categories here
        ]
      },
      "name": "category",
      "value": "Day"
    },
    {
      "bind": {
        "input": "range",
        "max": 10095,
        "min": 0,
        "step": 1
      },
      "name": "SamplingDepth",
      "value": 0,
      "update": "SamplingDepthValue"
    },
    # add a JS event listener to update these two signals together
    {
      "name": "SamplingDepthValue",
      "value": "",
       "bind": {
       "input": "number"
      }
    },

    {
      "name": "RotateLabels",
      "value": False,
      "bind" :{
        "input": "checkbox"
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
          {"data": "grouped", "field":"count"}
        ]
      },
      "nice": True,
      "range": "height"
    }
  ],
  "axes": [
    {
      "orient": "bottom",
      "scale": "xscale",
      "title": {
        "signal": "category"
      },
      "labelAngle": 45,
      "labelAlign":"left"
    },
    {
      "orient": "left",
      "scale": "yscale",
      "title": "Number of Samples"
    }
  ],
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
        "data": "grouped2"
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
          "align": {"value": "center"},
          "baseline": {"value": "top"},
          "fill": {"value": "#333"}
        },
        "update": {
          "x": {"scale": "xscale", "signal": "tooltip.groupbyCat", "band": 0.5},
          "y": {"scale": "yscale", "signal": "tooltip.count", "offset": -12},
          "text": {"signal": "tooltip.count"},
          "fillOpacity": [
            {"test": "datum === tooltip", "value": 0},
            {"value": 1}
          ]
        }
      }
    }
  ]
}
    vega_spec = json.dumps(spec)
    return vega_spec
    # TEMPLATES = pkg_resources.resource_filename('q2_feature_table', 'vega')
    # context['vega_spec'] = json.dumps(spec)
    #
    # copy_tree(os.path.join(TEMPLATES, 'vega'), output_dir)
    # index = os.path.join(TEMPLATES, 'index.html')
    # q2templates.render(index, output_dir, context=context)
