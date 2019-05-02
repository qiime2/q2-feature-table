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
        # create data in json format for Vega
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
            "autosize": {"contains": "content",
                         "type": "fit-x",
                         "resize": True},
            "width": 800,
            "data": [
                {
                  "transform": [
                    {
                      "as": "selectedCategory",
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
                        "selectedCategory"
                      ]
                    }
                  ],
                  "name": "groupedBySelectedCategory",
                  "source": "table"
                },
                {
                  "transform": [
                    {
                      "expr": "datum.frequency >= samplingDepth",
                      "type": "filter"
                    },
                    {
                      "type": "aggregate",
                      "groupby": [
                        "selectedCategory"
                      ]
                    }
                  ],
                  "name": "groupedBySelectedCategoryAndRetained",
                  "source": "table"
                }
              ],
            "signals": [
                    {"name": "chartHeight", "update": "width / 2.5"},
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
                      "name": "samplingDepth",
                      "value": 0
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
                        "data": "groupedBySelectedCategory"
                      },
                      "encode": {
                        "enter": {
                          "x": {
                            "scale": "xscale",
                            "field": "selectedCategory"
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
                          "tooltip": {
                            "signal": (
                                "{\"title\": datum.selectedCategory, "
                                "\"Samples Dropped\": datum.count}")
                          },
                          "fill": {
                            "value": "#D3D3D3"
                          }
                        }
                      }
                    },
                    {
                      "type": "rect",
                      "from": {
                        "data": "groupedBySelectedCategoryAndRetained"
                      },
                      "encode": {
                        "enter": {
                          "x": {
                            "scale": "xscale",
                            "field": "selectedCategory"
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
                          "tooltip": {
                            "signal": (
                              "{\"title\": datum.selectedCategory, "
                              "\"Samples Retained\": datum.count}")
                              },
                          "y": {
                            "scale": "yscale",
                            "field": "count"
                          },
                          "y2": {
                            "scale": "yscale",
                            "value": 0
                          },
                          "fill": {
                            "value": "#4682B4"
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
                            "value": "#000000"
                          }
                        },
                        "update": {
                          "x": {
                            "scale": "xscale",
                            "signal": "tooltip.selectedCategory",
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
                        "data": "groupedBySelectedCategory",
                        "field": "selectedCategory",
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
                            "data": "groupedBySelectedCategory",
                            "field": "count"
                          }
                        ]
                      },
                      "nice": True,
                      "range": [
                        {
                          "signal": "chartHeight"
                        },
                        0
                      ]
                    }
                  ],
                  "axes": [
                    {
                      "orient": "bottom",
                      "scale": "xscale",
                      "title": {
                        "signal": "category"
                      },
                      "labelAngle": 90,
                      "labelAlign": "left",
                      "labelBaseline": "middle"},
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
