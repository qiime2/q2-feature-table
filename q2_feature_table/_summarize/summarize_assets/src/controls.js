// ----------------------------------------------------------------------------
// Copyright (c) 2016-2017, QIIME 2 development team.
//
// Distributed under the terms of the Modified BSD License.
//
// The full license is in the file LICENSE, distributed with this software.
// ----------------------------------------------------------------------------

import * as d3 from 'd3';


const getRetainedSamples = (trows, slider) => {
  if (slider.node().value == 0) {
    // Explicitly handle 0, as this means no rarefaction (so all samples are retained).
    return trows.data().length;
  } else {
    return trows.data().reduce((i, j) => i + (+j[1] >= slider.node().value ? 1 : 0), 0);
  }
};


const calcSampleRetainment = (trows, slider) => {
  const retained = getRetainedSamples(trows, slider);
  return [retained.toLocaleString('en-US'), ((retained / trows.data().length) * 100).toFixed(2)];
};


const calcFeatureRetainment = (trows, slider) => {
  if (slider.node().value == 0) {
    // Explicitly handle 0, as this means no rarefaction (so all sequences are retained).
    const retained = d3.sum(trows.data(), d => d[1]);
    return [retained.toLocaleString('en-US'), (100).toFixed(2)];
  } else {
    const retained = getRetainedSamples(trows, slider) * slider.node().value;
    return [retained.toLocaleString('en-US'), ((retained / d3.sum(trows.data(), d => d[1])) * 100).toFixed(2)];
  }
};


const toggleStyles = (trows, slider) => (
  trows
    .attr('class', d => (d[1] < slider.node().value ? 'alert-danger' : ''))
);


const initializeControls = () => {
  const slider = d3.select('#slider');
  const sliderValue = d3.select('#slider-value');
  const trows = d3.select('tbody').selectAll('tr');
  const formGroup = d3.select('.form-group');
  const hiddenDepth = d3.select('#hidden-depths');

  const updateStats = (trows, slider) => {
    formGroup.selectAll('span')
      .data([['Retained ', `${calcFeatureRetainment(trows, slider).join(' (')}%)`, ' sequences in ', `${calcSampleRetainment(trows, slider).join(' (')}%)`, ' samples at the specifed sampling depth.']])
      .text(d => d.join(''))
      .insert('br');
  };

  sliderValue.on('input', () => {
    if (+sliderValue.node().value > +slider.node().max) {
      sliderValue.node().value = slider.node().max;
      slider.node().value = slider.node().max;
    }
    else if (sliderValue.node().value < 0) {
      sliderValue.node().value = 0;
      slider.node().value = slider.node().min;
    }
    else {
      slider.node().value = sliderValue.node().value;
    }
    updateStats(trows, slider);
    toggleStyles(trows, slider);
    hiddenDepth.attr('value', slider.node().value).dispatch('change');
  });

  sliderValue.on('change', () => {
    if (!sliderValue.node().value) {
      sliderValue.node().value = 0;
      slider.node().value = slider.node().min;
    }
    updateStats(trows, slider);
    toggleStyles(trows, slider);
    hiddenDepth.attr('value', slider.node().value).dispatch('change');
  });

  slider.on('input', () => {
    sliderValue.node().value = slider.node().value;
    toggleStyles(trows, slider);
    updateStats(trows, slider);
    hiddenDepth.attr('value', slider.node().value).dispatch('change');
  });

  formGroup
      .append('div')
    .attr('display', 'inline-block')
    .style('padding-top', '10px')
      .selectAll('span')
      .data([['Retained ', `${calcFeatureRetainment(trows, slider).join(' (')}%)`, ' sequences in ', `${calcSampleRetainment(trows, slider).join(' (')}%)`, ' samples at the specifed sampling depth.']])
      .enter()
    .append('span')
    .text(d => d.join(''))
    .insert('br');

  sliderValue.node().value = slider.node().value;
};

export default initializeControls;
