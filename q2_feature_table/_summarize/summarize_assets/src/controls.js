// ----------------------------------------------------------------------------
// Copyright (c) 2016-2017, QIIME 2 development team.
//
// Distributed under the terms of the Modified BSD License.
//
// The full license is in the file LICENSE, distributed with this software.
// ----------------------------------------------------------------------------

import * as d3 from 'd3';


const calcSampleLoss = (trows, slider) => {
  const lost = trows.data().reduce((i, j) => i + (j[1] < +slider.node().value ? 1 : 0), 0);
  return [lost, ((lost / trows.data().length) * 100).toFixed(2)];
};


const toggleStyles = (trows, slider) => (
  trows
    .attr('class', d => (+d[1] < +slider.node().value ? 'alert-danger' : ''))
);


const initializeControls = () => {
  const slider = d3.select('#slider');
  const sliderValue = d3.select('#slider-value');
  const trows = d3.select('tbody').selectAll('tr');
  const formGroup = d3.select('.form-group');
  const hiddenDepth = d3.select('#hidden-depths');


  sliderValue.on('input', () => {
    if (+sliderValue.node().value > +slider.node().max) {
      sliderValue.node().value = slider.node().max;
      slider.node().value = slider.node().max;
    } else {
      slider.node().value = sliderValue.node().value;
    }
    toggleStyles(trows, slider);
    hiddenDepth.attr('value', slider.node().value).dispatch('change');
  });

  sliderValue.on('change', () => {
    if (!sliderValue.node().value) {
      sliderValue.node().value = 0;
      slider.node().value = slider.node().min;
    }
    toggleStyles(trows, slider);
    hiddenDepth.attr('value', slider.node().value).dispatch('change');
  });

  slider.on('input', () => {
    sliderValue.node().value = slider.node().value;
    toggleStyles(trows, slider);
    formGroup.selectAll('span')
      .data(['Sample Loss: ', `${calcSampleLoss(trows, slider).join(' (')}%)`])
      .text(d => d);
    hiddenDepth.attr('value', slider.node().value).dispatch('change');
  });

  formGroup
      .append('div')
    .attr('display', 'inline-block')
    .style('padding-top', '10px')
      .selectAll('span')
    .data(['Sample Loss: ', `${calcSampleLoss(trows, slider).join(' (')}%)`])
      .enter()
    .append('span')
    .text(d => d);

  sliderValue.node().value = slider.node().value;
};

export default initializeControls;
