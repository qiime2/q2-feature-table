// ----------------------------------------------------------------------------
// Copyright (c) 2016-2017, QIIME 2 development team.
//
// Distributed under the terms of the Modified BSD License.
//
// The full license is in the file LICENSE, distributed with this software.
// ----------------------------------------------------------------------------

import * as d3 from 'd3';


const initializeControls = () => {
  const slider = d3.select('#slider');
  const sliderValue = d3.select('#slider-value');

  sliderValue.on('input', () => {
    if (+sliderValue.node().value > +slider.node().max) {
      sliderValue.node().value = slider.node().max;
      slider.node().value = slider.node().max;
    } else {
      slider.node().value = sliderValue.node().value;
    }
  });

  sliderValue.on('change', () => {
    if (!sliderValue.node().value) {
      sliderValue.node().value = 0;
      slider.node().value = slider.node().min;
    }
  });

  slider
    .on('input.slide', () => {
      sliderValue.node().value = slider.node().value;
      d3.select('tbody')
        .selectAll('tr')
        .attr('class', d => (+d[1] < +slider.node().value ? 'alert-danger' : ''));
    });

  sliderValue.node().value = slider.node().value;
};

export default initializeControls;
