// ----------------------------------------------------------------------------
// Copyright (c) 2016-2017, QIIME 2 development team.
//
// Distributed under the terms of the Modified BSD License.
//
// The full license is in the file LICENSE, distributed with this software.
// ----------------------------------------------------------------------------

import {
  axisBottom,
  axisLeft,
  scaleBand,
  scaleLinear,
  histogram,
  max,
  select,
  selectAll,
} from 'd3';
import naturalSort from 'javascript-natural-sort';

let DROPPED = [];

const buildBins = (data, x) => (histogram()(data.map(x)));


const addSampleMetadata = (sampleID) => {
  if (DROPPED.includes(sampleID)) {
    DROPPED = DROPPED.filter(i => i !== sampleID);
  }
};


const dropSampleMetadata = (sampleID) => {
  if (!DROPPED.includes(sampleID)) {
    DROPPED.push(sampleID);
  }
};


const updateChart = (svg, metadata, props) => {
  const option = select('select').node().value;
  const validKeys = Object.keys(metadata[option]).filter(key => !DROPPED.includes(key));
  const validData = validKeys.map(key => metadata[option][key]);
  const data = Object.keys(metadata[option]).map(key => metadata[option][key]);
  const setArray = [...new Set(data)];

  setArray.sort(naturalSort);

  svg.selectAll('g').remove();
  svg.selectAll('.bar').remove();
  svg.selectAll('.overlay-group').remove();

  const g = svg
    .append('g')
      .attr('transform', `translate(${props.margin.left}, ${props.margin.top})`);

  const x = scaleBand()
    .domain(setArray)
    .range([0, props.width])
    .paddingInner([0.1])
    .paddingOuter([0.3]);

  const originalBins = buildBins(data, x);
  const newBins = buildBins(validData, x);

  const y = scaleLinear()
    .domain([0, Math.max(max(originalBins, d => d.length), max(newBins, d => d.length))])
    .range([props.height, 0]);

  const bar = g.selectAll('.bar')
    .data(originalBins)
    .enter()
      .append('g')
        .attr('class', 'bar');

  g.append('g')
    .attr('class', 'axis axis--x')
    .attr('transform', `translate(0, ${props.height})`)
    .call(axisBottom(x));

  if (setArray.length > 5) {
    g.selectAll('text')
      .attr('y', 0)
      .attr('x', 9)
      .attr('dy', '.35em')
      .attr('transform', 'rotate(90)')
      .style('text-anchor', 'start');
  }

  bar.append('rect')
    .attr('class', 'original')
    .attr('fill', 'black')
    .attr('opacity', 0.15)
    .attr('x', d => d.x0)
    .attr('y', d => y(d.length))
    .attr('width', d => (d.x1 - d.x0 - 1 > 0 ? d.x1 - d.x0 - 1 : x.bandwidth() - 1))
    .attr('height', d => props.height - y(d.length));

  selectAll('.bar')
      .data(newBins)
    .append('g')
      .attr('class', 'overlay-group')
    .append('rect')
      .attr('class', 'overlay')
      .attr('fill', 'steelblue')
      .style('opacity', 1)
      .style('fill-opacity', 1)
      .attr('stroke', 'lightgray')
      .attr('x', d => d.x0)
      .attr('y', d => y(d.length))
      .attr('width', d => (d.x1 - d.x0 - 1 > 0 ? d.x1 - d.x0 - 1 : x.bandwidth() - 1))
      .attr('height', d => props.height - y(d.length));

  g.append('g')
    .attr('class', 'axis axis--y')
    .call(axisLeft(y));
};


const initializeChart = (metadata) => {
  const margin = { top: 10, right: 30, bottom: 30, left: 40 };
  const props = {
    margin,
    width: 960 - margin.left - margin.right,
    height: 500 - margin.top - margin.bottom,
  };

  const svg = select('#histogram')
    .append('svg')
      .attr('width', props.width + props.margin.left + props.margin.right)
      .attr('height', 600 + props.margin.top + props.margin.bottom);
  const selection = select('.form-group')
    .append('select')
      .attr('class', 'form-control')
      .on('change', () => updateChart(svg, metadata, props));
  selection.selectAll('option')
    .data(Object.keys(metadata))
    .enter()
      .append('option')
      .text(d => d);

  const callUpdate = () => {
    select('tbody')
      .selectAll('tr')
        .each(d => (
          +d[1] < +select('#slider').node().value ?
            dropSampleMetadata(d[0]) :
            addSampleMetadata(d[0])
        ));
    updateChart(svg, metadata, props);
  };

  select('#slider')
    .on('input.drop', callUpdate);
  select('#slider-value')
    .on('input.type', callUpdate);

  updateChart(svg, metadata, props);
};

export default initializeChart;
