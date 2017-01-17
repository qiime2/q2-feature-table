// ----------------------------------------------------------------------------
// Copyright (c) 2016-2017, QIIME 2 development team.
//
// Distributed under the terms of the Modified BSD License.
//
// The full license is in the file LICENSE, distributed with this software.
// ----------------------------------------------------------------------------

import * as d3 from 'd3';

// TODO: I don't like this...
let DROPPED = [];
let METADATA;
let PROPS;
let X;
let Y;

const updateChart = () => {
  const option = d3.select('select').node().value;
  const validKeys = Object.keys(METADATA[option]).filter(key => !DROPPED.includes(key));
  const data = validKeys.map(key => METADATA[option][key]);

  const bins = d3.histogram()
    .domain(X.domain())(data);

  // TODO: There has to be a better way than transforming each group downwards...
  d3.selectAll('.overlay-group')
    .data(bins)
    .attr('style', function transform(d) {
      const barHeight = d3.select(d3.select(this).node().parentNode).select('rect').node().getBBox().height;
      return `transform: translate(0, ${barHeight - (PROPS.height - Y(d.length))}px)`;
    })
    .select('.overlay')
    .attr('height', d => PROPS.height - Y(d.length));
};

export const addSampleMetadata = (sampleID) => {
  if (DROPPED.includes(sampleID)) {
    DROPPED = DROPPED.filter(i => i !== sampleID);
    updateChart();
  }
};

export const dropSampleMetadata = (sampleID) => {
  if (!DROPPED.includes(sampleID)) {
    DROPPED.push(sampleID);
    updateChart();
  }
};

const buildChart = (svg, metadata, props) => {
  const option = d3.select('select').node().value;

  const data = Object.keys(metadata[option]).map(key => metadata[option][key]);
  const isLinear = !isNaN(data.reduce((x, y) => x + (+y), 0));
  const setArray = [...new Set(data)];

  svg.selectAll('g').remove();
  svg.selectAll('rect').remove();

  const g = svg
    .append('g')
    .attr('transform', `translate(${props.margin.left}, ${props.margin.top})`);

  const minX = Math.min(...data);
  const maxX = Math.max(...data);

  const x = isLinear ?
    d3.scaleLinear()
      .domain([minX - ((maxX - minX) * 0.03), maxX]) :
    d3.scaleOrdinal()
      .domain(setArray);

  x.range([0, props.width]);

  const bins = d3.histogram()
    .domain(x.domain())(data);

  const y = d3.scaleLinear()
    .domain([0, d3.max(bins, d => d.length)])
    .range([props.height, 0]);

  const bar = g.selectAll('.bar')
    .data(bins)
    .enter()
      .append('g')
      .attr('class', 'bar')
      .attr('transform', d => `translate(${x(d.x0)}, ${y(d.length)})`);

  bar
    .append('rect')
    .attr('x', 1)
    .attr('width', x(bins[0].x1) - x(bins[0].x0) - 1)
    .attr('height', d => props.height - y(d.length));

  bar
    .append('g')
    .attr('class', 'overlay-group')
    .append('rect')
    .attr('class', 'overlay')
    .attr('x', 1)
    .attr('width', x(bins[0].x1) - x(bins[0].x0) - 1)
    .attr('height', d => props.height - y(d.length));

  g
    .append('g')
    .attr('class', 'axis axis--x')
    .attr('transform', `translate(0, ${props.height})`)
    .call(d3.axisBottom(x));

  g
    .append('g')
    .attr('class', 'axis axis--y')
    .call(d3.axisLeft(y));

  X = x;
  Y = y;
};

const initializeChart = (metadata) => {
  METADATA = metadata;
  const margin = { top: 10, right: 30, bottom: 30, left: 40 };
  const props = {
    margin,
    width: 960 - margin.left - margin.right,
    height: 500 - margin.top - margin.bottom,
  };
  PROPS = props;

  const svg = d3.select('#histogram')
    .append('svg')
    .attr('width', props.width + props.margin.left + props.margin.right)
    .attr('height', props.height + props.margin.top + props.margin.bottom);
  const select = d3.select('.form-group')
    .append('select')
    .attr('class', 'form-control')
    .on('change', () => buildChart(svg, metadata, props));
  select
    .selectAll('option')
    .data(Object.keys(metadata))
    .enter()
      .append('option')
      .text(d => d);

  buildChart(svg, metadata, props);
};

export default initializeChart;
