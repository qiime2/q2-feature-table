// ----------------------------------------------------------------------------
// Copyright (c) 2016-2017, QIIME 2 development team.
//
// Distributed under the terms of the Modified BSD License.
//
// The full license is in the file LICENSE, distributed with this software.
// ----------------------------------------------------------------------------

import * as d3 from 'd3';

let DROPPED = [];


const buildBins = (data, x, isLinear) => (
  isLinear ? d3.histogram().domain(x.domain())(data) : d3.histogram()(data.map(x))
);


const updateChart = (metadata, props, xScale, yScale) => {
  const option = d3.select('select').node().value;
  const validKeys = Object.keys(metadata[option]).filter(key => !DROPPED.includes(key));
  const data = validKeys.map(key => metadata[option][key]);
  const isLinear = !isNaN(data.reduce((x, y) => x + (+y), 0));

  const bins = buildBins(data, xScale, isLinear);

  // TODO: There has to be a better way than transforming each group downwards...
  d3.selectAll('.overlay-group')
    .data(bins)
    .attr('style', function transform(d) {
      const parent = d3.select(this).node().parentNode;
      const barHeight = d3.select(parent)
        .select('rect').node().getBBox().height;

      return `transform: translate(0, ${barHeight - (props.height - yScale(d.length))}px)`;
    })
  .select('.overlay')
    .attr('height', d => props.height - yScale(d.length));
};

export const addSampleMetadata = (sampleID) => {
  if (DROPPED.includes(sampleID)) {
    DROPPED = DROPPED.filter(i => i !== sampleID);
  }
};

export const dropSampleMetadata = (sampleID) => {
  if (!DROPPED.includes(sampleID)) {
    DROPPED.push(sampleID);
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
      .domain([minX - ((maxX - minX) * 0.03), maxX])
      .range([0, props.width]) :
    d3.scaleBand()
      .domain(setArray)
      .rangeRound([0, props.width]);

  const bins = buildBins(data, x, isLinear);

  const y = d3.scaleLinear()
    .domain([0, d3.max(bins, d => d.length)])
    .range([props.height, 0]);

  const bar = g.selectAll('.bar')
    .data(bins)
    .enter()
      .append('g')
        .attr('class', 'bar')
        .attr('transform', d => `translate(${isLinear ? x(d.x0) : d.x0}, ${y(d.length)})`);

  const xAxis = d3.axisBottom(x);
  const yAxis = d3.axisLeft(y);

  const calcBarWidth = (d) => {
    if (isLinear) {
      return Math.floor(x(x.ticks(xAxis.ticks()[0])[1]) - x(x.ticks(xAxis.ticks()[0])[0])) - 1;
    }
    return (d.x0 === d.x1) ? 50 : d.x1 - d.x0 - 1;
  };

  bar.append('rect')
    .attr('x', 1)
    .attr('width', calcBarWidth)
    .attr('height', d => props.height - y(d.length))
  .select(function backToBar() { return this.parentNode; })
    .append('g')
      .attr('class', 'overlay-group')
    .append('rect')
      .attr('class', 'overlay')
      .attr('x', 1)
      .attr('width', calcBarWidth)
      .attr('height', d => props.height - y(d.length));

  g.append('g')
    .attr('class', 'axis axis--x')
    .attr('transform', `translate(0, ${props.height})`)
    .call(xAxis);

  if (setArray.length > 5) {
    g.selectAll('text')
      .attr('y', 0)
      .attr('x', 9)
      .attr('dy', '.35em')
      .attr('transform', 'rotate(90)')
      .style('text-anchor', 'start');
  }


  g.append('g')
    .attr('class', 'axis axis--y')
    .call(yAxis);

  d3.select('#slider')
    .on('input.drop', () => {
      d3.select('tbody')
        .selectAll('tr')
          .each(d => (
            +d[1] < +d3.select('#slider').node().value ?
              dropSampleMetadata(d[0]) :
              addSampleMetadata(d[0])
          ));
      updateChart(metadata, props, x, y);
    });

  // If anything has been dropped already, update on redraw
  updateChart(metadata, props, x, y);
};


const initializeChart = (metadata) => {
  const margin = { top: 10, right: 30, bottom: 30, left: 40 };
  const props = {
    margin,
    width: 960 - margin.left - margin.right,
    height: 500 - margin.top - margin.bottom,
  };

  const svg = d3.select('#histogram')
    .append('svg')
      .attr('width', props.width + props.margin.left + props.margin.right)
      .attr('height', 600 + props.margin.top + props.margin.bottom);
  const select = d3.select('.form-group')
    .append('select')
      .attr('class', 'form-control')
      .on('change', () => buildChart(svg, metadata, props));
  select.selectAll('option')
    .data(Object.keys(metadata))
    .enter()
      .append('option')
      .text(d => d);

  buildChart(svg, metadata, props);
};

export default initializeChart;
