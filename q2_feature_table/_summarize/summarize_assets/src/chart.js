// ----------------------------------------------------------------------------
// Copyright (c) 2016-2017, QIIME 2 development team.
//
// Distributed under the terms of the Modified BSD License.
//
// The full license is in the file LICENSE, distributed with this software.
// ----------------------------------------------------------------------------

import {
  select,
  selectAll,
  histogram,
  scaleBand,
  scaleLinear,
  max,
  axisBottom,
  axisLeft,
} from 'd3';

let DROPPED = [];

const buildBins = (data, x) => (histogram()(data.map(x)));


const updateChart = (metadata, props, xScale, yScale) => {
  const option = select('select').node().value;
  const validKeys = Object.keys(metadata[option]).filter(key => !DROPPED.includes(key));
  const data = validKeys.map(key => metadata[option][key]);

  const bins = buildBins(data, xScale);

  selectAll('.overlay-group')
    .data(bins)
  .select('.overlay')
    .attr('y', d => yScale(d.length))
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
  const option = select('select').node().value;

  const data = Object.keys(metadata[option]).map(key => metadata[option][key]);
  const setArray = [...new Set(data)];

  svg.selectAll('g').remove();
  svg.selectAll('rect').remove();

  const g = svg
    .append('g')
      .attr('transform', `translate(${props.margin.left}, ${props.margin.top})`);

  const x = scaleBand()
    .domain(setArray)
    .rangeRound([0, props.width - props.margin.left - props.margin.right]);

  const bins = buildBins(data, x);

  const y = scaleLinear()
    .domain([0, max(bins, d => d.length)])
    .range([props.height, 0]);

  const bar = g.selectAll('.bar')
    .data(bins)
    .enter()
      .append('g')
        .attr('class', 'bar');
        // .attr('transform', d => `translate(${d.x0}, ${y(d.length)})`);

  const xAxis = axisBottom(x);
  const yAxis = axisLeft(y);

  bar.append('rect')
    .attr('fill', 'black')
    .attr('opacity', 0.15)
    .attr('x', d => d.x0)
    .attr('y', d => y(d.length))
    .attr('width', x.bandwidth() - 1)
    .attr('height', d => props.height - y(d.length))
  .select(function backToBar() { return this.parentNode; })
    .append('g')
      .attr('class', 'overlay-group')
    .append('rect')
      .attr('class', 'overlay')
      .attr('fill', 'steelblue')
      .attr('opacity', 1)
      .attr('stroke', 'lightgray')
      .attr('x', d => d.x0)
      .attr('y', d => y(d.length))
      .attr('width', x.bandwidth() - 1)
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

  select('#slider')
    .on('input.drop', () => {
      select('tbody')
        .selectAll('tr')
          .each(d => (
            +d[1] < +select('#slider').node().value ?
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

  const svg = select('#histogram')
    .append('svg')
      .attr('width', props.width + props.margin.left + props.margin.right)
      .attr('height', 600 + props.margin.top + props.margin.bottom);
  const selection = select('.form-group')
    .append('select')
      .attr('class', 'form-control')
      .on('change', () => buildChart(svg, metadata, props));
  selection.selectAll('option')
    .data(Object.keys(metadata))
    .enter()
      .append('option')
      .text(d => d);

  buildChart(svg, metadata, props);
};

export default initializeChart;
