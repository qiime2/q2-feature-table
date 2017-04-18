// ----------------------------------------------------------------------------
// Copyright (c) 2016-2017, QIIME 2 development team.
//
// Distributed under the terms of the Modified BSD License.
//
// The full license is in the file LICENSE, distributed with this software.
// ----------------------------------------------------------------------------

import naturalSort from 'javascript-natural-sort';
import * as d3 from 'd3';


const categoricalBin = (data, counts, x) => {
  const result = [];
  x.domain().forEach((cat) => {
    const catArray = Object.keys(data)
                      .filter(sampleID => data[sampleID] === cat)
                      .map(sampleID => counts[sampleID]);
    catArray.x0 = x(cat);
    result.push(catArray);
  });
  return result;
};


const updateChart = (bar, props, y, threshold) => {
  bar.selectAll('.overlay')
    .attr('y', d => y(d.filter(count => count >= threshold).length))
    .attr('height', d => props.height - y(d.filter(count => count >= threshold).length));
};

const buildChart = (svg, props, metadata, counts) => {
  svg.selectAll('*').remove();

  const threshold = d3.select('#slider').node().min;
  const option = d3.select('select').node().value;
  const data = metadata[option];
  const categories = [...new Set(Object.keys(data).map(key => data[key]))];

  categories.sort(naturalSort);

  svg
    .append('text')
    .attr('transform', 'rotate(-90)')
    .attr('x', 0 - (props.height / 2))
    .attr('y', 15)
    .attr('dy', '0em')
    .attr('font-size', '14px')
    .style('text-anchor', 'middle')
    .text('Number of Samples');

  const g = svg
      .append('g')
    .attr('transform', `translate(${props.margin.left}, ${props.margin.top})`);

  const x = d3
    .scaleBand()
    .domain(categories)
    .range([props.margin.left, props.width])
    .paddingInner([0.1])
    .paddingOuter([0.3]);

  const bins = categoricalBin(data, counts, x);

  const y = d3
    .scaleLinear()
    .domain([0, d3.max(bins, d => d.length)])
    .range([props.height, 0]);

  g.append('g')
    .attr('class', 'axis axis--x')
    .attr('transform', `translate(0, ${props.height})`)
    .call(d3.axisBottom(x));

  // vertical axis labels if there are a good amount of categories
  if (categories.length > 5) {
    g.selectAll('text')
      .attr('y', 0)
      .attr('x', 9)
      .attr('dy', '.35em')
      .attr('transform', 'rotate(90)')
      .style('text-anchor', 'start');
  }

  g.append('g')
    .attr('class', 'axis axis--y')
    .attr('transform', `translate(${props.margin.left}, 0)`)
    .call(d3.axisLeft(y));

  const bar = g
      .selectAll('.bar')
    .data(bins)
      .enter()
    .append('g')
    .attr('class', 'bar');

  // build backdrop
  bar
      .append('rect')
    .attr('class', 'original')
    .attr('fill', 'black')
    .attr('opacity', 0.15)
    .attr('x', d => d.x0)
    .attr('y', d => y(d.length))
    .attr('width', x.bandwidth() - 1)
    .attr('height', d => props.height - y(d.length));

  // build overlay
  bar
      .append('g')
    .attr('class', 'overlay-group')
      .append('rect')
    .attr('class', 'overlay')
    .attr('fill', 'steelblue')
    .style('opacity', 1)
    .style('fill-opacity', 1)
    .attr('stroke', 'lightgray')
    .attr('x', d => d.x0)
    .attr('y', d => y(d.filter(count => count >= threshold).length))
    .attr('width', x.bandwidth() - 1)
    .attr('height', d => props.height - y(d.filter(count => count >= threshold).length));

  const hiddenDepths = d3.select('#hidden-depths');
  hiddenDepths
    .on('change', () => updateChart(bar, props, y, +hiddenDepths.node().value));
  updateChart(bar, props, y, +hiddenDepths.node().value);
};


const initializeChart = (metadata, counts) => {
  d3.select('svg').remove();
  d3.select('.metadata-dropdown').remove();

  const margin = { top: 10, right: 30, bottom: 30, left: 30 };
  const width = d3.select('#histogram').node().offsetWidth * 0.75;
  const props = {
    margin,
    width: width - margin.left - margin.right,
    height: ((width * 9) / 16) - margin.top - margin.bottom };

  const svg = d3.select('#histogram')
      .append('svg')
    .attr('width', props.width + props.margin.left + props.margin.right)
    .attr('height', props.height + props.margin.top + props.margin.bottom + 150);

  const formGroup = d3.select('.form-group');
  const selection = formGroup
      .append('select')
    .attr('class', 'form-control metadata-dropdown')
    .on('change', () => buildChart(svg, props, metadata, counts));
  selection.selectAll('option')
    .data(Object.keys(metadata))
      .enter()
    .append('option')
    .text(d => d);
  formGroup
      .append('input')
    .attr('type', 'hidden')
    .attr('value', d3.select('#slider').node().value)
    .attr('id', 'hidden-depths');

  buildChart(svg, props, metadata, counts);
  svg.attr('display', 'block').style('margin', '0 auto');
};

export default initializeChart;
