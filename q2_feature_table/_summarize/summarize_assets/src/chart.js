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


const categoricalBin = (data, counts, x) => {
  const result = [];
  x.domain().forEach(d => result.push(Object.keys(data).filter(e => data[e] === d).map(id => counts[id])));
  for (let i = 0; i < x.domain().length; i += 1) {
    result[i].x0 = x(x.domain()[i]);
  }
  return result
};


const updateChart = (svg, metadata, props, counts, threshold) => {
  const option = select('select').node().value;

  const data = metadata[option];
  const categories = [...new Set(Object.keys(data).map(key => data[key]))];

  categories.sort(naturalSort);


  svg.selectAll('g').remove();
  svg.selectAll('.bar').remove();
  svg.selectAll('.overlay-group').remove();

  const g = svg
    .append('g')
      .attr('transform', `translate(${props.margin.left}, ${props.margin.top})`);

  const x = scaleBand()
    .domain(categories)
    .range([0, props.width])
    .paddingInner([0.1])
    .paddingOuter([0.3]);

  const bins = categoricalBin(data, counts, x);

  const y = scaleLinear()
    .domain([0, max(bins, d => d.length)])
    .range([props.height, 0]);

  const bar = g.selectAll('.bar')
    .data(bins)
    .enter()
      .append('g')
        .attr('class', 'bar');

  g.append('g')
    .attr('class', 'axis axis--x')
    .attr('transform', `translate(0, ${props.height})`)
    .call(axisBottom(x));

  if (categories.length > 5) {
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
    .attr('width', x.bandwidth() - 1)
    .attr('height', d => props.height - y(d.length));

  selectAll('.bar')
      .data(bins)
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

  g.append('g')
    .attr('class', 'axis axis--y')
    .call(axisLeft(y));
};


const initializeChart = (metadata, counts) => {
  select('svg').remove();
  select('.metadata-dropdown').remove();

  const margin = { top: 10, right: 30, bottom: 30, left: 30 };
  const width = select('#histogram').node().offsetWidth * 0.75;
  const props = {
    margin,
    width: width - margin.left - margin.right,
    height: ((width * 9) / 16) - margin.top - margin.bottom,
  };
  let threshold = 0;

  const svg = select('#histogram')
    .append('svg')
      .attr('width', props.width + props.margin.left + props.margin.right)
      .attr('height', props.height + props.margin.top + props.margin.bottom + 150);
  const selection = select('.form-group')
    .append('select')
      .attr('class', 'form-control metadata-dropdown')
      .on('change', () => updateChart(svg, metadata, props, counts, threshold));
  selection.selectAll('option')
    .data(Object.keys(metadata))
    .enter()
      .append('option')
      .text(d => d);

  const callUpdate = () => {
    threshold = +select('#slider').node().value;
    updateChart(svg, metadata, props, counts, threshold);
  };

  select('#slider')
    .on('input.drop', callUpdate);
  select('#slider-value')
    .on('change.drop', callUpdate);
  select('#slider-value')
    .on('input.type', callUpdate);

  updateChart(svg, metadata, props, counts, threshold);
  svg.attr('display', 'block').style('margin', '0 auto');
};

export default initializeChart;
