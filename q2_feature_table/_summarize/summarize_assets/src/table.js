// ----------------------------------------------------------------------------
// Copyright (c) 2016-2017, QIIME 2 development team.
//
// Distributed under the terms of the Modified BSD License.
//
// The full license is in the file LICENSE, distributed with this software.
// ----------------------------------------------------------------------------

import * as d3 from 'd3';

export const myTableData = [];

const initializeTable = (counts) => {
  const table = d3.select('#table').append('table');
  const header = table.append('thead').append('tr');
  const tableBody = table.append('tbody');
  const formGroup = d3.select('.form-group').append('span');

  Object.keys(counts).forEach(key => myTableData.push([key, counts[key]]));

  table.attr('class', 'table table-hover');
  header.selectAll('th')
    .data(['Sample ID', 'Feature Count'])
    .enter()
      .append('th')
        .text(d => d);
  const rows = tableBody.selectAll('tr')
    .data(myTableData)
    .enter()
      .append('tr');
  rows.selectAll('td')
    .data(d => d)
    .enter()
      .append('td')
      .text(d => d);


  const calcSampleLoss = () => {
    const lost = myTableData.reduce((i, j) => i + (j[1] < +d3.select('#slider').node().value ? 1 : 0), 0);
    return [lost, ` (${((lost / myTableData.length) * 100).toFixed(2)}%)`];
  };

  formGroup
    .append('div')
  .selectAll('span')
    .data(['Sample Loss: ', ...calcSampleLoss()])
    .enter()
      .append('span')
      .text(d => d);

  d3.select('#slider')
    .on('input.calc', () => {
      formGroup.selectAll('span')
        .data(['Sample Loss: ', ...calcSampleLoss()])
        .text(d => d);
    });
};

export default initializeTable;
