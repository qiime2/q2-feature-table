// ----------------------------------------------------------------------------
// Copyright (c) 2016-2019, QIIME 2 development team.
//
// Distributed under the terms of the Modified BSD License.
//
// The full license is in the file LICENSE, distributed with this software.
// ----------------------------------------------------------------------------

import * as d3 from 'd3';

const initializeTable = (counts) => {
  const table = d3.select('#table')
    .append('table')
    .attr('class', 'table table-hover');

  const myTableData = [];
  Object.keys(counts).forEach(key => myTableData.push([key, counts[key]]));

  table
    .append('thead')
    .append('tr')
      .selectAll('th')
    .data(['Sample ID', 'Sequence Count'])
      .enter()
    .append('th')
    .text(d => d);

  const rows = table
    .append('tbody')
      .selectAll('tr')
    .data(myTableData)
      .enter()
    .append('tr');

  rows
      .selectAll('td')
    .data(d => d)
      .enter()
    .append('td')
    .text(d => d.toLocaleString('en-US'));
};

export default initializeTable;
