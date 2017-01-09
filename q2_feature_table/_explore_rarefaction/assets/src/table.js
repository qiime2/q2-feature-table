import * as d3 from 'd3';

const chart = (data, metadata, counts) => {
  const table = d3.select('#table').append('table');
  const header = table.append('thead').append('tr');
  const tableBody = table.append('tbody');

  const myTableData = [];
  Object.keys(counts).forEach(key => myTableData.push([key, counts[key]]));

  table
    .attr('class', 'table table-striped table-hover');
  header
    .selectAll('th')
    .data(['Sample ID', 'Feature Count'])
    .enter()
    .append('th')
    .text(d => d);
  const rows = tableBody
    .selectAll('tr')
    .data(myTableData)
    .enter()
    .append('tr');
  const cells = rows  // eslint-disable-line no-unused-vars
    .selectAll('td')
    .data(d => d)
    .enter()
    .append('td')
    .text(d => d);
};

export default chart;
