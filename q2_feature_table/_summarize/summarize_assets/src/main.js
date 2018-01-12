// ----------------------------------------------------------------------------
// Copyright (c) 2016-2018, QIIME 2 development team.
//
// Distributed under the terms of the Modified BSD License.
//
// The full license is in the file LICENSE, distributed with this software.
// ----------------------------------------------------------------------------

import initializeControls from './controls';
import initializeTable from './table';
import initializeChart from './chart';

export const init = (metadata, counts) => {
  initializeTable(counts);
  if (Object.keys(metadata).length) {
    initializeChart(metadata, counts);
  }
  initializeControls();
};
