// ----------------------------------------------------------------------------
// Copyright (c) 2016-2017, QIIME 2 development team.
//
// Distributed under the terms of the Modified BSD License.
//
// The full license is in the file LICENSE, distributed with this software.
// ----------------------------------------------------------------------------

import initializeControls from './controls';
import initializeTable from './table';
import initializeHistogram from './hist';

const init = (table, metadata, counts) => {
  initializeControls();
  initializeTable(counts);
  initializeHistogram(metadata);
};

export default init;
