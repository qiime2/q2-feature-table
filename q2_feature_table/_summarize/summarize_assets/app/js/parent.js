// ----------------------------------------------------------------------------
// Copyright (c) 2016-2017, QIIME 2 development team.
//
// Distributed under the terms of the Modified BSD License.
//
// The full license is in the file LICENSE, distributed with this software.
// ----------------------------------------------------------------------------

var frame;
var frameSrc = {
  'overview-tab': './overview.html',
  'sample-detail-tab': './sample-frequency-detail.html',
  'feature-detail-tab': './feature-frequency-detail.html'
};

function frameLoad(event) {
  frame.height = `${event.data + 50}px`;
}

function toggleClass() {
  if (document.querySelector('.active').id != this.id) {
    document.querySelector('.active').className = '';
    this.className = 'active';
    frame.height = '0px';
    frame.src = frameSrc[this.id];
  }
}

function init() {
  var overviewTab = document.getElementById('overview-tab');
  var sampleTab = document.getElementById('sample-detail-tab');
  var featureTab = document.getElementById('feature-detail-tab');
  frame = document.getElementById('tab-frame');

  overviewTab.addEventListener('click', toggleClass);
  sampleTab.addEventListener('click', toggleClass);
  featureTab.addEventListener('click', toggleClass);
}

document.addEventListener('DOMContentLoaded', init);
window.addEventListener('message', frameLoad);
