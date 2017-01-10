import * as d3 from 'd3';


const init = () => {
  const slider = d3.select('#slider');
  const sliderValue = d3.select('#slider-value');

  sliderValue.on('input', () => {
    if (parseInt(sliderValue.node().value, 10) > parseInt(slider.node().max, 10)) {
      sliderValue.node().value = slider.node().max;
      slider.node().value = slider.node().max;
    } else {
      slider.node().value = sliderValue.node().value;
    }
  });

  sliderValue.on('change', () => {
    if (!sliderValue.node().value) {
      sliderValue.node().value = 0;
      slider.node().value = slider.node().min;
    }
  });


  slider
    .on('input', () => {
      sliderValue.node().value = slider.node().value;
      d3.select('tbody')
        .selectAll('tr')
        .attr('class', d => (+d[1] < slider.node().value ? 'alert-danger' : ''));
    });

  sliderValue.node().value = slider.node().value;
};

export default init;
