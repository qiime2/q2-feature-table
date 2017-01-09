import { select } from 'd3';
import { chart } from './chart';

const slider = select('#slider');
const sliderValue = select('#slider-value');

slider.on('input', () => {
  sliderValue.node().value = slider.node().value;
});

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

sliderValue.node().value = slider.node().value;
chart(t, m);
