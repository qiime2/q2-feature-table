webpackJsonp([0],[
/* 0 */
/***/ function(module, exports, __webpack_require__) {

	'use strict';

	var _d = __webpack_require__(1);

	var _chart = __webpack_require__(2);

	var slider = (0, _d.select)('#slider');
	var sliderValue = (0, _d.select)('#slider-value');

	slider.on('input', function () {
	  sliderValue.node().value = slider.node().value;
	});

	sliderValue.on('input', function () {
	  if (parseInt(sliderValue.node().value, 10) > parseInt(slider.node().max, 10)) {
	    sliderValue.node().value = slider.node().max;
	    slider.node().value = slider.node().max;
	  } else {
	    slider.node().value = sliderValue.node().value;
	  }
	});

	sliderValue.on('change', function () {
	  if (!sliderValue.node().value) {
	    sliderValue.node().value = 0;
	    slider.node().value = slider.node().min;
	  }
	});

	sliderValue.node().value = slider.node().value;
	(0, _chart.chart)(t, m);

/***/ },
/* 1 */,
/* 2 */
/***/ function(module, exports, __webpack_require__) {

	'use strict';

	Object.defineProperty(exports, "__esModule", {
	  value: true
	});
	exports.chart = undefined;

	var _d = __webpack_require__(1);

	var d3 = _interopRequireWildcard(_d);

	function _interopRequireWildcard(obj) { if (obj && obj.__esModule) { return obj; } else { var newObj = {}; if (obj != null) { for (var key in obj) { if (Object.prototype.hasOwnProperty.call(obj, key)) newObj[key] = obj[key]; } } newObj.default = obj; return newObj; } }

	var chart = exports.chart = function chart(table, metadata) {
	  console.log(table);
	  console.log(metadata);
	};

/***/ }
]);