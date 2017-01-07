var path = require('path');
var webpack = require('webpack');

module.exports = {
  entry: {
    app: './src/main.js',
    // vendor: ['thenby', 'd3', 'd3-scale-chromatic', 'natural-sort']
  },
  plugins: [
    new webpack.optimize.CommonsChunkPlugin('vendor', 'app/vendor.bundle.js'),
    new webpack.optimize.UglifyJsPlugin({
      compress: { warnings: false }
    }),
    new webpack.NoErrorsPlugin(),
  ],
  output: {
    path: __dirname,
    filename: 'app/bundle.js'
  },
  module: {
    loaders: [
      {
        loader: 'babel-loader',
        exclude: /node_modules/,
        query: {
          presets: ['es2015']
        }
      }
    ]
  },
};
