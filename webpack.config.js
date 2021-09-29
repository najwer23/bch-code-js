
const path = require('path')
const HtmlWebpackPlugin = require('html-webpack-plugin');

module.exports = {
  entry: {
    index: path.resolve(__dirname, './assets/js/index.js'),
    style: path.resolve(__dirname, './assets/js/style.js') 
  },
  output: {
    path: path.resolve(__dirname, './build'),
    filename: '[name].bunde.js'
  },
  plugins: [
    new HtmlWebpackPlugin({
      template: path.resolve(__dirname, "index.html"),
      favicon: "./assets/img/f.png",
    })
  ],
  module: {
    rules: [
      {
        test: /\.js?$/,
        use: {
          loader: 'babel-loader',
          options: {
            presets: ['@babel/preset-env']
          }
        }
      },
      {
        test: /\.(woff|woff2|eot|ttf|otf)$/,
        use: ['file-loader']
      },
      {
        test: /\.(sass|scss|css)$/,
        use: ['style-loader', 'css-loader', 'sass-loader']
      },
      {
        test: /\.(?:ico|gif|svg|png|jpg|jpeg)$/i,
        type: 'asset/resource',
      },
    ]
  },
}










// const path = require('path')

// module.exports = {
//   entry: {
//     app: path.resolve(__dirname, './assets/js/index.js'),
//     styleMain: path.resolve(__dirname, './assets/js/style.js'),
//   },
//   output: {
//     path: path.resolve(__dirname, './build'),
//     filename: '[name].js'
//   },
//   module: {
//     rules: [
//       {
//         test: /\.(sass|scss|css)$/,
//         use: ['style-loader', 'css-loader', 'sass-loader']
//       },
//       {
//         test: /\.(woff|woff2|eot|ttf|otf)$/,
//         use: ['file-loader']
//       },
//       {
//         test: /\.(png|svg|jpg|jpeg|gif)$/,
//         use: ['file-loader']
//       }
//     ]
//   },
//   plugins: [
//     new HtmlWebpackPlugin({
//       template: path.resolve(__dirname, 'index.html'),
//       inject: false 
//     })
//   ]
// }

// // https://github.com/webpack/webpack/issues/1732