# High-Performance Stock Market Data Analysis System
## Overview
This project is a high-performance Stock Market Data Analysis System that processes historical stock data using SIMD (AVX), multithreading (pthreads), and optimized cache efficient data structures. It calculates essential stock metrics such as:
- Moving Averages (SMA, EMA)
- Bollinger Bands
- MACD Indicator
- Stock Price Statistics (Mean, Variance, Skewness, Kurtosis)
- Volume Spikes Detection
- Flash Crash Detection

By leveraging vectorized operations (AVX) and multithreading (pthreads), the system is optimized for low-latency data processing.

## Features
- Efficiently Parsing CSVs
- High-performance calculations using SIMD (AVX)
- Multithreaded execution with pthreads for parallel processing
- Optimized cache utilization for **low-latency** execution.

## Implemented Indicators 📈

1. **Moving Average (MA)**
   $MA_t = \frac{1}{N} \sum_{i=0}^{N-1} P_{t-i}$

2. **Mean Absolute Error (MAE)**
   $MAE = \frac{1}{N} \sum_{i=1}^{N} |P_i - \hat{P}_i|$

3. **Moving Average Convergence Divergence (MACD)**
   $MACD_t = EMA_{12}(P_t) - EMA_{26}(P_t)$

4. **Bollinger Bands**
   $BB_{upper} = MA_t + k \cdot \sigma_t$
   $BB_{lower} = MA_t - k \cdot \sigma_t$

## Installation
### Dependencies

Ensure you have the following installed on your system:
- GCC (G++) with support for AVX & pthreads
- C++17 or later
- CMake (optional but recommended)
#### Clone the Repository
```
git clone https://github.com/Xeeshan85/High-Performance-Stock-Market-Data-Analysis.git
cd High-Performance-Stock-Market-Data-Analysis
```
### Building the Project
**Using CMake**
Create a Build Directory
```
mkdir build
cd build
```
Generate and Compile
```
cmake ..
make
```
Run
```
./stock_analyzer
```

### Usage Example
To analyze stock data, the program:
- Reads stock metadata (`stock_metadata.csv`).
- Loads individual stock files (`AAPL.csv`, `GOOG.csv`, etc.).
- Computes various stock indicators.
- Detects volume spikes and flash crashes.

### License
This project is open-source under the MIT License.
