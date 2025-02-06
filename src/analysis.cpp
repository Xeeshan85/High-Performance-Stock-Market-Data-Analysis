#include "../include/analysis.h"
#include <immintrin.h>
#include <cmath>


void computeRollingAverages(std::vector<double> &data, std::vector<double>& avgs, int window_size) {
    int n = data.size();
    if (n < window_size) return;

    avgs.resize(n);
    
    // Initial window sum
    double window_sum = 0;
    for (int i = 0; i < window_size; i++) {
        window_sum += data[i];
    }
    
    // First window average
    avgs[window_size - 1] = window_sum / window_size;

    for (int i = window_size; i < n; i++) {
        window_sum = (window_sum - data[i - window_size]) + data[i];
        avgs[i] = window_sum / window_size;
    }
}

void computeBollingerBands(std::vector<double> &data, std::vector<double> &upperBand, std::vector<double> &lowerBand, int window_size) {
    int n = data.size();
    if (n < window_size) return;

    std::vector<double> avgs(n);
    upperBand.resize(n); lowerBand.resize(n);
    
    double window_sum = 0, window_sum_sq = 0;
    for (int i = 0; i < window_size; i++) {
        window_sum += data[i];
        window_sum_sq += data[i] * data[i];
    }
    
    double avg = window_sum / window_size;
    double variance = (window_sum_sq / window_size) - (avg * avg);
    double stddev = sqrt(variance);
    
    avgs[window_size - 1] = avg;
    upperBand[window_size - 1] = avg + (2 * stddev);
    lowerBand[window_size - 1] = avg - (2 * stddev);
    
    // Slide window
    for (int i = window_size; i < n; i++) {
        window_sum = (window_sum - data[i - window_size]) + data[i];
        window_sum_sq = (window_sum_sq - (data[i - window_size] * data[i - window_size])) + (data[i] * data[i]);
        
        avg = window_sum / window_size;
        variance = (window_sum_sq / window_size) - (avg * avg);
        stddev = sqrt(variance);
        
        avgs[i] = avg;
        upperBand[i] = avg + (2 * stddev);
        lowerBand[i] = avg - (2 * stddev);
    }
}

void computeEMA(std::vector<double> &data, std::vector<double>& ema, int period) {
    int n = data.size();
    if (n == 0 || period <= 0) return;

    ema.resize(n);
    double alpha = 2.0/(period + 1);

    ema[0]=data[0];

    for (int i = 0; i < n; i++) {
        ema[i] = alpha * data[i] + (1 - alpha) * ema[i - 1];
    }
}

void computeMACD(std::vector<double> &data, std::vector<double> &macd, int shortPeriod, int longPeriod) {
    int n = data.size();
    std::vector<double> shortEMA, longEMA;

    computeEMA(data, shortEMA, shortPeriod);
    computeEMA(data, longEMA, longPeriod);

    macd.resize(n, 0.0);
    int width = 4;
    int end = n - (n % width); // Alligning End

    for (int i = 0; i < end; i+= width) {
        __m256d shortVec = _mm256_loadu_pd(&shortEMA[i]);
        __m256d longVec = _mm256_loadu_pd(&longEMA[i]);

        __m256d res = _mm256_sub_pd(shortVec, longVec); // Subtract
        _mm256_storeu_pd(&macd[i], res); // Store result
    }

    for (int i = end; i < n; i++) { // Remaining elements
        macd[i] = shortEMA[i] - longEMA[i];
    }
}

// Without SIMD

// void computeStatistics(const std::vector<double> &data, double &mean, double &variance, double &stddev, double &skewness, double &kurtosis) {
//     int n = data.size();
//     if (n < 2) return;
//     // Mean
//     double sum = 0;
//     for (double val : data) sum += val;
//     mean = sum / n;
//     // Variance & std
//     double sumSq = 0;
//     for (double val : data) sumSq += (val - mean) * (val - mean);
//     variance = sumSq / (n - 1);
//     stddev = sqrt(variance);
//     // Skewness & kurtosis
//     double sumCubed = 0, sumQuartic = 0;
//     for (double val : data) {
//         double diff = val - mean;
//         sumCubed += diff * diff * diff;
//         sumQuartic += diff * diff * diff * diff;
//     }
//     skewness = (sumCubed / n) / (stddev * stddev * stddev);
//     kurtosis = (sumQuartic / n) / (variance * variance) - 3.0; // Excess Kurtosis
// }

// With SIMD; Brought program time from 3.16% to 0.63%
void computeStatistics(const std::vector<double> &data, double &mean, double &variance, double &stddev, double &skewness, double &kurtosis) {
    int n = data.size();
    if (n < 2) return;

    const int width = 4;
    int end = n - (n % width);

    // Mean
    __m256d sum_vec = _mm256_setzero_pd();
    for (int i = 0; i < end; i += width) {
        __m256d data_vec = _mm256_loadu_pd(&data[i]);
        sum_vec = _mm256_add_pd(sum_vec, data_vec);
    }

    double sum = 0, temp[4];
    _mm256_storeu_pd(temp, sum_vec);
    for (int i = 0; i < 4; i++) sum += temp[i];
    
    // remaining elements
    for (int i = end; i < n; i++) {
        sum += data[i];
    }
    mean = sum / n;

    // Create mean std::vector
    __m256d mean_vec = _mm256_set1_pd(mean);
    
    // variance
    __m256d sum_sq_vec = _mm256_setzero_pd();
    __m256d sum_cube_vec = _mm256_setzero_pd();
    __m256d sum_quart_vec = _mm256_setzero_pd();
    
    for (int i = 0; i < end; i += width) {
        __m256d data_vec = _mm256_loadu_pd(&data[i]);
        __m256d diff_vec = _mm256_sub_pd(data_vec, mean_vec);
        
        // squares
        __m256d sq_vec = _mm256_mul_pd(diff_vec, diff_vec);
        sum_sq_vec = _mm256_add_pd(sum_sq_vec, sq_vec);
        
        // cubes
        __m256d cube_vec = _mm256_mul_pd(sq_vec, diff_vec);
        sum_cube_vec = _mm256_add_pd(sum_cube_vec, cube_vec);
        
        // 4rth power
        __m256d quart_vec = _mm256_mul_pd(sq_vec, sq_vec);
        sum_quart_vec = _mm256_add_pd(sum_quart_vec, quart_vec);
    }
    
    // Reduce std::vectors to scalars
    double sumSq = 0, sumCubed = 0, sumQuartic = 0;
    _mm256_storeu_pd(temp, sum_sq_vec);
    for (int i = 0; i < 4; i++) sumSq += temp[i];
    
    _mm256_storeu_pd(temp, sum_cube_vec);
    for (int i = 0; i < 4; i++) sumCubed += temp[i];
    
    _mm256_storeu_pd(temp, sum_quart_vec);
    for (int i = 0; i < 4; i++) sumQuartic += temp[i];
    
    // Handle remaining elements
    for (int i = end; i < n; i++) {
        double diff = data[i] - mean;
        double sq = diff * diff;
        sumSq += sq;
        sumCubed += sq * diff;
        sumQuartic += sq * sq;
    }

    // Calculate final statistics
    variance = sumSq / (n - 1);
    stddev = sqrt(variance);
    skewness = (sumCubed / n) / (stddev * stddev * stddev);
    kurtosis = (sumQuartic / n) / (variance * variance) - 3.0;
}

void detectVolumeSpikes(const std::vector<int> &volume, std::vector<int> &spikeIndices, double threshold) {
    int n = volume.size();
    if (n < 2) return;

    double mean = 0, stddev = 0;
    for (int v : volume) mean += v;
    mean /= n;
    for (int v : volume) stddev += (v - mean) * (v - mean);
    stddev = sqrt(stddev / (n - 1));

    // Identify spikes
    for (int i = 0; i < n; i++) {
        double z_score = (volume[i] - mean) / stddev;
        if (z_score > threshold) spikeIndices.push_back(i);
    }
}

void detectFlashCrashes(const std::vector<double> &closePrices, std::vector<int> &crashIndices, double dropThreshold) {
    int n = closePrices.size();
    if (n < 2) return;

    for (int i = 1; i < n; i++) {
        double percentageDrop = ((closePrices[i - 1] - closePrices[i]) / closePrices[i - 1]) * 100.0;
        if (percentageDrop >= dropThreshold) {
            crashIndices.push_back(i);
        }
    }
}
