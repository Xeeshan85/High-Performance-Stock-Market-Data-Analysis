#include <iostream> // For input output
#include <fstream>  // For reading csvs
#include <sstream>  // For string stream
#include <chrono>   // For measuring time
#include <iomanip>
#include <immintrin.h>  // For SIMD
#include <pthread.h>    // For pthreads
#include <unordered_map>
#include <cmath>    // For NaN
// #include <omp.h>
#include <vector>
#include <bits/std_thread.h>
using namespace std;


struct StockData {
    string name;
    vector<int> date, volume;
    vector<double> open, high, low, close;
};

struct StockMetrics {
    vector<double> rollingAverages;
    vector<double> upperBand;
    vector<double> lowerBand;
    double mean;
    double variance;
    double stddev;
    double skewness;
    double kurtosis;
    vector<int> volumeSpikes;
    vector<int> crashPoints;
};

struct ThreadData {
    vector<StockData>* stockDataList;
    vector<StockMetrics>* metricsResults;
    int start;
    int end;
    pthread_mutex_t* mutex;
};


int stringToDate(const string& dateStr) {
    return stoi(dateStr.substr(0, 4)) * 10000 +
           stoi(dateStr.substr(5, 2)) * 100 +
           stoi(dateStr.substr(8, 2));
}

vector<string> readStockMetaData(string fileName) {
    ifstream fileIn;
    fileIn.open(fileName);
    if (!fileIn) {
        cerr << "Error opening file." << endl;
    }

    string line;
    vector<string> stocks;
    getline(fileIn, line); // for header line
    while (getline(fileIn, line)) {
        stringstream ss(line);
        string cell, s;

        getline(ss, cell, ',');
        getline(ss, cell, ',');
        s = cell;
        getline(ss, cell, ',');
        if (cell == "N")
            stocks.push_back(s);
    }
    fileIn.close();

    return stocks;
}

void readStockFile(StockData& stockData) {
    ifstream fileIn(stockData.name);
    if (!fileIn) {
        cerr << "Error opening file: " << stockData.name << endl;
        return;
    }

    string line;
    getline(fileIn, line); // Skip header

    while (getline(fileIn, line)) {
        stringstream ss(line);
        string cell;
        try {
            // Handling date
            getline(ss, cell, ',');
            try {
                stockData.date.push_back(stringToDate(cell));
            } catch (const exception& e) {
                cerr << "Error parsing date in line: " << line << " | " << e.what() << endl;
                continue;
            }

            // convert string to double safely
            auto safeStod = [](const string& str) -> double {
                try {
                    return str.empty() ? NAN : stod(str);
                } catch (...) {
                    return NAN; // Assign NAN on failure
                }
            };

            double opn, hig, lo, clos, vol;
            if (!getline(ss, cell, ',') || isnan(opn = safeStod(cell))) continue;
            if (!getline(ss, cell, ',') || isnan(hig = safeStod(cell))) continue;
            if (!getline(ss, cell, ',') || isnan(lo = safeStod(cell))) continue;
            if (!getline(ss, cell, ',') || isnan(clos = safeStod(cell))) continue;
            if (!getline(ss, cell, ',') || isnan(vol = safeStod(cell))) continue;

            stockData.open.push_back(opn);
            stockData.high.push_back(hig);
            stockData.low.push_back(lo);
            stockData.close.push_back(clos);
            stockData.volume.push_back(vol);

        } catch (const exception& e) {
            cerr << "Unexpected error while parsing line: " << line << " | " << e.what() << endl;
        }
    }

    fileIn.close();
}


void computeRollingAverages(vector<double> &data, vector<double>& avgs, int window_size=5) {
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

void computeBollingerBands(vector<double> &data, vector<double> &upperBand, vector<double> &lowerBand, int window_size=5) {
    int n = data.size();
    if (n < window_size) return;

    vector<double> avgs(n);
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

void computeEMA(vector<double> &data, vector<double>& ema, int period) {
    int n = data.size();
    if (n == 0 || period <= 0) return;

    ema.resize(n);
    double alpha = 2.0/(period + 1);

    ema[0]=data[0];

    for (int i = 0; i < n; i++) {
        ema[i] = alpha * data[i] + (1 - alpha) * ema[i - 1];
    }
}

void computeMACD(vector<double> &data, vector<double> &macd, int shortPeriod, int longPeriod) {
    int n = data.size();
    vector<double> shortEMA, longEMA;

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

// void computeStatistics(const vector<double> &data, double &mean, double &variance, double &stddev, double &skewness, double &kurtosis) {
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
void computeStatistics(const vector<double> &data, double &mean, double &variance, double &stddev, double &skewness, double &kurtosis) {
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

    // Create mean vector
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
    
    // Reduce vectors to scalars
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

void detectVolumeSpikes(const vector<int> &volume, vector<int> &spikeIndices, double threshold=2.0) {
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

void detectFlashCrashes(const vector<double> &closePrices, vector<int> &crashIndices, double dropThreshold = 5.0) {
    int n = closePrices.size();
    if (n < 2) return;

    for (int i = 1; i < n; i++) {
        double percentageDrop = ((closePrices[i - 1] - closePrices[i]) / closePrices[i - 1]) * 100.0;
        if (percentageDrop >= dropThreshold) {
            crashIndices.push_back(i);
        }
    }
}


void* computeMetrics(void* arg) {
    ThreadData* data = (ThreadData*)arg;
    
    for (int i = data->start; i < data->end; i++) {
        StockMetrics metrics;
        
        computeRollingAverages(data->stockDataList->at(i).close, metrics.rollingAverages, 5);
        computeBollingerBands(data->stockDataList->at(i).close, metrics.upperBand, metrics.lowerBand, 5);
        computeStatistics(data->stockDataList->at(i).close, metrics.mean, metrics.variance, 
                         metrics.stddev, metrics.skewness, metrics.kurtosis);
        detectVolumeSpikes(data->stockDataList->at(i).volume, metrics.volumeSpikes);
        detectFlashCrashes(data->stockDataList->at(i).close, metrics.crashPoints, 5.0);

        pthread_mutex_lock(data->mutex);
        data->metricsResults->at(i) = metrics;
        pthread_mutex_unlock(data->mutex);
    }
    return nullptr;
}

vector<StockMetrics> parallelProcessing(vector<StockData>& stockDataList) {
    int num_threads = std::thread::hardware_concurrency();
    vector<pthread_t> threads(num_threads);
    vector<ThreadData> threadData(num_threads);
    vector<StockMetrics> metricsResults(stockDataList.size());
    pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
    
    int chunk_size = stockDataList.size() / num_threads;
    int remainder = stockDataList.size() % num_threads;
    
    int start = 0;
    for (int i = 0; i < num_threads; i++) {
        threadData[i].stockDataList = &stockDataList;
        threadData[i].metricsResults = &metricsResults;
        threadData[i].mutex = &mutex;
        threadData[i].start = start;
        threadData[i].end = start + chunk_size + (i < remainder ? 1 : 0);
        start = threadData[i].end;
        
        pthread_create(&threads[i], nullptr, computeMetrics, (void*)&threadData[i]);
    }
    
    for (int i = 0; i < num_threads; i++) {
        pthread_join(threads[i], nullptr);
    }
    
    pthread_mutex_destroy(&mutex);
    return metricsResults;
}


int main() {
    auto start = chrono::high_resolution_clock::now();

    vector<string> stocks = readStockMetaData("stocks_metadata.csv");

    // ######### Reading Stocks data #############
    int totalFiles = stocks.size();
    vector<StockData> stockDataList(totalFiles);
    unordered_map<string, StockData*> stockDataMap;
    // #pragma omp parallel for
    for (int i = 0; i < totalFiles; i++) {
        stockDataList[i].name = "stocks/" + stocks[i] + ".csv";
        readStockFile(stockDataList[i]);
    
        // #pragma omp critical
        stockDataMap[stocks[i]] = &stockDataList[i];
    }
    cout << stockDataList.size() << " " << totalFiles << endl;


    // Calculating Metrics
    vector<StockMetrics> results = parallelProcessing(stockDataList);
    unordered_map<string, StockMetrics> stockMetricsMap;
    for (int i = 0; i < totalFiles; i++) {
        stockMetricsMap[stocks[i]] = results[i];
    }
    cout << "Metrics size: " << results.size() << endl;

    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed = end - start;
    cout << "Time taken: " << elapsed.count() << " seconds" << endl;

    return 0;
}