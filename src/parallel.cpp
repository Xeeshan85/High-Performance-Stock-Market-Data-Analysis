#include "../include/parallel.h"
#include "../include/analysis.h"
#include <pthread.h>
#include <bits/std_thread.h>

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

std::vector<StockMetrics> parallelProcessing(std::vector<StockData>& stockDataList) {
    int num_threads = std::thread::hardware_concurrency();
    std::vector<pthread_t> threads(num_threads);
    std::vector<ThreadData> threadData(num_threads);
    std::vector<StockMetrics> metricsResults(stockDataList.size());
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
