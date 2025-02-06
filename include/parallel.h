#ifndef PARALLEL_H
#define PARALLEL_H

#pragma once
#include <vector>
#include <pthread.h>
#include "../include/stockData.h"

struct ThreadData {
    std::vector<StockData>* stockDataList;
    std::vector<StockMetrics>* metricsResults;
    int start;
    int end;
    pthread_mutex_t* mutex;
};

void* computeMetrics(void* arg);
std::vector<StockMetrics> parallelProcessing(std::vector<StockData>& stockDataList);

#endif