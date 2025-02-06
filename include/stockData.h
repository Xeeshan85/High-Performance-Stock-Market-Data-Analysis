#ifndef STOCKDATA_H
#define STOCKDATA_H

#pragma once
#include <iostream>
#include <string>
#include <vector>
#include "../include/stockData.h"

struct StockData {
    std::string name;
    std::vector<int> date, volume;
    std::vector<double> open, high, low, close;
};

struct StockMetrics {
    std::vector<double> rollingAverages;
    std::vector<double> upperBand;
    std::vector<double> lowerBand;
    double mean;
    double variance;
    double stddev;
    double skewness;
    double kurtosis;
    std::vector<int> volumeSpikes;
    std::vector<int> crashPoints;
};


#endif