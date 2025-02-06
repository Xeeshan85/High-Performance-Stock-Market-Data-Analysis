#ifndef ANALYSIS_H
#define ANALYSIS_H

#pragma once 
#include <vector>
#include "../include/stockData.h"

void computeRollingAverages(std::vector<double> &data, std::vector<double>& avgs, int window_size=5);
void computeBollingerBands(std::vector<double> &data, std::vector<double> &upperBand, std::vector<double> &lowerBand, int window_size=5);
void computeEMA(std::vector<double> &data, std::vector<double>& ema, int period);
void computeMACD(std::vector<double> &data, std::vector<double> &macd, int shortPeriod, int longPeriod);
void computeStatistics(const std::vector<double> &data, double &mean, double &variance, double &stddev, double &skewness, double &kurtosis);
void detectVolumeSpikes(const std::vector<int> &volume, std::vector<int> &spikeIndices, double threshold=2.0);
void detectFlashCrashes(const std::vector<double> &closePrices, std::vector<int> &crashIndices, double dropThreshold = 5.0);


#endif