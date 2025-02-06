#include "../include/fileIO.h"
#include "../include/parallel.h"
#include "../include/stockData.h"
#include <unordered_map>
#include <chrono>
#include <iostream>

int main() {
    auto start = std::chrono::high_resolution_clock::now();

    std::vector<std::string> stocks;
    readStockMetaData(stocks, "../stocks_metadata.csv");

    int totalFiles = stocks.size();
    std::vector<StockData> stockDataList(totalFiles);
    std::unordered_map<std::string, StockData*> stockDataMap;

    for (int i = 0; i < totalFiles; i++) {
        stockDataList[i].name = "../stocks/" + stocks[i] + ".csv";
        readStockFile(stockDataList[i]);
        stockDataMap[stocks[i]] = &stockDataList[i];
    }
    std::cout << stockDataList.size() << " " << totalFiles << std::endl;

    std::vector<StockMetrics> results = parallelProcessing(stockDataList);
    std::unordered_map<std::string, StockMetrics> stockMetricsMap;
    for (int i = 0; i < totalFiles; i++) {
        stockMetricsMap[stocks[i]] = results[i];
    }
    std::cout << "Metrics size: " << results.size() << std::endl;

    // std::string s;
    // std::cout << "Enter Stock Symbol: ";
    // getline(std::cin, s);
    // if (stockMetricsMap[s]) {

    // }


    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Time taken: " << elapsed.count() << " seconds" << std::endl;

    return 0;
}