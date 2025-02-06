#ifndef FILEIO_H
#define FILEIO_H

#pragma once
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include "../include/stockData.h"

int stringToDate(const std::string& dateStr);
void readStockMetaData(std::vector<std::string>& stocks, std::string fileName);
void readStockFile(StockData& stockData);

#endif
