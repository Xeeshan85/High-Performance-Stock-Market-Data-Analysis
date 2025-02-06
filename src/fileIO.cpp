#include "../include/fileIO.h"
#include <cmath>

int stringToDate(const std::string& dateStr) {
    return stoi(dateStr.substr(0, 4)) * 10000 +
           stoi(dateStr.substr(5, 2)) * 100 +
           stoi(dateStr.substr(8, 2));
}

void readStockMetaData(std::vector<std::string>& stocks, std::string fileName) {
    std::ifstream fileIn;
    fileIn.open(fileName);
    if (!fileIn) {
        std::cerr << "Error opening file." << std::endl;
    }

    std::string line;
    getline(fileIn, line); // for header line
    while (getline(fileIn, line)) {
        std::stringstream ss(line);
        std::string cell, s;

        getline(ss, cell, ',');
        getline(ss, cell, ',');
        s = cell;
        getline(ss, cell, ',');
        if (cell == "N")
            stocks.push_back(s);
    }
    fileIn.close();

    // return stocks;
}

void readStockFile(StockData& stockData) {
    std::ifstream fileIn(stockData.name);
    if (!fileIn) {
        std::cerr << "Error opening file: " << stockData.name << std::endl;
        return;
    }

    std::string line;
    getline(fileIn, line); // Skip header

    while (getline(fileIn, line)) {
        std::stringstream ss(line);
        std::string cell;
        try {
            // Handling date
            getline(ss, cell, ',');
            try {
                stockData.date.push_back(stringToDate(cell));
            } catch (const std::exception& e) {
                std::cerr << "Error parsing date in line: " << line << " | " << e.what() << std::endl;
                continue;
            }

            // convert std::string to double safely
            auto safeStod = [](const std::string& str) -> double {
                try {
                    return str.empty() ? NAN : stod(str);
                } catch (...) {
                    return NAN; // Assign NAN on failure
                }
            };

            double opn, hig, lo, clos, vol;
            if (!getline(ss, cell, ',') || std::isnan(opn = safeStod(cell))) continue;
            if (!getline(ss, cell, ',') || std::isnan(hig = safeStod(cell))) continue;
            if (!getline(ss, cell, ',') || std::isnan(lo = safeStod(cell))) continue;
            if (!getline(ss, cell, ',') || std::isnan(clos = safeStod(cell))) continue;
            if (!getline(ss, cell, ',') || std::isnan(vol = safeStod(cell))) continue;

            stockData.open.push_back(opn);
            stockData.high.push_back(hig);
            stockData.low.push_back(lo);
            stockData.close.push_back(clos);
            stockData.volume.push_back(vol);

        } catch (const std::exception& e) {
            std::cerr << "Unexpected error while parsing line: " << line << " | " << e.what() << std::endl;
        }
    }

    fileIn.close();
}
