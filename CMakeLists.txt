cmake_minimum_required(VERSION 3.20)
project(StockAnalysis)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Enable AVX/AVX2 optimizations
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -O2 -Wall -Wextra")

# Include directories
include_directories(include)

# Source files
file(GLOB SOURCES "src/*.cpp")

# Create executable
add_executable(stock_analyzer ${SOURCES})
