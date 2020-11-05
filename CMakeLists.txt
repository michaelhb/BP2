cmake_minimum_required(VERSION 3.9) # modern cmake pls

project(BubbleProfiler2
    VERSION 0.1
    DESCRIPTION "The best bounce solver in the universe"
    LANGUAGES CXX
)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake)

set(CMAKE_CXX_STANDARD 17)

find_package(CASADI REQUIRED)
find_package(gnuplot-iostream REQUIRED)

# Should only contain library code
add_subdirectory(src)

# Executable code lives here
add_subdirectory(apps)

