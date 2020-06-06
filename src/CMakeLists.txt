add_library(CasadiBounceSolver CasadiBounceSolver.cpp)

target_include_directories(CasadiBounceSolver PUBLIC ../include ${GNUPLOT_IOS_INCLUDE_DIR})
target_link_libraries(CasadiBounceSolver PRIVATE ${GNUPLOT_IOS_LIBRARIES})

