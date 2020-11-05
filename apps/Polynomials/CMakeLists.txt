add_executable(Polynomials Polynomials.cpp)
target_link_libraries(Polynomials CasadiBounceSolver ${CASADI_LIBRARIES})
set_target_properties(Polynomials PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

add_executable(ActionDeltaPlot ActionDeltaPlot.cpp)
target_link_libraries(ActionDeltaPlot CasadiBounceSolver ${CASADI_LIBRARIES})
set_target_properties(ActionDeltaPlot PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)
configure_file(ActionDeltaPlot.py ${CMAKE_SOURCE_DIR}/bin COPYONLY)
configure_file(ct_potential.py ${CMAKE_SOURCE_DIR}/bin COPYONLY)

add_executable(StandardMultiField StandardMultiField.cpp)
target_link_libraries(StandardMultiField CasadiBounceSolver ${CASADI_LIBRARIES})
set_target_properties(StandardMultiField PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

set(GEN_PLOT ${CMAKE_SOURCE_DIR}/bin/ActionDeltaPlot)
set(GEN_PLOT_CT "python3 -W ignore ${CMAKE_SOURCE_DIR}/bin/ct_potential.py")
set(SHOW_PLOT ${CMAKE_SOURCE_DIR}/bin/ActionDeltaPlot.py)
set(OUTPUT_DIR ${CMAKE_BINARY_DIR})

configure_file(actiondelta_gridpoints.sh ${CMAKE_SOURCE_DIR}/bin)
configure_file(compare_to_ct.sh ${CMAKE_SOURCE_DIR}/bin)
