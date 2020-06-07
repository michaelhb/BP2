add_executable(Polynomials Polynomials.cpp)
target_link_libraries(Polynomials CasadiBounceSolver ${CASADI_LIBRARIES})
set_target_properties(Polynomials PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)
