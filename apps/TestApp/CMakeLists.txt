add_executable(TestApp TestApp.cpp)
target_link_libraries(TestApp CasadiBounceSolver ${CASADI_LIBRARIES})
target_include_directories(TestApp PUBLIC ${CASADI_INCLUDE_DIR})
set_target_properties(TestApp PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)
