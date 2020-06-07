add_executable(TestApp TestApp.cpp)
target_link_libraries(TestApp CasadiBounceSolver)
set_target_properties(TestApp PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)
