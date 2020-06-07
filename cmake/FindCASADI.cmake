# Find CasADi header
find_path(CASADI_INCLUDE_DIR
  casadi/casadi.hpp
  HINTS $ENV{CASADI_PREFIX}/include
)

# Find CasADi shared library
find_library(CASADI_LIBRARY NAMES casadi
    HINTS ${CASADI_INCLUDE_DIR}/../lib $ENV{CASADI_PREFIX}/lib
)

message("CASADI_INCLUDE_DIR: ${CASADI_INCLUDE_DIR}")
message("CASADI_LIBRARY: ${CASADI_LIBRARY}")

if (CASADI_INCLUDE_DIR AND CASADI_LIBRARY)
    # Create a target to represent the library
    add_library(casadi SHARED IMPORTED)

    # Attach the shared library to the target
    set_target_properties(casadi PROPERTIES IMPORTED_LOCATION ${CASADI_LIBRARY})

    # Attach the header to the target
    set_target_properties(casadi PROPERTIES INTERFACE_INCLUDE_DIRECTORIES ${CASADI_INCLUDE_DIR})
else()
    message(FATAL_ERROR "Could not find CasADi")
endif()
