SET(LIBRARY GNUPLOT_IOS)
SET(LIBRARYNAME Gnuplot-iostream)
FIND_PACKAGE(Boost COMPONENTS iostreams REQUIRED)
FIND_PATH(${LIBRARY}_INCLUDE_DIR gnuplot-iostream/gnuplot-iostream.h # insert the header file name here
    PATHS 
    /opt/local/include
    $ENV{HOME}/build/include
    )

set(${LIBRARY}_LIBRARIES ${Boost_LIBRARIES})
set(${LIBRARY}_INCLUDE_DIR ${LIBRARY}_INCLUDE_DIR ${Boost_INCLUDE_DIR})

if (${LIBRARY}_INCLUDE_DIR AND ${LIBRARY}_LIBRARIES)# AND ${LIBRARY}_LIBRARIES)
    set (${LIBRARY}_FOUND true)
endif(${LIBRARY}_INCLUDE_DIR AND ${LIBRARY}_LIBRARIES)# AND ${LIBRARY}_LIBRARIES)


if (${LIBRARY}_FOUND)
    if (NOT ${LIBRARY}_FIND_QUIETLY)
        message(STATUS "Found ${LIBRARYNAME}: ${${LIBRARY}_LIBRARIES}")
    endif(NOT ${LIBRARY}_FIND_QUIETLY)
else (${LIBRARY}_FOUND)
   if (${LIBRARY}_FIND_REQUIRED)
      message(FATAL_ERROR "Could not find ${LIBRARYNAME}")
    endif(${LIBRARY}_FIND_REQUIRED)
endif(${LIBRARY}_FOUND) 