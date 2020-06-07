# Find the gnuplot-ios header
find_path(GNUPLOT_IOS_INCLUDE_DIR gnuplot-iostream.h # insert the header file name here
    PATHS 
    /opt/local/include
    $ENV{HOME}/build/include
    ) 

if(GNUPLOT_IOS_INCLUDE_DIR)
    # Create a target to represent the header
    add_library(gnuplot_ios INTERFACE IMPORTED)

    # Add the header to the target
    set_target_properties(gnuplot_ios PROPERTIES INTERFACE_INCLUDE_DIRECTORIES ${GNUPLOT_IOS_INCLUDE_DIR})
    message("GNUPLOT_IOS_INCLUDE_DIR: ${GNUPLOT_IOS_INCLUDE_DIR}")
else()
    message(FATAL_ERROR "Could not find gnuplot_ios")
endif()