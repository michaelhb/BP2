find_package(Boost COMPONENTS iostreams regex system REQUIRED) 

add_library(CasadiBounceSolver CasadiBounceSolver.cpp BouncePath.cpp)

target_include_directories(CasadiBounceSolver PUBLIC ../include)

# NOTE - the boost libs are dependencies of gnuplot_ios. I'd prefer 
# to add them to the gnuplot_ios target, but I don't know how to do
# that on a header-only imported library...
target_link_libraries(CasadiBounceSolver 
    PUBLIC casadi 
    PRIVATE gnuplot_ios Boost::iostreams Boost::regex Boost::system
)

