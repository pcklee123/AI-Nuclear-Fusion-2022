cmake_minimum_required(VERSION 3.0.0)

project(TS3)

# Set the C compiler
#set(CMAKE_C_COMPILER "C:\msys64\ucrt64\bin\gcc.exe")
#set(CMAKE_C_COMPILER "C:\msys64\ucrt64\bin\gcc.exe")

# Set the CXX compiler
#set(CMAKE_CXX_COMPILER "C:\msys64\ucrt64\bin\g++.exe")

#include(CTest)
#enable_testing()

# Find and link FFTW
find_package(FFTW REQUIRED)
include_directories(${FFTW_INCLUDE_DIRS})
link_directories(${FFTW_LIBRARY_DIRS})
add_definitions(${FFTW_DEFINITIONS})
set(LIBS ${LIBS} ${FFTW_LIBRARIES})

# Find and link GSL
find_package(GSL REQUIRED)
include_directories(${GSL_INCLUDE_DIRS})
link_directories(${GSL_LIBRARY_DIRS})
set(LIBS ${LIBS} ${GSL_LIBRARIES})

# Find and link VTK
find_package(VTK REQUIRED)
include(${VTK_USE_FILE})
set(LIBS ${LIBS} ${VTK_LIBRARIES})

# Find and link OpenCL
find_package(OpenCL REQUIRED)
include_directories(${OpenCL_INCLUDE_DIRS})
link_directories(${OpenCL_LIBRARY_DIRS})
set(LIBS ${LIBS} ${OpenCL_LIBRARIES})

# Enable OpenMP
find_package(OpenMP REQUIRED)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set(LIBS ${LIBS} ${OpenMP_CXX_LIBRARIES})


# Add all your source files
set(SOURCES 
    TS3.cpp 
    calc_trilin_constants.cpp
    calcEBV_FFT.cpp
    calcU.cpp
    changedt.cpp
    cl_code.cpp
    generate.cpp
    get_densityfields.cpp
    save.cpp
    sel_part_print.cpp
    utils.cpp
)

# Create an executable target and link the libraries
add_executable(my_executable ${SOURCES})
set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)

# Link the libraries
target_link_libraries(my_executable ${LIBS})