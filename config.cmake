# Defaultly perform IPO
set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE CACHE STRING "PERFORM IPO")

# Change PATH_TO_IRL_ROOT to your own path to IRL (where the main CMakeLists.txt is)
set(IRL_ROOT_LOCATION "PATH_TO_IRL_ROOT" CACHE PATH "Path to IRL root" FORCE)
set(IRL_INSTALL_LOCATION "${IRL_ROOT_LOCATION}/install/Release" CACHE PATH "Path to installed IRL files" FORCE)

# Change PATH_TO_R3D_ROOT to your own path to R3D
set(R3D_LOCATION "PATH_TO_R3D_ROOT" CACHE PATH "Path to R3D root" FORCE)

# Change PATH_TO_VOFTOOLS_ROOT to your own path to VOFTools
set(VOFTOOLS_LOCATION "PATH_TO_VOFTOOLS_ROOT" CACHE PATH "Path to VOFTOOLS root" FORCE)

# Mark as Release build
set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build Type")

# Set compilers you will use if different from provided GNU default
set(CMAKE_CXX_COMPILER "g++" CACHE STRING "C++ Compiler")
set(CMAKE_C_COMPILER "gcc" CACHE STRING "C Compiler")
set(CMAKE_Fortran_COMPILER "gfortran" CACHE STRING "Fortran Compiler")

# Compiler flags for C++/C/Fortran. -fopenmp needed for timing
set(CMAKE_CXX_FLAGS "-O3 -DNDEBUG -DNDEBUG_PERF -march=native -fopenmp"
    CACHE STRING "C++ compile flags")

set(CMAKE_C_FLAGS "-O3 -DNDEBUG -march=native -fopenmp"
    CACHE STRING "C compile flags")

set(CMAKE_Fortran_FLAGS "-O3 -ffree-line-length-none -march=native -fopenmp"
    CACHE STRING "Fortran compile flags")

