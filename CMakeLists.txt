# Almost all CMake files should start with this
# You should always specify a range with the newest
# and oldest tested versions of CMake. This will ensure
# you pick up the best policies.
cmake_minimum_required(VERSION 3.11)

# This is your project statement. You should always list languages;
# Listing the version is nice here since it sets lots of useful variables
project(TIMING_COMP VERSION 0.1
	LANGUAGES CXX C Fortran)

set(CMAKE_CXX_STANDARD 14)

### Require out-of-source builds
file(TO_CMAKE_PATH "${PROJECT_BINARY_DIR}/CMakeLists.txt" LOC_PATH)
if(EXISTS "${LOC_PATH}")
    message(FATAL_ERROR "You cannot build in a source directory (or any directory with a CMakeLists.txt file). Please make a build subdirectory.")
endif()

set(TIMING_SOURCE "${PROJECT_SOURCE_DIR}/src")


# Add libraries
add_library(irl_functions STATIC)
set_target_properties(irl_functions
    PROPERTIES
    COMPILE_FLAGS "${CMAKE_CXX_FLAGS}"
    CXX_STANDARD 14)
target_include_directories(irl_functions PRIVATE "./")
target_link_libraries(irl_functions PUBLIC "${IRL_INSTALL_LOCATION}/lib/libirl.a")
target_link_libraries(irl_functions PUBLIC "${IRL_INSTALL_LOCATION}/absl/lib/libabsl_all.a")
target_link_libraries(irl_functions PUBLIC "-lpthread")
if(APPLE)
  target_link_libraries(irl_functions PUBLIC "-framework CoreFoundation")
endif()
target_include_directories(irl_functions PUBLIC "${IRL_ROOT_LOCATION}")
target_include_directories(irl_functions SYSTEM PUBLIC "${IRL_ROOT_LOCATION}/external/abseil-cpp/")

add_library(irl_gvm_functions STATIC)
set_target_properties(irl_gvm_functions
    PROPERTIES
    COMPILE_FLAGS "${CMAKE_CXX_FLAGS}"
    CXX_STANDARD 14)
target_include_directories(irl_gvm_functions PRIVATE "./")
target_link_libraries(irl_gvm_functions PUBLIC "${IRL_INSTALL_LOCATION}/lib/libirl.a")
target_link_libraries(irl_gvm_functions PUBLIC "${IRL_INSTALL_LOCATION}/absl/lib/libabsl_all.a")
target_link_libraries(irl_gvm_functions PUBLIC "-lpthread")
if(APPLE)
  target_link_libraries(irl_gvm_functions PUBLIC "-framework CoreFoundation")
endif()
target_include_directories(irl_gvm_functions PUBLIC "${IRL_ROOT_LOCATION}")
target_include_directories(irl_gvm_functions SYSTEM PUBLIC "${IRL_ROOT_LOCATION}/external/abseil-cpp/")

add_library(r3d_functions STATIC)
set_target_properties(r3d_functions
    PROPERTIES
    COMPILE_FLAGS "${CMAKE_C_FLAGS}"
    LANGUAGE "C")
target_include_directories(r3d_functions PRIVATE "./")    
target_include_directories(r3d_functions PUBLIC "${R3D_LOCATION}")
target_link_libraries(r3d_functions PUBLIC "${R3D_LOCATION}/libr3d.a")



add_library(voftools_functions STATIC)
set_target_properties(voftools_functions
    PROPERTIES
    COMPILE_FLAGS "${CMAKE_FORTRAN_FLAGS}"
    LANGUAGE "Fortran")
target_include_directories(voftools_functions PRIVATE "./")    
target_include_directories(voftools_functions PUBLIC "${VOFTOOLS_LOCATION}")
target_link_libraries(voftools_functions PUBLIC "${VOFTOOLS_LOCATION}/libvoftools.a")


# Add timing executable
add_executable(timing_comp)
set_target_properties(timing_comp
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}"
    COMPILE_FLAGS "${CMAKE_CXX_FLAGS}"
    INTERPROCEDURAL_OPTIMIZATION True
    CXX_STANDARD 14)


target_include_directories(timing_comp PRIVATE "./")
target_link_libraries(timing_comp PRIVATE irl_functions)
target_link_libraries(timing_comp PRIVATE irl_gvm_functions)
target_link_libraries(timing_comp PRIVATE r3d_functions)
target_link_libraries(timing_comp PRIVATE voftools_functions)

add_subdirectory("${TIMING_SOURCE}")

