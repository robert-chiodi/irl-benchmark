set(CURRENT_DIR "${CMAKE_CURRENT_SOURCE_DIR}")
target_sources(voftools_functions PRIVATE "${CURRENT_DIR}/voftools_wrapper.f90")
target_sources(voftools_functions PRIVATE "${CURRENT_DIR}/voftools_intersection.f90")
target_sources(voftools_functions PRIVATE "${CURRENT_DIR}/voftools_intersection.h")

