cmake_minimum_required(VERSION 3.21)

project(coolercpp
    VERSION 0.4.0
    DESCRIPTION "C++ library compatible with the cooler Python API and file format"
    LANGUAGES CXX)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

option(COOLERCPP_BUILD_TESTS "Build the doctest unit tests" ${PROJECT_IS_TOP_LEVEL})
option(COOLERCPP_BUILD_HARNESS "Build the Python equivalence harness driver"
       ${PROJECT_IS_TOP_LEVEL})

if(PROJECT_IS_TOP_LEVEL)
    set(CMAKE_CXX_STANDARD 20)
    set(CMAKE_CXX_STANDARD_REQUIRED ON)
    set(CMAKE_CXX_EXTENSIONS OFF)
    if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
        set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
    endif()
    # Installed binaries find libhdf5 through the path it was linked from, so
    # nothing needs LD_LIBRARY_PATH.
    set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
endif()

# Native dependencies: the HDF5 C library and zlib.
find_package(HDF5 REQUIRED COMPONENTS C)
# The threaded chunk writer and reader (fastio) filter chunks themselves.
find_package(ZLIB REQUIRED)
find_package(Threads REQUIRED)

add_library(coolercpp
    src/api.cpp
    src/create.cpp
    src/dtype.cpp
    src/fastio.cpp
    src/fileops.cpp
    src/h5.cpp
    src/json.cpp
    src/numpy_compat.cpp
    src/parallel.cpp
    src/rangequery.cpp
    src/rangequery_api.cpp
    src/region.cpp
    src/table.cpp
    src/util.cpp)
add_library(coolercpp::coolercpp ALIAS coolercpp)

target_compile_features(coolercpp PUBLIC cxx_std_20)
target_include_directories(coolercpp
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src)
# HDF5 is an implementation detail: no public header includes hdf5.h.
target_link_libraries(coolercpp PRIVATE HDF5::HDF5 ZLIB::ZLIB PUBLIC Threads::Threads)
set_target_properties(coolercpp PROPERTIES
    POSITION_INDEPENDENT_CODE ON
    VERSION ${PROJECT_VERSION}
    SOVERSION ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR})
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(coolercpp PRIVATE -Wall -Wextra -Wpedantic)
endif()

configure_file(include/coolercpp/version.hpp.in
    ${CMAKE_CURRENT_BINARY_DIR}/generated/coolercpp/version.hpp @ONLY)
target_include_directories(coolercpp PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/generated>)

# ---- install and package config ----
install(TARGETS coolercpp EXPORT coolercppTargets
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
install(DIRECTORY include/coolercpp DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    FILES_MATCHING PATTERN "*.hpp")
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/generated/coolercpp/version.hpp
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/coolercpp)
install(EXPORT coolercppTargets
    NAMESPACE coolercpp::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/coolercpp)
configure_package_config_file(cmake/coolercppConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/coolercppConfig.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/coolercpp)
write_basic_package_version_file(${CMAKE_CURRENT_BINARY_DIR}/coolercppConfigVersion.cmake
    COMPATIBILITY SameMinorVersion)
install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/coolercppConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/coolercppConfigVersion.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/coolercpp)

if(COOLERCPP_BUILD_HARNESS)
    add_subdirectory(harness)
endif()

if(COOLERCPP_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()
