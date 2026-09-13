cmake_minimum_required(VERSION 3.21)
project(coolercpp_consumer LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_BUILD_RPATH_USE_LINK_PATH TRUE)

if(CONSUME STREQUAL "find_package")
    find_package(coolercpp 0.2 REQUIRED)
elseif(CONSUME STREQUAL "FetchContent")
    include(FetchContent)
    FetchContent_Declare(coolercpp SOURCE_DIR "${COOLERCPP_SOURCE}")
    FetchContent_MakeAvailable(coolercpp)
else()
    message(FATAL_ERROR "CONSUME must be find_package or FetchContent")
endif()

add_executable(consumer main.cpp)
target_link_libraries(consumer PRIVATE coolercpp::coolercpp)
