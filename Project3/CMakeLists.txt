project(testProject)
cmake_minimum_required(VERSION 3.1)

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)
set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/lib)

include_directories(include)

add_library(test src/test.cpp)

add_executable(main src/main.cpp)

target_link_libraries(main test)
