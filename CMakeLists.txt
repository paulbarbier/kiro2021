cmake_minimum_required(VERSION 3.14)
project(kiro2021)

find_package(nlohmann_json 3.2.0 REQUIRED)
find_package(ortools CONFIG REQUIRED)

add_executable(kiro main.cpp instance.h instance.cpp)
target_link_libraries(kiro ortools::ortools)