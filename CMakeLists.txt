cmake_minimum_required(VERSION 3.0)
project(kiro2021)

find_package(nlohmann_json 3.2.0 REQUIRED)

add_executable(kiro main.cpp instance.h instance.cpp)
