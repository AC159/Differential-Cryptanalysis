cmake_minimum_required(VERSION 3.23)
project(Soen321DifferentialCryptanalysis)

set(CMAKE_CXX_STANDARD 20)

add_executable(Soen321DifferentialCryptanalysis main.cpp wes-key-7-main.cpp)

add_executable(wes wes/wes.c)
