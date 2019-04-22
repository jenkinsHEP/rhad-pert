cmake_minimum_required(VERSION 3.13.4)
enable_language(Fortran)
project(rhad-pert)
set(CMAKE_CXX_STANDARD 11)

file(GLOB SOURCES "src/core/*.cpp" "src/core/rhad-1.01/*.f")

#main program
add_executable(rhad-pert ${SOURCES} ../src/main.cpp)