cmake_minimum_required(VERSION 3.12)
project(lsopt)

set(CMAKE_CXX_STANDARD 17)

add_executable(lsopt main.cpp hyperdual.h optimizer.h)