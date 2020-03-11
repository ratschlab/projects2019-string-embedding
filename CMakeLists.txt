cmake_minimum_required(VERSION 3.13)

project(string_ebmedding)

set(CMAKE_CXX_STANDARD 17)

add_executable(test_embed test_embedding.cpp)
add_executable(gen_seq test_gen_seq.cpp)
add_library(testlib command_group.cpp)
