cmake_minimum_required(VERSION 3.13)

project(string_ebmedding)

set(CMAKE_CXX_STANDARD 17)

add_executable(test_embed test_embedding.cpp command_group.cpp)
add_executable(gen_seq gen_seqs.cpp )
add_library(testlib STATIC command_group.cpp)
