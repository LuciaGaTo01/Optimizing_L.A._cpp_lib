cmake_minimum_required(VERSION 3.0.0)
project(linalg LANGUAGES CXX VERSION 0.1.0)

set(CMAKE_CXX_STANDARD 20)

include_directories(linalg
                    PRIVATE
                        ~/Escritorio/TFG/cloned_repos/Optimizing_L.A._cpp_lib/linpack/LA_library/LA_lib
)

add_executable(linalg main.cpp)
