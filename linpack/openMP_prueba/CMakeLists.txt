cmake_minimum_required(VERSION 3.0.0)
project(linalg LANGUAGES CXX VERSION 0.1.0)

set(CMAKE_CXX_STANDARD 20)

include_directories(linalg
                    PRIVATE
                        ~/Escritorio/TFG/cloned_repos/Optimizing_L.A._cpp_lib/linpack/openMP_prueba/LA_lib
)

find_package(OpenMP REQUIRED)
add_executable(linalg main.cpp)
target_link_libraries(linalg PRIVATE OpenMP::OpenMP_CXX)
