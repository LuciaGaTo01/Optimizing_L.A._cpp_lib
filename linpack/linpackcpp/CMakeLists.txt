cmake_minimum_required(VERSION 3.0.0)
project(latest VERSION 0.1.0)

set(CMAKE_CXX_STANDARD 20)

include_directories(linpack
                    PRIVATE
                        <linear-algebra-installation-path>
)

add_executable(linpack main.cpp)