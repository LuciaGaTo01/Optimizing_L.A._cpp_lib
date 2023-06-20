cmake_minimum_required(VERSION 3.0.0)
project(latest VERSION 0.1.0)

set(CMAKE_CXX_STANDARD 20)

include(FetchContent)

FetchContent_Declare(linalg
    GIT_REPOSITORY "https://github.com/BobSteagall/wg21"
    GIT_SHALLOW ON
)

FetchContent_MakeAvailable(linalg)

include_directories(
  ${FETCHCONTENT_BASE_DIR}/linalg-src/include
  ${FETCHCONTENT_BASE_DIR}/mdspan-src/include
)

add_executable(linalg main.c)