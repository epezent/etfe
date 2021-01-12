cmake_minimum_required(VERSION 3.0.0)
project(fft VERSION 0.1.0)

include(FetchContent) 
FetchContent_Declare(mahi-gui GIT_REPOSITORY https://github.com/mahilab/mahi-gui.git) 
FetchContent_MakeAvailable(mahi-gui)

FetchContent_Declare(iir GIT_REPOSITORY https://github.com/epezent/iir1) 
FetchContent_MakeAvailable(iir)

add_library(kissfft "kissfft/kiss_fft.h" "kissfft/kiss_fftr.h" "kissfft/kiss_fft.c" "kissfft/kiss_fftr.c")
target_include_directories(kissfft PUBLIC "kissfft")

add_executable(filter_toy "filter_toy.cpp" "ETFE.hpp")
target_link_libraries(filter_toy mahi::gui kissfft iir::iir_static)

add_executable(etfe "etfe_cl.cpp" "ETFE.hpp")
target_link_libraries(etfe mahi::util kissfft)