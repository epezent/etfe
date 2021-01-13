cmake_minimum_required(VERSION 3.0.0)
project(ETFE VERSION 0.1.0)

add_library(etfe "include/ETFE.hpp" "kissfft/kiss_fft.h" "kissfft/kiss_fftr.h" "kissfft/kiss_fft.c" "kissfft/kiss_fftr.c")
target_include_directories(etfe PUBLIC "include" "kissfft")

if (CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)

    include(FetchContent) 
    FetchContent_Declare(mahi-gui GIT_REPOSITORY https://github.com/mahilab/mahi-gui.git) 
    FetchContent_MakeAvailable(mahi-gui)

    FetchContent_Declare(iir GIT_REPOSITORY https://github.com/berndporr/iir1) 
    FetchContent_MakeAvailable(iir)

    add_executable(filter_toy "examples/filter_toy.cpp" "examples/filter_toy.rc")
    target_link_libraries(filter_toy mahi::gui etfe iir::iir_static)

    add_executable(etfe_cl "examples/etfe_cl.cpp")
    target_link_libraries(etfe_cl mahi::util etfe)
    set_target_properties(etfe_cl PROPERTIES OUTPUT_NAME "etfe")

endif()
