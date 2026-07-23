cmake_minimum_required(VERSION 3.10)
project(SynchrotronTracker VERSION 1.0 LANGUAGES CXX)

# --- Стандарт C++ ---
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# --- OpenMP (требуется для распараллеливания по частицам) ---
find_package(OpenMP REQUIRED)

# --- Бинарники ---
# Все реализации в *.h, поэтому здесь только .cpp с main().
add_executable(SynchrotronTracker   src/main.cpp)
add_executable(SynchrotronOptimizer src/optimizer.cpp)
add_executable(MatrixTracker        src/matrix_tracker_main.cpp)   # матричное ядро (matrix_core.h)

foreach(target SynchrotronTracker SynchrotronOptimizer MatrixTracker)
    target_include_directories(${target} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
    target_link_libraries     (${target} PUBLIC OpenMP::OpenMP_CXX)

    if(MSVC)
        target_compile_options(${target} PRIVATE /W4 /openmp /O2)
    else()
        target_compile_options(${target} PRIVATE -Wall -Wextra -O3 -fopenmp)
    endif()
endforeach()

message(STATUS "OpenMP found: ${OpenMP_CXX_FOUND}")
message(STATUS "Targets: SynchrotronTracker, SynchrotronOptimizer")

# workspace/ остаётся в корне проекта; .exe находит его сам через
# findWorkspace() (см. src/main.cpp / src/optimizer.cpp).
