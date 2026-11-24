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

# --- Python-расширение synchrotronsim._engine (nanobind), шаг 0 docs/LIBRARY_ARCHITECTURE.md ---
# Отдельная папка сборки, чтобы не трогать кэш build/:
#   cmake -S . -B build-py -G "MinGW Makefiles" -DCMAKE_BUILD_TYPE=Release -DSS_BUILD_PYTHON=ON
#         -DPython_EXECUTABLE=<python.exe>
#   cmake --build build-py --target _engine
# Модуль кладётся прямо в пакет synchrotronsim/ (импорт без установки).
option(SS_BUILD_PYTHON "Собрать Python-расширение synchrotronsim._engine (nanobind)" OFF)
if(SS_BUILD_PYTHON)
    find_package(Python 3.10 COMPONENTS Interpreter Development.Module REQUIRED)
    execute_process(COMMAND "${Python_EXECUTABLE}" -m nanobind --cmake_dir
                    OUTPUT_STRIP_TRAILING_WHITESPACE OUTPUT_VARIABLE SS_NANOBIND_CMAKE)
    find_package(nanobind CONFIG REQUIRED PATHS "${SS_NANOBIND_CMAKE}" NO_DEFAULT_PATH)
    nanobind_add_module(_engine NB_STATIC src/py_engine.cpp)
    target_include_directories(_engine PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
    if(MINGW)
        # рантайм MinGW (libgcc, libstdc++, libgomp, winpthread) — статически: Python ≥ 3.8 не ищет
        # зависимые DLL по PATH, а модуль должен грузиться без C:\mingw64\bin в окружении.
        # OpenMP::OpenMP_CXX тянет импорт-библиотеку libgomp.dll.a, поэтому OpenMP здесь — флагом
        # драйвера: с -static он берёт libgomp.a
        target_compile_options(_engine PRIVATE -O3 -fopenmp)
        target_link_options(_engine PRIVATE -fopenmp -static -static-libgcc -static-libstdc++)
    elseif(MSVC)
        target_link_libraries(_engine PRIVATE OpenMP::OpenMP_CXX)
        target_compile_options(_engine PRIVATE /openmp /O2)
    else()
        target_link_libraries(_engine PRIVATE OpenMP::OpenMP_CXX)
        target_compile_options(_engine PRIVATE -O3 -fopenmp)
    endif()
    set_target_properties(_engine PROPERTIES
        LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/synchrotronsim
        LIBRARY_OUTPUT_DIRECTORY_RELEASE ${CMAKE_CURRENT_SOURCE_DIR}/synchrotronsim)
    message(STATUS "Python-расширение: ${Python_EXECUTABLE} (nanobind ${SS_NANOBIND_CMAKE})")
endif()

message(STATUS "OpenMP found: ${OpenMP_CXX_FOUND}")
message(STATUS "Targets: SynchrotronTracker, SynchrotronOptimizer")

# workspace/ остаётся в корне проекта; .exe находит его сам через
# findWorkspace() (см. src/main.cpp / src/optimizer.cpp).
