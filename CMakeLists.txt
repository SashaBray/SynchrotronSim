cmake_minimum_required(VERSION 3.10)
project(SynchrotronTracker VERSION 1.0 LANGUAGES CXX)

# 1. Настройка стандарта C++ (используем 17 для современных функций)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# 2. Поиск и настройка OpenMP для многопоточности
find_package(OpenMP REQUIRED)

# 3. Указываем исходные файлы проекта
# Добавьте сюда .cpp файлы, когда вынесете реализацию из заголовков
set(SOURCES 
    src/main.cpp
    # src/FieldMap.cpp (если создадите)
    # src/AnsysConverter.cpp
)

# 4. Создание исполняемого файла
add_executable(${PROJECT_NAME} ${SOURCES})

# 5. Подключаем OpenMP к нашему проекту
target_link_libraries(${PROJECT_NAME} PUBLIC OpenMP::OpenMP_CXX)

# 6. Дополнительные настройки для разных ОС
if(MSVC)
    # Настройки для Visual Studio
    target_compile_options(${PROJECT_NAME} PRIVATE /W4 /openmp)
else()
    # Настройки для GCC/Clang
    target_compile_options(${PROJECT_NAME} PRIVATE -Wall -Wextra -O3 -fopenmp)
endif()

# 7. Указываем, где лежат заголовочные файлы (.h)
target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)

# Сообщение для пользователя при конфигурации
message(STATUS "OpenMP found: ${OpenMP_CXX_FOUND}")


# Копирование папки данных в бинарную директорию при сборке
add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
    "${CMAKE_CURRENT_SOURCE_DIR}/workspace"
    "$<TARGET_FILE_DIR:${PROJECT_NAME}>/workspace"
)

