cmake_minimum_required(VERSION 3.0.0)
project(py_matrix VERSION 0.1.0 LANGUAGES C CXX)

include(CTest)
enable_testing()


find_package(CUDA)
if(CUDA_FOUND)
    add_compile_definitions(CUDA_SUPPORTED)
endif()
add_compile_definitions(TEST_TEST)


include_directories(inc)

find_package(Python3 REQUIRED COMPONENTS Development)

add_library(py_matrix SHARED src/MatrixMisc.cpp src/Matrix.cpp src/py_matrix.cpp)
target_link_libraries(py_matrix PRIVATE Python3::Python)
target_compile_options(py_matrix PRIVATE -w -Wfatal-errors -O3 -Wall -std=c++17 -fPIC)
target_include_directories(py_matrix PRIVATE ${Python3_INCLUDE_DIRS})
set_target_properties(py_matrix PROPERTIES PREFIX "" OUTPUT_NAME "py_matrix${CMAKE_SHARED_LIBRARY_SUFFIX}" SUFFIX "" LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR})

if(CUDA_FOUND)
    # target_include_directories(py_matrix PRIVATE ${CUDA_INCLUDE_DIRS})
    target_link_libraries(py_matrix PRIVATE ${CUDA_LIBRARIES})
    set_target_properties(py_matrix PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
endif()

find_program(STUBGEN_EXECUTABLE
    NAMES stubgen
    PATH_SUFFIXES bin
)

if(NOT STUBGEN_EXECUTABLE)
    message(WARNING "Stub generation tool 'stubgen' from 'mypy' package not found. Stub generation will be skipped.
Install 'mypy' package to generate the stub")
else()
    add_custom_target(generate_stub ALL
        COMMAND ${STUBGEN_EXECUTABLE} -m py_matrix -o ${CMAKE_SOURCE_DIR}
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        COMMENT "Generating stub"
        VERBATIM
    )

    add_dependencies(generate_stub py_matrix)
endif()

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
