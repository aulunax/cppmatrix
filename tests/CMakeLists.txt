enable_testing()

add_executable(CoreTest core_features.cc)
add_executable(OperatorsTest basic_operators.cc)
add_executable(AdvOperationsTest advanced_operations.cc)

target_link_libraries(CoreTest GTest::gtest_main py_matrix)
target_link_libraries(OperatorsTest GTest::gtest_main py_matrix)
target_link_libraries(AdvOperationsTest GTest::gtest_main py_matrix)

include_directories(${CMAKE_SOURCE_DIR}/inc)

option(gtest_build_tests "Build all of gtest's own tests." OFF)
option(gtest_build_samples "Build gtest's sample programs." OFF)

include(GoogleTest)
gtest_discover_tests(CoreTest)
gtest_discover_tests(OperatorsTest)
gtest_discover_tests(AdvOperationsTest)