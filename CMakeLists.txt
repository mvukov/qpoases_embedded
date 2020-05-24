cmake_minimum_required(VERSION 2.8.3)

cmake_policy(SET CMP0042 NEW)

project(qpoases_embedded)

find_package(catkin QUIET)
find_package(benchmark QUIET)

if (catkin_FOUND)
  catkin_package(
    INCLUDE_DIRS include
    LIBRARIES ${PROJECT_NAME}
  )
endif()

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

option(QPOASES_DEBUG "Enables debug messages." OFF)

# Add -Wconversion -Wsign-conversion -Wshadow.
add_compile_options(-std=c++17 -Wall -Wextra -Wpedantic -Wfloat-equal -Werror)

include_directories(include)

add_library(${PROJECT_NAME} SHARED
  src/Bounds.cpp
  src/Constraints.cpp
  src/Indexlist.cpp
  src/MessageHandling.cpp
  src/QProblem.cpp
  src/QProblemB.cpp
  src/SubjectTo.cpp
)
if (QPOASES_DEBUG)
  target_compile_definitions(${PROJECT_NAME} PRIVATE QPOASES_DEBUG)
endif()

if (catkin_FOUND)
  install(DIRECTORY include/${PROJECT_NAME}/
    DESTINATION ${CATKIN_PACKAGE_INCLUDE_DESTINATION}
  )
  install(TARGETS ${PROJECT_NAME}
    ARCHIVE DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
    LIBRARY DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
    RUNTIME DESTINATION ${CATKIN_PACKAGE_BIN_DESTINATION}
  )
endif()

set(GENERATED_INCLUDE_DIR ${PROJECT_BINARY_DIR}/include)
set(GENERATED_TEST_DATA_DIR ${GENERATED_INCLUDE_DIR}/${PROJECT_NAME})

if (benchmark_FOUND)
  set(CHAIN_MASS_NMPC_BENCHMARK_DATA_FILE
    ${GENERATED_TEST_DATA_DIR}/chain_mass_nmpc_benchmark_data.h)
  set(HANGING_CHAIN_BENCHMARK_DATA_FILE
    ${GENERATED_TEST_DATA_DIR}/hanging_chain_benchmark_data.h)

  add_custom_command(
    OUTPUT ${HANGING_CHAIN_BENCHMARK_DATA_FILE}
    COMMAND python3 ${CMAKE_CURRENT_SOURCE_DIR}/test/hanging_chain.py
    ARGS ${HANGING_CHAIN_BENCHMARK_DATA_FILE}
    DEPENDS
      ${CMAKE_CURRENT_SOURCE_DIR}/test/common.py
      ${CMAKE_CURRENT_SOURCE_DIR}/test/hanging_chain.py
    COMMENT "Exporting hanging chain QP benchmark data vectors."
    VERBATIM
  )

  add_custom_command(
    OUTPUT ${CHAIN_MASS_NMPC_BENCHMARK_DATA_FILE}
    COMMAND python3 ${CMAKE_CURRENT_SOURCE_DIR}/test/chain_mass_nmpc.py
    ARGS ${CHAIN_MASS_NMPC_BENCHMARK_DATA_FILE}
    DEPENDS
      ${CMAKE_CURRENT_SOURCE_DIR}/test/common.py
      ${CMAKE_CURRENT_SOURCE_DIR}/test/chain_mass_nmpc.py
    COMMENT "Exporting hanging chain QP benchmark data vectors."
    VERBATIM
  )

  add_executable(benchmark_hanging_chain
    ${HANGING_CHAIN_BENCHMARK_DATA_FILE}
    test/benchmark_hanging_chain.cpp)
  target_include_directories(benchmark_hanging_chain
    PRIVATE benchmark::benchmark ${GENERATED_INCLUDE_DIR})
  target_link_libraries(benchmark_hanging_chain
    ${PROJECT_NAME} benchmark::benchmark)

  add_executable(benchmark_chain_mass_nmpc
    ${CHAIN_MASS_NMPC_BENCHMARK_DATA_FILE}
    test/benchmark_chain_mass_nmpc.cpp)
  target_include_directories(benchmark_chain_mass_nmpc
    PRIVATE benchmark::benchmark ${GENERATED_INCLUDE_DIR})
  target_link_libraries(benchmark_chain_mass_nmpc
    ${PROJECT_NAME} benchmark::benchmark)
endif()

if (CATKIN_ENABLE_TESTING)
  add_library(utils src/Utils.cpp)
  target_compile_definitions(utils PRIVATE QPOASES_DEBUG)

  set(HANGING_CHAIN_DATA_FILE
    ${GENERATED_TEST_DATA_DIR}/hanging_chain_test_data.h)
  set(CHAIN_MASS_NMPC_DATA_FILE
    ${GENERATED_TEST_DATA_DIR}/chain_mass_nmpc_test_data.h)

  add_custom_command(
    OUTPUT ${HANGING_CHAIN_DATA_FILE}
    COMMAND python3 ${CMAKE_CURRENT_SOURCE_DIR}/test/hanging_chain.py
    ARGS ${HANGING_CHAIN_DATA_FILE}
    DEPENDS
      ${CMAKE_CURRENT_SOURCE_DIR}/test/common.py
      ${CMAKE_CURRENT_SOURCE_DIR}/test/hanging_chain.py
    COMMENT "Exporting hanging chain QP test data vectors."
    VERBATIM
  )

  add_custom_command(
    OUTPUT ${CHAIN_MASS_NMPC_DATA_FILE}
    COMMAND python3 ${CMAKE_CURRENT_SOURCE_DIR}/test/chain_mass_nmpc.py
    ARGS ${CHAIN_MASS_NMPC_DATA_FILE}
    DEPENDS
      ${CMAKE_CURRENT_SOURCE_DIR}/test/common.py
      ${CMAKE_CURRENT_SOURCE_DIR}/test/chain_mass_nmpc.py
    COMMENT "Exporting chain-mass NMPC QP test data vectors."
    VERBATIM
  )

  catkin_add_gtest(e2e_tests
    ${CHAIN_MASS_NMPC_DATA_FILE}
    ${HANGING_CHAIN_DATA_FILE}
    test/test_chain_mass_nmpc.cpp
    test/test_example1.cpp
    test/test_example1b.cpp
    test/test_hanging_chain.cpp
    test/test_main.cpp)
  target_compile_definitions(e2e_tests PRIVATE QPOASES_DEBUG)
  target_include_directories(e2e_tests PRIVATE ${GENERATED_INCLUDE_DIR})
  target_link_libraries(e2e_tests ${PROJECT_NAME} utils)
endif()
