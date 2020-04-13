cmake_minimum_required(VERSION 2.8.3)

cmake_policy(SET CMP0042 NEW)

project(qpoases_embedded)

find_package(catkin QUIET)

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

if (CATKIN_ENABLE_TESTING)
  add_library(utils src/Utils.cpp)
  target_compile_definitions(utils PRIVATE QPOASES_DEBUG)

  catkin_add_gtest(e2e_tests
    test/test_example1.cpp
    test/test_example1b.cpp
    test/test_main.cpp)
  target_compile_definitions(e2e_tests PRIVATE QPOASES_DEBUG)
  target_link_libraries(e2e_tests ${PROJECT_NAME} utils)
endif()
