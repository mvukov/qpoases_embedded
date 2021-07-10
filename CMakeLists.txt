cmake_minimum_required(VERSION 2.8.3)

cmake_policy(SET CMP0042 NEW)

project(qpoases_embedded)

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

install(DIRECTORY include/${PROJECT_NAME}/
  DESTINATION include/${PROJECT_NAME}
  PATTERN include/qpoases_embedded/Utils.* EXCLUDE
)
install(TARGETS ${PROJECT_NAME}
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION bin
)
