cmake_minimum_required(VERSION 2.8)

project(4DOFReg)

set(CMAKE_BUILD_TYPE Release)


include_directories(include)

file(GLOB_RECURSE SOURCES "src/*.cpp")

add_executable(${PROJECT_NAME} ${SOURCES})

#set your eigen3 include dir
include_directories(/usr/include/eigen3)

# Find PCL
find_package( PCL 1.8 REQUIRED )

if( PCL_FOUND )
  # [C/C++]>[General]>[Additional Include Directories]
  include_directories( ${PCL_INCLUDE_DIRS} )

  # [C/C++]>[Preprocessor]>[Preprocessor Definitions]
  add_definitions( ${PCL_DEFINITIONS} )
  
  # For Use Not PreCompiled Features 
  #add_definitions( -DPCL_NO_PRECOMPILE )

  # [Linker]>[General]>[Additional Library Directories]
  link_directories( ${PCL_LIBRARY_DIRS} )

  # [Linker]>[Input]>[Additional Dependencies]
  target_link_libraries(${PROJECT_NAME} ${PCL_LIBRARIES} )
endif()


# Find cblas
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")
find_package(CBLAS REQUIRED)
include_directories(${CBLAS_INCLUDE_DIR})
target_link_libraries(${PROJECT_NAME} ${CBLAS_LIBRARIES})

## Find super4pcs(used after installing super4pcs), install super4pcs library and uncomment if you want to use it.
##set this to the cmake sub-folder of your S4PCS installation dir. And change the filename of "Super4PCSConfig.cmake" to "FindSuper4PCS.cmake"
#list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/S4PCS/Super4PCS/build/install/lib/cmake")
#find_package(Super4PCS REQUIRED)
#include_directories(${Super4PCS_INCLUDE_DIR})
#link_directories(${Super4PCS_LIB_DIR})
#target_link_libraries(${PROJECT_NAME} ${Super4PCS_LIBRARIES})
