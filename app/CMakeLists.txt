# Author: J. Xia
# Date: Aug 8, 2024

############
### Apps ###
############

# Find the ROOT package
find_package(ROOT REQUIRED)
# Include the ROOT headers
include_directories(${ROOT_INCLUDE_DIRS})

# Add an executable target for the tool
add_executable(check_uniform_voxel check_uniform_voxel.cc)

set( WCSIM_LIBS_TO_EXPORT WCSimCore WCSimRoot ${ROOT_LIBRARIES})

target_link_libraries(check_uniform_voxel PRIVATE ${WCSIM_LIBS_TO_EXPORT})
target_include_directories(check_uniform_voxel PRIVATE ${CMAKE_SOURCE_DIR}/src)
