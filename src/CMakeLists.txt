

set(COMMON_INCLUDES ${PROJECT_SOURCE_DIR}/include)
set(MAIN_SOURCE "${PROJECT_SOURCE_DIR}/src/Main.cpp")
include_directories(${COMMON_INCLUDES})


file(GLOB SRCS "${PROJECT_SOURCE_DIR}/src/*.cpp")
list(REMOVE_ITEM SRCS "${MAIN_SOURCE}")

##################################
# add boost library
##################################
# set(Boost_USE_STATIC_LIBS ON) 
FIND_PACKAGE(Boost REQUIRED COMPONENTS  program_options iostreams system filesystem)
message(${Boost_PROGRAM_OPTIONS_LIBRARIES})
include_directories(${Boost_INCLUDE_DIRS})

FIND_PACKAGE(Eigen3 REQUIRED)
FIND_PACKAGE(Armadillo REQUIRED)
FIND_PACKAGE(Threads REQUIRED)
FIND_PACKAGE(OpenMP REQUIRED)
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
find_package(HDF5 REQUIRED COMPONENTS C CXX)
include_directories(${HDF5_INCLUDE_DIRS})
# find_package(Zstd REQUIRED)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
# set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/)
###########

# PLINK2.0 pgen library;
INCLUDE_DIRECTORIES("submodule/plink-ng/2.0")
INCLUDE_DIRECTORIES("submodule/plink-ng/2.0/simde")
##################################
# Compile local library
##################################
# add_compile_options(-Xclang -fopenmp)
## build library based on library type
add_library(${CUR_LIB_NAME} ${CUR_LIB_TYPE} ${SRCS}  )
## Link the local library to other libraries
target_link_libraries( ${CUR_LIB_NAME}
    ${Boost_PROGRAM_OPTIONS_LIBRARIES} 
    Eigen3::Eigen
    Boost::iostreams
    Boost::filesystem
    # julia
    # ${JULIA_LIBRARY}
    # Pgenlib
    OpenMP::OpenMP_CXX
    # ${OpenMP_CXX_FLAGS}
    ${CMAKE_THREAD_LIBS_INIT}
    ${ARMADILLO_LIBRARIES}
    ${BLAS_LIBRARIES}
    ${LAPACK_LIBRARIES}
    ${HDF5_LIBRARIES}
    )

##################################
# Compile the executable
##################################
add_executable(${targetName}64 ${MAIN_SOURCE})

# Link the executable to  libraries
target_link_libraries(${targetName}64 
    ${CUR_LIB_NAME}
    # ${Boost_PROGRAM_OPTIONS_LIBRARIES}
    # Eigen3::Eigen
    # ${OpenMP_CXX_FLAGS}
    # ${ARMADILLO_LIBRARIES}
    ${BLAS_LIBRARIES}
    ${LAPACK_LIBRARIES}
    ${HDF5_LIBRARIES} 
    # ${CMAKE_THREAD_LIBS_INIT} 
    )

##################################
#### How and what to install

# install(TARGETS ${CUR_LIB_NAME} LIBRARY DESTINATION lib ARCHIVE DESTINATION lib)
# ## install(TARGETS ${targetName}64 RUNTIME DESTINATION bin)
# install(TARGETS ${targetName}64  DESTINATION ${CMAKE_BINARY_DIR})
# install(DIRECTORY ../include/${targetName} DESTINATION include)


