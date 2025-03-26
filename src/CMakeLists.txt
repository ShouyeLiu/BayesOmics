

set(COMMON_INCLUDES ${PROJECT_SOURCE_DIR}/include)
set(MAIN_SOURCE "${PROJECT_SOURCE_DIR}/src/Main.cpp")
include_directories(${COMMON_INCLUDES})


file(GLOB SRCS "${PROJECT_SOURCE_DIR}/src/*.cpp")
list(REMOVE_ITEM SRCS "${MAIN_SOURCE}")

##################################
# set(STATLIB_DIR "/Users/uqsliu28/.local/lib/stats/include/")
# include_directories(${STATLIB_DIR})
# set(GCEM_DIR "/Users/uqsliu28/.local/lib/gcem/include/")
# include_directories(${GCEM_DIR})
## use julia 
## Check operating system
# if(APPLE)
#     # macOS specific Julia path
#     set(JULIA_INCLUDE "/opt/homebrew/Cellar/julia/1.9.4/include/julia")
#     set(JULIA_LIBRARY_DIR "/opt/homebrew/Cellar/julia/1.9.4/lib/")
#     set(JULIA_LIBRARY "/opt/homebrew/Cellar/julia/1.9.4/lib/libjulia.dylib") 
# elseif(UNIX AND NOT APPLE)
#     # Linux specific Julia path
#     set(JULIA_INCLUDE "/QRISdata/Q4062/softwares/local/miniconda3/include/julia")
#     set(JULIA_LIBRARY_DIR "/QRISdata/Q4062/softwares/local/miniconda3/lib/")
#     set(JULIA_LIBRARY "/QRISdata/Q4062/softwares/local/miniconda3/lib/libjulia.so") 
# endif()
# find_package(Julia REQUIRED)
include_directories(${JULIA_INCLUDE})
link_directories(${JULIA_LIBRARY_DIR})
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
# find_package(Zstd REQUIRED)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
# set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/)
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
    OpenMP::OpenMP_CXX
    # ${OpenMP_CXX_FLAGS}
    ${CMAKE_THREAD_LIBS_INIT}
    ${ARMADILLO_LIBRARIES}
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
    # ${CMAKE_THREAD_LIBS_INIT} 
    )

##################################
#### How and what to install

# install(TARGETS ${CUR_LIB_NAME} LIBRARY DESTINATION lib ARCHIVE DESTINATION lib)
# ## install(TARGETS ${targetName}64 RUNTIME DESTINATION bin)
# install(TARGETS ${targetName}64  DESTINATION ${CMAKE_BINARY_DIR})
# install(DIRECTORY ../include/${targetName} DESTINATION include)


