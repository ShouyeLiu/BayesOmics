set(BENCHMARK_ENABLE_TESTING OFF CACHE BOOL "Suppressing benchmark's tests" FORCE)

FetchContent_Declare(
  googlebenchmark
  GIT_REPOSITORY    https://github.com/google/benchmark.git
  GIT_TAG           d572f4777349d43653b21d6c2fc63020ab326db2 # version 1.7.1
)

FetchContent_MakeAvailable(
  googlebenchmark
)

file(GLOB BENCHE_SRC_LIST ${PROJECT_SOURCE_DIR}/bench/*.cc)

set(COMMON_INCLUDES ${PROJECT_SOURCE_DIR}/include/)
include_directories(${COMMON_INCLUDES})

add_executable(bench_${targetName} ${BENCHE_SRC_LIST})

target_link_libraries(bench_${targetName}
  benchmark
  ${CUR_LIB_NAME}
)
