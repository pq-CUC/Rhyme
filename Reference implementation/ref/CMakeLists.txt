cmake_minimum_required(VERSION 3.18)
project(ccc LANGUAGES C CXX) # CXX for the google test

enable_testing() # Enables running `ctest`

set(CMAKE_C_STANDARD 11)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/libs/)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/libs/)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin/)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(ccc_SRCS
  src/poly.c
  src/ntt.c
  src/encoding.c 
  src/reduce.c
  src/polyvec.c
  src/polymat.c
  src/sampler.c
  src/packing.c
  src/sign.c
  src/aes256ctr.c
)

set(ccc_FIPS202_SRCS
  ${ccc_SRCS}
  src/symmetric-shake.c
)
set(FIPS202_SRCS src/fips202.c)

if(MSVC)
  set(C_FLAGS /nologo /O2 /W4 /wd4146 /wd4244)
else()
  set(C_FLAGS -O3 -fomit-frame-pointer -march=native -Wall -Wextra -Wpedantic)
  set(C_TEST_FLAGS -g -fsanitize=address,leak,undefined,pointer-compare,pointer-subtract)
endif()

find_package(OpenSSL REQUIRED)

include_directories(include)
include_directories(api)
link_directories(libs)

add_library(fips202 SHARED ${FIPS202_SRCS})
add_library(RNG SHARED ${PROJECT_SOURCE_DIR}/src/randombytes.c)
target_compile_options(RNG PRIVATE -O3 -fomit-frame-pointer -fPIC)
target_link_libraries(RNG PUBLIC OpenSSL::Crypto)

set(LINK_TEST_FLAGS -fsanitize=address,leak,undefined,pointer-compare,pointer-subtract)

########################
#### WITH SANITIZER ####
########################
# ccc 2 SHAKE ONLY
set(LIB_NAME2_TEST ${PROJECT_NAME}2_TEST)
add_library(${LIB_NAME2_TEST} SHARED ${ccc_FIPS202_SRCS})
target_compile_definitions(${LIB_NAME2_TEST} PUBLIC ccc_MODE=2)
target_compile_options(${LIB_NAME2_TEST} PRIVATE ${C_FLAGS} ${C_TEST_FLAGS})
target_link_libraries(${LIB_NAME2_TEST} INTERFACE fips202)
target_link_libraries(${LIB_NAME2_TEST} PUBLIC RNG)
target_link_options(${LIB_NAME2_TEST} PUBLIC ${LINK_TEST_FLAGS})

# ccc 3 SHAKE ONLY
set(LIB_NAME3_TEST ${PROJECT_NAME}3_TEST)
add_library(${LIB_NAME3_TEST} SHARED ${ccc_FIPS202_SRCS})
target_compile_definitions(${LIB_NAME3_TEST} PUBLIC ccc_MODE=3)
target_compile_options(${LIB_NAME3_TEST} PRIVATE ${C_FLAGS} ${C_TEST_FLAGS})
target_link_libraries(${LIB_NAME3_TEST} INTERFACE fips202)
target_link_libraries(${LIB_NAME3_TEST} PUBLIC RNG)
target_link_options(${LIB_NAME3_TEST} PUBLIC ${LINK_TEST_FLAGS})

# ccc 5 SHAKE ONLY
set(LIB_NAME5_TEST ${PROJECT_NAME}5_TEST)
add_library(${LIB_NAME5_TEST} SHARED ${ccc_FIPS202_SRCS})
target_compile_definitions(${LIB_NAME5_TEST} PUBLIC ccc_MODE=5)
target_compile_options(${LIB_NAME5_TEST} PRIVATE ${C_FLAGS} ${C_TEST_FLAGS})
target_link_libraries(${LIB_NAME5_TEST} INTERFACE fips202)
target_link_libraries(${LIB_NAME5_TEST} PUBLIC RNG)
target_link_options(${LIB_NAME5_TEST} PUBLIC ${LINK_TEST_FLAGS})

###########################
#### WITHOUT SANITIZER ####
###########################
# ccc 2 SHAKE ONLY
set(LIB_NAME2 ${PROJECT_NAME}2)
add_library(${LIB_NAME2} SHARED ${ccc_FIPS202_SRCS})
target_compile_definitions(${LIB_NAME2} PUBLIC ccc_MODE=2)
target_compile_options(${LIB_NAME2} PRIVATE ${C_FLAGS})
target_link_libraries(${LIB_NAME2} INTERFACE fips202)
target_link_libraries(${LIB_NAME2} PUBLIC RNG)

# ccc 3 SHAKE ONLY
set(LIB_NAME3 ${PROJECT_NAME}3)
add_library(${LIB_NAME3} SHARED ${ccc_FIPS202_SRCS})
target_compile_definitions(${LIB_NAME3} PUBLIC ccc_MODE=3)
target_compile_options(${LIB_NAME3} PRIVATE ${C_FLAGS})
target_link_libraries(${LIB_NAME3} INTERFACE fips202)
target_link_libraries(${LIB_NAME3} PUBLIC RNG)

# ccc 5 SHAKE ONLY
set(LIB_NAME5 ${PROJECT_NAME}5)
add_library(${LIB_NAME5} SHARED ${ccc_FIPS202_SRCS})
target_compile_definitions(${LIB_NAME5} PUBLIC ccc_MODE=5)
target_compile_options(${LIB_NAME5} PRIVATE ${C_FLAGS})
target_link_libraries(${LIB_NAME5} INTERFACE fips202)
target_link_libraries(${LIB_NAME5} PUBLIC RNG)

add_subdirectory(test)
add_subdirectory(kat)
add_subdirectory(benchmark)
