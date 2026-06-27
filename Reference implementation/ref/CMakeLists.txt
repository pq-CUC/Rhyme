cmake_minimum_required(VERSION 3.18)
project(rhyme LANGUAGES C)

enable_testing()
set(CMAKE_C_STANDARD 11)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/libs/)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/libs/)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(rhyme_SRCS
  src/poly.c
  src/ntt.c
  src/ntt_tables.c
  src/sampler.c
  src/encoding.c
  src/packing.c
  src/sign.c
  src/zpntt.c
  src/aes256ctr.c
  src/symmetric-shake.c
  src/keygen/kg_main.c
  src/keygen/kg_solver.c
  src/keygen/kg_zint.c
  src/keygen/kg_ntt.c
  src/keygen/kg_primes.c
)
set(FIPS202_SRCS src/fips202.c)

if(MSVC)
  set(C_FLAGS /nologo /O2 /W4)
else()
  set(C_FLAGS -O3 -fomit-frame-pointer -maes -mssse3 -Wall -Wextra)
  set(C_TEST_FLAGS -g -fsanitize=address,leak,undefined)
  set(LINK_TEST_FLAGS -fsanitize=address,leak,undefined)
endif()

include_directories(include)
include_directories(src/keygen)

add_library(fips202 STATIC ${FIPS202_SRCS})
add_library(rhyme_rng STATIC src/randombytes.c)

# one library per mode, plus a sanitizer variant
foreach(MODE 128 256 384 512)
  set(LIB rhyme${MODE})
  add_library(${LIB} STATIC ${rhyme_SRCS})
  target_compile_definitions(${LIB} PUBLIC RHYME_MODE=${MODE})
  target_compile_options(${LIB} PRIVATE ${C_FLAGS})
  target_link_libraries(${LIB} PUBLIC fips202 m)

  if(NOT MSVC)
    add_library(${LIB}_TEST STATIC ${rhyme_SRCS})
    target_compile_definitions(${LIB}_TEST PUBLIC RHYME_MODE=${MODE})
    target_compile_options(${LIB}_TEST PRIVATE ${C_FLAGS} ${C_TEST_FLAGS})
    target_link_libraries(${LIB}_TEST PUBLIC fips202 m)
    target_link_options(${LIB}_TEST PUBLIC ${LINK_TEST_FLAGS})
  endif()
endforeach()

add_subdirectory(test)
add_subdirectory(kat)
add_subdirectory(benchmark)
