
set(kat_srcs
    kat_sign.c
)

find_package(OpenSSL REQUIRED)

# --- Rank-2 SIGN (CCC Mode 2) ---
set(KAT_SIGN_LIB_NAME2 ${PROJECT_NAME}2-kat-sign)
add_library(${KAT_SIGN_LIB_NAME2} SHARED ${kat_srcs})
target_include_directories(${KAT_SIGN_LIB_NAME2} PUBLIC ../include)
target_compile_definitions(${KAT_SIGN_LIB_NAME2} PUBLIC ccc_MODE=2)
target_compile_options(${KAT_SIGN_LIB_NAME2} PRIVATE ${C_TEST_FLAGS}) 
target_link_libraries(${KAT_SIGN_LIB_NAME2} PUBLIC ${LIB_NAME2_TEST})
target_link_options(${KAT_SIGN_LIB_NAME2} PUBLIC ${LINK_TEST_FLAGS}) 


add_executable(${KAT_SIGN_LIB_NAME2}-main PQCgenKAT_sign.c)
target_link_libraries(${KAT_SIGN_LIB_NAME2}-main ${KAT_SIGN_LIB_NAME2} ${LIB_NAME2_TEST} OpenSSL::Crypto)

# --- Rank-3 SIGN (CCC Mode 3) ---
set(KAT_SIGN_LIB_NAME3 ${PROJECT_NAME}3-kat-sign)
add_library(${KAT_SIGN_LIB_NAME3} SHARED ${kat_srcs})
target_include_directories(${KAT_SIGN_LIB_NAME3} PUBLIC ../include)
target_compile_definitions(${KAT_SIGN_LIB_NAME3} PUBLIC ccc_MODE=3)
target_compile_options(${KAT_SIGN_LIB_NAME3} PRIVATE ${C_TEST_FLAGS}) 
target_link_libraries(${KAT_SIGN_LIB_NAME3} PUBLIC ${LIB_NAME3_TEST})
target_link_options(${KAT_SIGN_LIB_NAME3} PUBLIC ${LINK_TEST_FLAGS}) 

add_executable(${KAT_SIGN_LIB_NAME3}-main PQCgenKAT_sign.c)
target_link_libraries(${KAT_SIGN_LIB_NAME3}-main ${KAT_SIGN_LIB_NAME3} ${LIB_NAME3_TEST} OpenSSL::Crypto)

# --- Rank-5 SIGN (CCC Mode 5) ---
set(KAT_SIGN_LIB_NAME5 ${PROJECT_NAME}5-kat-sign)
add_library(${KAT_SIGN_LIB_NAME5} SHARED ${kat_srcs})
target_include_directories(${KAT_SIGN_LIB_NAME5} PUBLIC ../include)
target_compile_definitions(${KAT_SIGN_LIB_NAME5} PUBLIC ccc_MODE=5)
target_compile_options(${KAT_SIGN_LIB_NAME5} PRIVATE ${C_TEST_FLAGS}) 
target_link_libraries(${KAT_SIGN_LIB_NAME5} PUBLIC ${LIB_NAME5_TEST})
target_link_options(${KAT_SIGN_LIB_NAME5} PUBLIC ${LINK_TEST_FLAGS})  

add_executable(${KAT_SIGN_LIB_NAME5}-main PQCgenKAT_sign.c)
target_link_libraries(${KAT_SIGN_LIB_NAME5}-main ${KAT_SIGN_LIB_NAME5} ${LIB_NAME5_TEST} OpenSSL::Crypto)