cmake_minimum_required(VERSION 3.16)
project(ckks_report_c VERSION 1.0 LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(SEAL 4.1 REQUIRED)

add_library(ckks_api STATIC ckks_api.cpp)
target_include_directories(ckks_api PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(ckks_api PRIVATE SEAL::seal)

if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/keygen_center.cpp")
	add_executable(keygen_center keygen_center.cpp)
	target_link_libraries(keygen_center PRIVATE SEAL::seal)
endif()

if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/encrypt_data.cpp")
	add_executable(encrypt_data encrypt_data.cpp)
	target_link_libraries(encrypt_data PRIVATE SEAL::seal)
endif()

if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/decrypt_data.cpp")
	add_executable(decrypt_data decrypt_data.cpp)
	target_link_libraries(decrypt_data PRIVATE SEAL::seal)
endif()

add_executable(keygen_center_c keygen_center_c.c)
target_link_libraries(keygen_center_c PRIVATE ckks_api)

add_executable(encrypt_data_c encrypt_data_c.c)
target_link_libraries(encrypt_data_c PRIVATE ckks_api)

add_executable(decrypt_data_c decrypt_data_c.c)
target_link_libraries(decrypt_data_c PRIVATE ckks_api)

add_executable(benchmark_c benchmark_c.c)
target_link_libraries(benchmark_c PRIVATE ckks_api)
