cmake_minimum_required(VERSION 3.12)
set(ODEPACK_VERSION 0.1.0)
project(odepack LANGUAGES Fortran VERSION ${ODEPACK_VERSION})

find_package(LAPACK REQUIRED)


option(BUILDING_DEMO "Build demo" ON)

# include path
set(_ODEPACK_INCLUDE_DIR ${CMAKE_CURRENT_LIST_DIR}/include)
# src files
set(
    _ODEPACK_SRC
    src/dlsoda/dlsoda.f
    src/dlsoda/opkda1.f
    src/dlsoda/opkda2.f
    src/dverk/dverk.f
)
add_library(odepack SHARED ${_ODEPACK_SRC})
target_compile_options(
    odepack
    PRIVATE
    -std=legacy
    -fpic
    -fno-underscoring
)
target_link_libraries(
    odepack
    PRIVATE
    LAPACK::LAPACK
)

target_include_directories(
    odepack
    PRIVATE
    $<BUILD_INTERFACE:${_ODEPACK_INCLUDE_DIR}>
)
set(ODEPACK_LIBRARIES odepack)

# interface
add_library(odepack::odepack ALIAS odepack)
set_property(
    TARGET 
    odepack 
    PROPERTY 
    INTERFACE_INCLUDE_DIRECTORIES 
    $<BUILD_INTERFACE:${_ODEPACK_INCLUDE_DIR}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_PREFIX}/include/odepack>
)
set_property(
    TARGET 
    odepack
    PROPERTY 
    INTERFACE_LINK_LIBRARIES
    odepack
)
unset(_ODEPACK_INCLUDE_DIR)

if (BUILDING_DEMO)
    add_subdirectory(demo)
endif()

# # Tests and benchmark
# if (BUILD_TESTING)
#     message(STATUS "[odepack] Building test")
#     add_subdirectory(test)
# endif()

# if (BUILD_BENCHMARK)
#     message(STATUS "[odepack] Building benchmark")
#     add_subdirectory(benchmark)
# endif()

# 将配置文件中的变量设置为实际值
configure_file(odepack-config.cmake.in odepack-config.cmake @ONLY)

# export odepack
export(TARGETS odepack NAMESPACE odepack:: FILE odepackTargets.cmake)

# install config
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/odepack-config.cmake"
        DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/odepack")

# install target
install(EXPORT odepackTargets
        NAMESPACE odepack::
        DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/odepack")

# install header
install(DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/include
        DESTINATION "${CMAKE_INSTALL_PREFIX}"
        FILES_MATCHING PATTERN "*.h")

# install libraries
install(TARGETS odepack
        EXPORT odepackTargets
        ARCHIVE DESTINATION "${CMAKE_INSTALL_PREFIX}/lib"
        LIBRARY DESTINATION "${CMAKE_INSTALL_PREFIX}/lib"
        RUNTIME DESTINATION "${CMAKE_INSTALL_PREFIX}/bin"
        INCLUDES DESTINATION "${CMAKE_INSTALL_PREFIX}/include")
