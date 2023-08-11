project(odepack-demo LANGUAGES CXX)

add_executable(
    odepack-demo
    ${CMAKE_CURRENT_SOURCE_DIR}/test.cc
)
target_link_libraries(
    odepack-demo
    odepack::odepack
)