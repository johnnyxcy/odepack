project(odepack-demo LANGUAGES CXX)

add_executable(
    odepack-demo
    ${CMAKE_CURRENT_SOURCE_DIR}/main.cc
)
target_link_libraries(
    odepack-demo
    odepack::odepack
)