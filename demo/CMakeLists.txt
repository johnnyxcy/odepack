project(odepack-test LANGUAGES CXX)

add_executable(
    test-odepack
    ${CMAKE_CURRENT_SOURCE_DIR}/test.cc
)
target_link_libraries(
    test-odepack
    odepack::odepack
)
add_test(
    NAME test-odepack
    COMMAND $<TARGET_FILE:test-odepack>
)
