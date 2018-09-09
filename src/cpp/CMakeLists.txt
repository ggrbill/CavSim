project (CavSim)

file (GLOB ${PROJECT_NAME}_sources ${CMAKE_SOURCE_DIR}/src/cpp/*.cpp)
add_executable (${PROJECT_NAME} ${${PROJECT_NAME}_sources})
