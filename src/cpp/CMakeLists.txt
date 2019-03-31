project (CavSim)

add_library(_${PROJECT_NAME}
            STATIC
            Cavity.cpp
            Cavity.hpp
            solver.cpp
            solver.hpp)

add_executable (${PROJECT_NAME} main.cpp)

target_link_libraries(${PROJECT_NAME} _${PROJECT_NAME})

install(TARGETS ${PROJECT_NAME} _${PROJECT_NAME} 
        DESTINATION ${ARTIFACTS_DIR})

#file (GLOB ${PROJECT_NAME}_sources ${CMAKE_SOURCE_DIR}/src/cpp/*.cpp)
