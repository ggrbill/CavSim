project (CavSim)

add_library(${PROJECT_NAME}
            STATIC
            Cavity.cpp
            Cavity.hpp
            solver.cpp
            solver.hpp
            IO.cpp
            IO.hpp
            Legacy/Structures.cpp
            Legacy/Structures.hpp)

add_executable (${PROJECT_NAME}Old main.cpp)

target_link_libraries(${PROJECT_NAME}Old ${PROJECT_NAME})

install(TARGETS ${PROJECT_NAME}Old ${PROJECT_NAME} 
        DESTINATION ${ARTIFACTS_DIR})
