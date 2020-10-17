project (CavSim)

add_library(${PROJECT_NAME}
            STATIC
            Cavity.cpp
            Cavity.hpp
            classic/solver.cpp
            classic/solver.hpp
            classic/numeric.cpp
            classic/numeric.hpp
            classic/IO.cpp
            classic/IO.hpp
            classic/Structures.cpp
            classic/Structures.hpp
            classic/WUDS.cpp
            classic/WUDS.hpp
            classic/PRIMEcorrection.cpp
            classic/PRIMEcorrection.hpp
            classic/MassEquation.cpp
            classic/MassEquation.hpp
            classic/MomentumEquation.cpp
            classic/MomentumEquation.hpp)

add_executable (${PROJECT_NAME}Old main.cpp)

target_link_libraries(${PROJECT_NAME}Old ${PROJECT_NAME})

install(TARGETS ${PROJECT_NAME}Old ${PROJECT_NAME} 
        DESTINATION ${ARTIFACTS_DIR})
