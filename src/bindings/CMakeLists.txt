project (_CavSim)

find_package(pybind11 REQUIRED)
pybind11_add_module(
        ${PROJECT_NAME}
        CavityPython.cpp)

# Add include dirs of CavSim Library
target_include_directories(
        ${PROJECT_NAME} 
        PRIVATE 
                ../cpp)

# Link CavSim library to pybind11 module
target_link_libraries(
        ${PROJECT_NAME} 
        PRIVATE 
                CavSim)

install(TARGETS 
        ${PROJECT_NAME} 
        DESTINATION 
        ${ARTIFACTS_DIR})