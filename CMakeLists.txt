cmake_minimum_required(VERSION 3.20)

project(charged_particle_sim
    VERSION 0.1.0
    LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_executable(charged_particle_sim
    src/main.cpp
    src/Fields.cpp
    src/Integrator.cpp
)

target_include_directories(charged_particle_sim
    PRIVATE
        ${PROJECT_SOURCE_DIR}/include
)

if(MSVC)
    target_compile_options(charged_particle_sim PRIVATE /W4)
else()
    target_compile_options(charged_particle_sim PRIVATE -Wall -Wextra -Wpedantic)
endif()
