cmake_minimum_required(VERSION 3.25)

if(POLICY CMP0091)
    cmake_policy(SET CMP0091 NEW)
endif()

project(GameProject)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(MSVC)
    add_compile_options(/MP)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /WX- /W4 /O2")
    add_compile_options(/bigobj)
    set(CMAKE_MSVC_DEBUG_INFORMATION_FORMAT "Embedded")
    set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>DLL")
endif()

set(BUILD_STATIC_LIB ON)
set(TARGET_BUILD_PLATFORM "windows")
set(PLATFORM "Windows")
add_compile_definitions(_DEBUG)
add_definitions(-D_DEBUG)
add_definitions(-D_ITERATOR_DEBUG_LEVEL=2)
add_definitions(-DLIBTYPE=STATIC)
add_compile_definitions(IMGUI_DEFINE_MATH_OPERATORS)
option(MSVC_CRT_STATIC "Link primesieve.lib with /MT instead of the default /MD" OFF)
set(PLATFORM_CMAKELISTS "${CMAKE_CURRENT_SOURCE_DIR}/CMakeLists.txt")
set(CMAKE_BUILD_TYPE Debug CACHE STRING "Build type")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /Ob0 /Od /RTC1")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /wd4244")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /C2220")
set(PLAZA_ITERATOR_LEVEL 2)
set(ENGINE_SOURCE_DIR "C:/Plaza/Dev/PlazaEngine")

file(GLOB_RECURSE GAME_SOURCES CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/Assets/*.cpp" "${CMAKE_CURRENT_SOURCE_DIR}/Assets/*.h" "${CMAKE_CURRENT_SOURCE_DIR}/Scripts/*.cpp" "${CMAKE_CURRENT_SOURCE_DIR}/Scripts/*.h")
add_library(GameLib SHARED ${GAME_SOURCES})
target_precompile_headers(GameLib PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/Assets/Scripts/pch.h")
target_compile_definitions(GameLib PRIVATE ENGINE_IMPORT)

target_compile_definitions(GameLib PUBLIC PX_PHYSX_STATIC_LIB)
target_compile_definitions(GameLib PUBLIC EDITOR_MODE)
target_compile_definitions(GameLib PUBLIC COMPILING_GAME_DLL)
target_link_options(GameLib PUBLIC "/INCREMENTAL /PDBALTPATH:$<TARGET_PDB_FILE_NAME:GameLib>")
link_directories("${CMAKE_SOURCE_DIR}/bin")

set_property(TARGET GameLib PROPERTY CXX_STANDARD 20)

file(GLOB_RECURSE ENGINE_SOURCES
    "${ENGINE_SOURCE_DIR}/Core/*.cpp"
    "${ENGINE_SOURCE_DIR}/Core/*.h"
    "${ENGINE_SOURCE_DIR}/Modules/*.cpp"
    "${ENGINE_SOURCE_DIR}/Modules/*.h"
)

target_sources(GameLib PRIVATE ${ENGINE_SOURCES})

target_include_directories(GameLib PRIVATE
    "${ENGINE_SOURCE_DIR}/Core"
    "${ENGINE_SOURCE_DIR}/Modules"
)

include(FetchContent)
FetchContent_Declare(
    PlazaEngineCoreLibs
    SOURCE_DIR "${ENGINE_SOURCE_DIR}/"
)

FetchContent_MakeAvailable(PlazaEngineCoreLibs)

target_link_libraries(GameLib PUBLIC PlazaEngineCore PlazaEngineCoreLibs)

target_include_directories(GameLib PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/"
    "${CMAKE_CURRENT_SOURCE_DIR}/Assets/Scripts/"
    "${CMAKE_CURRENT_SOURCE_DIR}/Scripts/"
)

set(OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin")
add_custom_command(TARGET GameLib POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E make_directory "${OUTPUT_DIRECTORY}"
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
    "$<TARGET_FILE:GameLib>"
    "$<TARGET_PDB_FILE:GameLib>"
    "${OUTPUT_DIRECTORY}"
)

add_custom_target(CopyGameLib ALL DEPENDS GameLib)
add_custom_target(DeleteGamePDB ALL DEPENDS GameLib)

add_custom_command(TARGET GameLib
    PRE_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy "${CMAKE_CURRENT_SOURCE_DIR}/bin/PlazaEngineCore.lib" "${CMAKE_CURRENT_BINARY_DIR}/"
    COMMENT "Copying PlazaEngineCore.lib to game project folder")
