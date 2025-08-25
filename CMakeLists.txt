cmake_minimum_required (VERSION 3.5)
project (Engine)

set(CMAKE_CXX_STANDARD 26)
set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Release>:Release>")

find_package(OpenGL REQUIRED)
#find_package(spdlog REQUIRED)

if( CMAKE_BINARY_DIR STREQUAL CMAKE_SOURCE_DIR )
    message( FATAL_ERROR "Please select another Build Directory ! (and give it a clever name, like bin_Visual2012_64bits/)" )
endif()
if( CMAKE_SOURCE_DIR MATCHES " " )
	message( "Your Source Directory contains spaces. If you experience problems when compiling, this can be the cause." )
endif()
if( CMAKE_BINARY_DIR MATCHES " " )
	message( "Your Build Directory contains spaces. If you experience problems when compiling, this can be the cause." )
endif()



# Compile external dependencies 
add_subdirectory (external)

# On Visual 2005 and above, this module can set the debug working directory
#cmake_policy(SET CMP0026 OLD)
add_compile_options(-Wno-literal-suffix)


include_directories(
	external/glfw-3.1.2/include/
	external/glm-0.9.7.1/
	external/glew-2.1.0/include/
	external/spdlog-1.x/include/spdlog
	.
)

set(ALL_LIBS
	${OPENGL_LIBRARY}
	glfw
	GLEW_1130
	spdlog
)

add_definitions(
	-DTW_STATIC
	-DTW_NO_LIB_PRAGMA
	-DTW_NO_DIRECT3D
	-DGLEW_STATIC
	-D_CRT_SECURE_NO_WARNINGS
	-DSPDLOG_COMPILED_LIB
)


set(CORE_UTILS
		src/utils/filesUtils.cpp
		src/utils/mat.cpp
		src/utils/func.cpp
		src/utils/flows.cpp
		src/utils/integralTransforms.cpp
		src/utils/elemFunc.cpp
		src/utils/logging.cpp
		src/utils/solvers.cpp
#		src/utils/symbolic.cpp
)

set(CORE_GEOMETRY
        src/geometry/complexGeo.cpp
        src/geometry/discreteGeometry.cpp
#        src/geometry/hyperbolic.cpp
        src/geometry/planarGeometry.cpp
        src/geometry/smoothParametric.cpp
        src/geometry/smoothImplicit.cpp
		src/geometry/pde.cpp
		src/geometry/pdeDiscrete.cpp
		src/geometry/sph.cpp
)

set(CORE_PHYSICS
        src/physics/solidMeshes.cpp
		src/physics/rigid.cpp
#		src/physics/generic.cpp
)

set(CORE_ENGINE
		src/engine/indexedRendering.cpp
#		src/engine/bufferUtils.cpp
		src/engine/specific.cpp
		src/engine/glslUtils.cpp
		src/engine/renderingUtils.cpp
		src/engine/interface.cpp
		src/engine/dynamicalSystems.cpp
		src/engine/sdf-rendering/SDFObjects.cpp
		src/engine/sdf-rendering/SDFRendering.cpp
)


add_executable(zerodiv
		src/render-projects/zerodiv.cpp
		${CORE_UTILS}
		${CORE_GEOMETRY}
		${CORE_PHYSICS}
		${CORE_ENGINE}
)

target_link_libraries(zerodiv
        ${ALL_LIBS}
)
# Xcode and Visual working directories
set_target_properties(zerodiv PROPERTIES XCODE_ATTRIBUTE_CONFIGURATION_BUILD_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src/")
#create_target_launcher(zerodiv WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/src/")
