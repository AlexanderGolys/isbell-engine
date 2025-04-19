# CMake entry point
cmake_minimum_required (VERSION 3.0)
project (Tutorials)

find_package(OpenGL REQUIRED)


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
cmake_policy(SET CMP0026 OLD)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/external/rpavlik-cmake-modules-fe2273")
include(CreateLaunchers)
include(MSVCMultipleProcessCompile) # /MP

if(INCLUDE_DISTRIB)
	add_subdirectory(distrib)
endif(INCLUDE_DISTRIB)

# c++ 23
set(CMAKE_CXX_STANDARD 23)



include_directories(
	external/AntTweakBar-1.16/include/
	external/glfw-3.1.2/include/
	external/glm-0.9.7.1/
	external/glew-1.13.0/include/
	external/assimp-3.0.1270/include/
	external/bullet-2.81-rev2613/src/
	.
)

set(ALL_LIBS
	${OPENGL_LIBRARY}
	glfw
	GLEW_1130
)

add_definitions(
	-DTW_STATIC
	-DTW_NO_LIB_PRAGMA
	-DTW_NO_DIRECT3D
	-DGLEW_STATIC
	-D_CRT_SECURE_NO_WARNINGS
)








add_executable(guitar
        src/guitar.cpp
		src/fundamentals/macros.hpp
        src/fundamentals/filesUtils.cpp
        src/fundamentals/mat.cpp
        src/fundamentals/func.cpp
        src/common/indexedRendering.cpp
        src/common/specific.cpp
        src/common/glsl_utils.cpp
        src/common/renderingUtils.cpp
        src/geometry/complexGeo.cpp
        src/geometry/discreteGeometry.cpp
#        src/geometry/hyperbolic.cpp
        src/geometry/planarGeometry.cpp
        src/geometry/smoothParametric.cpp
        src/geometry/smoothImplicit.cpp
        src/physics/solidMeshes.cpp
		src/physics/rigid.cpp
		src/fundamentals/flows.cpp
		src/fundamentals/flows.hpp
		src/fundamentals/modules.hpp
        src/fundamentals/integralTransforms.hpp
        src/fundamentals/integralTransforms.cpp
		src/geometry/pde.hpp
		src/geometry/pde.cpp
		src/common/interface.cpp
		src/common/dynamicalSystems.cpp
		src/fundamentals/elemFunc.cpp
		src/fundamentals/elemFunc.hpp
)

target_link_libraries(guitar
        ${ALL_LIBS}
)
# Xcode and Visual working directories
set_target_properties(guitar PROPERTIES XCODE_ATTRIBUTE_CONFIGURATION_BUILD_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src/")
create_target_launcher(guitar WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/src/")





add_executable(guitarCyg
        src/stringCyg.cpp
		src/fundamentals/macros.hpp
        src/fundamentals/filesUtils.cpp
        src/fundamentals/mat.cpp
        src/fundamentals/func.cpp
        src/common/indexedRendering.cpp
        src/common/specific.cpp
        src/common/glsl_utils.cpp
        src/common/renderingUtils.cpp
        src/geometry/complexGeo.cpp
        src/geometry/discreteGeometry.cpp
#        src/geometry/hyperbolic.cpp
        src/geometry/planarGeometry.cpp
        src/geometry/smoothParametric.cpp
        src/geometry/smoothImplicit.cpp
        src/physics/solidMeshes.cpp
		src/physics/rigid.cpp
		src/fundamentals/flows.cpp
		src/fundamentals/flows.hpp
		src/fundamentals/modules.hpp
        src/fundamentals/integralTransforms.hpp
        src/fundamentals/integralTransforms.cpp
		src/geometry/pde.hpp
		src/geometry/pde.cpp
		src/geometry/pde_dicrete.hpp
		src/geometry/pde_dicrete.cpp
		src/common/interface.cpp
		src/common/dynamicalSystems.cpp
		src/fundamentals/elemFunc.cpp
		src/fundamentals/elemFunc.hpp
		src/geometry/pde_dicrete.cpp
		src/geometry/pde_dicrete.hpp
)

target_link_libraries(guitarCyg
        ${ALL_LIBS}
)
# Xcode and Visual working directories
set_target_properties(guitarCyg PROPERTIES XCODE_ATTRIBUTE_CONFIGURATION_BUILD_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src/")
create_target_launcher(guitarCyg WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/src/")

#
#add_executable(harmonic
#        src/harmonicCompletion.cpp
#		src/fundamentals/macros.hpp
#        src/fundamentals/filesUtils.cpp
#        src/fundamentals/mat.cpp
#        src/fundamentals/func.cpp
#        src/common/indexedRendering.cpp
#        src/common/specific.cpp
#        src/common/glsl_utils.cpp
#        src/common/renderingUtils.cpp
#        src/geometry/complexGeo.cpp
#        src/geometry/discreteGeometry.cpp
##        src/geometry/hyperbolic.cpp
#        src/geometry/planarGeometry.cpp
#        src/geometry/smoothParametric.cpp
#        src/geometry/smoothImplicit.cpp
#        src/physics/solidMeshes.cpp
#		src/physics/rigid.cpp
#		src/fundamentals/flows.cpp
#		src/fundamentals/flows.hpp
#		src/fundamentals/modules.hpp
#        src/fundamentals/integralTransforms.hpp
#        src/fundamentals/integralTransforms.cpp
#		src/geometry/pde.hpp
#		src/geometry/pde.cpp
#
#		src/common/interface.cpp
#		src/common/dynamicalSystems.cpp
#		src/fundamentals/elemFunc.cpp
#		src/fundamentals/elemFunc.hpp
#)
#
#target_link_libraries(harmonic
#        ${ALL_LIBS}
#)
## Xcode and Visual working directories
#set_target_properties(harmonic PROPERTIES XCODE_ATTRIBUTE_CONFIGURATION_BUILD_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src/")
#create_target_launcher(harmonic WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/src/")
#
#
add_executable(tests1
        src/.tests/discrete_func_tests.cpp
		src/fundamentals/macros.hpp
        src/fundamentals/filesUtils.cpp
        src/fundamentals/mat.cpp
        src/fundamentals/func.cpp
        src/common/indexedRendering.cpp
        src/common/specific.cpp
        src/common/glsl_utils.cpp
        src/common/renderingUtils.cpp
        src/geometry/complexGeo.cpp
        src/geometry/discreteGeometry.cpp
#        src/geometry/hyperbolic.cpp
        src/geometry/planarGeometry.cpp
        src/geometry/smoothParametric.cpp
        src/geometry/smoothImplicit.cpp
        src/physics/solidMeshes.cpp
		src/physics/rigid.cpp
		src/fundamentals/flows.cpp
		src/fundamentals/flows.hpp
		src/fundamentals/modules.hpp
        src/fundamentals/integralTransforms.hpp
        src/fundamentals/integralTransforms.cpp
		src/geometry/pde.hpp
		src/geometry/pde.cpp
		src/geometry/pde_dicrete.hpp
		src/geometry/pde_dicrete.cpp
		src/common/interface.cpp
		src/common/dynamicalSystems.cpp
		src/fundamentals/elemFunc.cpp
		src/fundamentals/elemFunc.hpp
)

target_link_libraries(tests1
        ${ALL_LIBS}
)
# Xcode and Visual working directories
set_target_properties(tests1 PROPERTIES XCODE_ATTRIBUTE_CONFIGURATION_BUILD_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src/")
create_target_launcher(tests1 WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/src/")




add_executable(gabor
        src/gabor.cpp
		src/fundamentals/macros.hpp
        src/fundamentals/filesUtils.cpp
        src/fundamentals/mat.cpp
        src/fundamentals/func.cpp
        src/common/indexedRendering.cpp
        src/common/specific.cpp
        src/common/glsl_utils.cpp
        src/common/renderingUtils.cpp
        src/geometry/complexGeo.cpp
        src/geometry/discreteGeometry.cpp
#        src/geometry/hyperbolic.cpp
        src/geometry/planarGeometry.cpp
        src/geometry/smoothParametric.cpp
        src/geometry/smoothImplicit.cpp
        src/physics/solidMeshes.cpp
		src/physics/rigid.cpp
		src/fundamentals/flows.cpp
		src/fundamentals/flows.hpp
		src/fundamentals/modules.hpp
        src/fundamentals/integralTransforms.hpp
        src/fundamentals/integralTransforms.cpp
		src/geometry/pde.hpp
		src/geometry/pde.cpp
		src/geometry/pde_dicrete.hpp
		src/geometry/pde_dicrete.cpp
		src/common/interface.cpp
		src/common/dynamicalSystems.cpp
		src/fundamentals/elemFunc.cpp
		src/fundamentals/elemFunc.hpp
		src/geometry/pde_dicrete.cpp
		src/geometry/pde_dicrete.hpp
)

target_link_libraries(gabor
        ${ALL_LIBS}
)
# Xcode and Visual working directories
set_target_properties(gabor PROPERTIES XCODE_ATTRIBUTE_CONFIGURATION_BUILD_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src/")
create_target_launcher(gabor WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/src/")

#
