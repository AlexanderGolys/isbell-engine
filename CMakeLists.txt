cmake_minimum_required (VERSION 3.0)
project (Engine)

set(CMAKE_CXX_STANDARD 26)
set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Release>:Release>")

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
add_compile_options(-Wno-literal-suffix)

if(INCLUDE_DISTRIB)
	add_subdirectory(distrib)
endif(INCLUDE_DISTRIB)

# c++ 26



include_directories(
#	external/AntTweakBar-1.16/include/
	external/glfw-3.1.2/include/
	external/glm-0.9.7.1/
	external/glew-1.13.0/include/
#	external/assimp-3.0.1270/include/
#	external/bullet-2.81-rev2613/src/
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
		src/render-projects/guitar.cpp
		src/utils/macros.hpp
		src/utils/filesUtils.cpp
		src/utils/mat.cpp
		src/utils/func.cpp
		src/engine/indexedRendering.cpp
		src/engine/specific.cpp
		src/engine/glslUtils.cpp
		src/engine/renderingUtils.cpp
        src/geometry/complexGeo.cpp
        src/geometry/discreteGeometry.cpp
#        src/geometry/hyperbolic.cpp
        src/geometry/planarGeometry.cpp
        src/geometry/smoothParametric.cpp
        src/geometry/smoothImplicit.cpp
        src/physics/solidMeshes.cpp
		src/physics/rigid.cpp
		src/utils/flows.cpp
		src/utils/flows.hpp
		src/utils/modules.hpp
		src/utils/integralTransforms.hpp
		src/utils/integralTransforms.cpp
		src/geometry/pde.hpp
		src/geometry/pde.cpp
		src/engine/interface.cpp
		src/engine/dynamicalSystems.cpp
		src/utils/elemFunc.cpp
		src/utils/elemFunc.hpp
)

target_link_libraries(guitar
        ${ALL_LIBS}
)
# Xcode and Visual working directories
set_target_properties(guitar PROPERTIES XCODE_ATTRIBUTE_CONFIGURATION_BUILD_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src/")
create_target_launcher(guitar WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/src/")





add_executable(guitarCyg
		src/render-projects/stringCyg.cpp
		src/utils/macros.hpp
		src/utils/filesUtils.cpp
		src/utils/mat.cpp
		src/utils/func.cpp
		src/engine/indexedRendering.cpp
		src/engine/specific.cpp
		src/engine/glslUtils.cpp
		src/engine/renderingUtils.cpp
        src/geometry/complexGeo.cpp
        src/geometry/discreteGeometry.cpp
#        src/geometry/hyperbolic.cpp
        src/geometry/planarGeometry.cpp
        src/geometry/smoothParametric.cpp
        src/geometry/smoothImplicit.cpp
        src/physics/solidMeshes.cpp
		src/physics/rigid.cpp
		src/utils/flows.cpp
		src/utils/flows.hpp
		src/utils/modules.hpp
		src/utils/integralTransforms.hpp
		src/utils/integralTransforms.cpp
		src/geometry/pde.hpp
		src/geometry/pde.cpp
		src/geometry/pdeDiscrete.hpp
		src/geometry/pdeDiscrete.cpp
		src/engine/interface.cpp
		src/engine/dynamicalSystems.cpp
		src/utils/elemFunc.cpp
		src/utils/elemFunc.hpp
		src/geometry/pdeDiscrete.cpp
		src/geometry/pdeDiscrete.hpp
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
#		src/utils/macros.hpp
#        src/utils/filesUtils.cpp
#        src/utils/mat.cpp
#        src/utils/func.cpp
#        src/engine/indexedRendering.cpp
#        src/engine/specific.cpp
#        src/engine/glslUtils.cpp
#        src/engine/renderingUtils.cpp
#        src/geometry/complexGeo.cpp
#        src/geometry/discreteGeometry.cpp
##        src/geometry/hyperbolic.cpp
#        src/geometry/planarGeometry.cpp
#        src/geometry/smoothParametric.cpp
#        src/geometry/smoothImplicit.cpp
#        src/physics/solidMeshes.cpp
#		src/physics/rigid.cpp
#		src/utils/flows.cpp
#		src/utils/flows.hpp
#		src/utils/modules.hpp
#        src/utils/integralTransforms.hpp
#        src/utils/integralTransforms.cpp
#		src/geometry/pde.hpp
#		src/geometry/pde.cpp
#
#		src/engine/interface.cpp
#		src/engine/dynamicalSystems.cpp
#		src/utils/elemFunc.cpp
#		src/utils/elemFunc.hpp
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
		src/tests/discreteFuncTests.cpp
		src/utils/macros.hpp
		src/utils/filesUtils.cpp
		src/utils/mat.cpp
		src/utils/func.cpp
		src/engine/indexedRendering.cpp
		src/engine/specific.cpp
		src/engine/glslUtils.cpp
		src/engine/renderingUtils.cpp
        src/geometry/complexGeo.cpp
        src/geometry/discreteGeometry.cpp
#        src/geometry/hyperbolic.cpp
        src/geometry/planarGeometry.cpp
        src/geometry/smoothParametric.cpp
        src/geometry/smoothImplicit.cpp
        src/physics/solidMeshes.cpp
		src/physics/rigid.cpp
		src/utils/flows.cpp
		src/utils/flows.hpp
		src/utils/modules.hpp
		src/utils/integralTransforms.hpp
		src/utils/integralTransforms.cpp
		src/geometry/pde.hpp
		src/geometry/pde.cpp
		src/geometry/pdeDiscrete.hpp
		src/geometry/pdeDiscrete.cpp
		src/engine/interface.cpp
		src/engine/dynamicalSystems.cpp
		src/utils/elemFunc.cpp
		src/utils/elemFunc.hpp
)

target_link_libraries(tests1
        ${ALL_LIBS}
)
# Xcode and Visual working directories
set_target_properties(tests1 PROPERTIES XCODE_ATTRIBUTE_CONFIGURATION_BUILD_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src/")
create_target_launcher(tests1 WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/src/")




add_executable(gabor
		src/render-projects/gabor.cpp
		src/utils/macros.hpp
		src/utils/filesUtils.cpp
		src/utils/mat.cpp
		src/utils/func.cpp
		src/engine/indexedRendering.cpp
		src/engine/specific.cpp
		src/engine/glslUtils.cpp
		src/engine/renderingUtils.cpp
        src/geometry/complexGeo.cpp
        src/geometry/discreteGeometry.cpp
#        src/geometry/hyperbolic.cpp
        src/geometry/planarGeometry.cpp
        src/geometry/smoothParametric.cpp
        src/geometry/smoothImplicit.cpp
        src/physics/solidMeshes.cpp
		src/physics/rigid.cpp
		src/utils/flows.cpp
		src/utils/flows.hpp
		src/utils/modules.hpp
		src/utils/integralTransforms.hpp
		src/utils/integralTransforms.cpp
		src/geometry/pde.hpp
		src/geometry/pde.cpp
		src/geometry/pdeDiscrete.hpp
		src/geometry/pdeDiscrete.cpp
		src/engine/interface.cpp
		src/engine/dynamicalSystems.cpp
		src/utils/elemFunc.cpp
		src/utils/elemFunc.hpp
		src/geometry/pdeDiscrete.cpp
		src/geometry/pdeDiscrete.hpp
)

target_link_libraries(gabor
        ${ALL_LIBS}
)
# Xcode and Visual working directories
set_target_properties(gabor PROPERTIES XCODE_ATTRIBUTE_CONFIGURATION_BUILD_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src/")
create_target_launcher(gabor WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/src/")




add_executable(zero
		src/render-projects/zero.cpp
		src/utils/macros.hpp
		src/utils/filesUtils.cpp
		src/utils/mat.cpp
		src/utils/func.cpp
		src/engine/indexedRendering.cpp
		src/engine/specific.cpp
		src/engine/glslUtils.cpp
		src/engine/renderingUtils.cpp
        src/geometry/complexGeo.cpp
        src/geometry/discreteGeometry.cpp
#        src/geometry/hyperbolic.cpp
        src/geometry/planarGeometry.cpp
        src/geometry/smoothParametric.cpp
        src/geometry/smoothImplicit.cpp
        src/physics/solidMeshes.cpp
		src/physics/rigid.cpp
		src/utils/flows.cpp
		src/utils/flows.hpp
		src/utils/modules.hpp
		src/utils/integralTransforms.hpp
		src/utils/integralTransforms.cpp
		src/geometry/pde.hpp
		src/geometry/pde.cpp
		src/geometry/pdeDiscrete.hpp
		src/geometry/pdeDiscrete.cpp
		src/engine/interface.cpp
		src/engine/dynamicalSystems.cpp
		src/utils/elemFunc.cpp
		src/utils/elemFunc.hpp
		src/geometry/pdeDiscrete.cpp
		src/geometry/pdeDiscrete.hpp
)

target_link_libraries(zero
        ${ALL_LIBS}
)
# Xcode and Visual working directories
set_target_properties(zero PROPERTIES XCODE_ATTRIBUTE_CONFIGURATION_BUILD_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src/")
create_target_launcher(zero WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/src/")

#

#


add_executable(fluids
		src/render-projects/fluids.cpp
		src/utils/macros.hpp
		src/utils/filesUtils.cpp
		src/utils/mat.cpp
		src/utils/func.cpp
		src/engine/indexedRendering.cpp
		src/engine/specific.cpp
		src/engine/glslUtils.cpp
		src/engine/renderingUtils.cpp
        src/geometry/complexGeo.cpp
        src/geometry/discreteGeometry.cpp
        src/geometry/sph.cpp
        src/geometry/planarGeometry.cpp
        src/geometry/smoothParametric.cpp
        src/geometry/smoothImplicit.cpp
        src/physics/solidMeshes.cpp
		src/physics/rigid.cpp
		src/utils/flows.cpp
		src/utils/flows.hpp
		src/utils/modules.hpp
		src/utils/integralTransforms.hpp
		src/utils/integralTransforms.cpp
		src/geometry/pde.hpp
		src/geometry/pde.cpp
		src/geometry/pdeDiscrete.hpp
		src/geometry/pdeDiscrete.cpp
		src/engine/interface.cpp
		src/engine/dynamicalSystems.cpp
		src/utils/elemFunc.cpp
		src/utils/elemFunc.hpp
		src/geometry/pdeDiscrete.cpp
		src/geometry/pdeDiscrete.hpp
)

target_link_libraries(fluids
        ${ALL_LIBS}
)
# Xcode and Visual working directories
set_target_properties(fluids PROPERTIES XCODE_ATTRIBUTE_CONFIGURATION_BUILD_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src/")
create_target_launcher(fluids WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/src/")

#

#

add_executable(foliatedTorus
		src/render-projects/foliatedTorus.cpp
		src/utils/macros.hpp
		src/utils/filesUtils.cpp
		src/utils/mat.cpp
		src/utils/func.cpp
		src/engine/indexedRendering.cpp
		src/engine/specific.cpp
		src/engine/glslUtils.cpp
		src/engine/renderingUtils.cpp
        src/geometry/complexGeo.cpp
        src/geometry/discreteGeometry.cpp
#        src/geometry/hyperbolic.cpp
        src/geometry/planarGeometry.cpp
        src/geometry/smoothParametric.cpp
        src/geometry/smoothImplicit.cpp
        src/physics/solidMeshes.cpp
		src/physics/rigid.cpp
		src/utils/flows.cpp
		src/utils/flows.hpp
		src/utils/modules.hpp
		src/utils/integralTransforms.hpp
		src/utils/integralTransforms.cpp
		src/geometry/pde.hpp
		src/geometry/pde.cpp
		src/geometry/pdeDiscrete.hpp
		src/geometry/pdeDiscrete.cpp
		src/engine/interface.cpp
		src/engine/dynamicalSystems.cpp
		src/utils/elemFunc.cpp
		src/utils/elemFunc.hpp
		src/geometry/pdeDiscrete.cpp
		src/geometry/pdeDiscrete.hpp
)

target_link_libraries(foliatedTorus
        ${ALL_LIBS}
)
# Xcode and Visual working directories
set_target_properties(foliatedTorus PROPERTIES XCODE_ATTRIBUTE_CONFIGURATION_BUILD_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src/")
create_target_launcher(foliatedTorus WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/src/")

#

add_executable(zerodiv
		src/render-projects/zerodiv.cpp
		src/utils/macros.hpp
		src/utils/filesUtils.cpp
		src/utils/mat.cpp
		src/utils/func.cpp
		src/engine/indexedRendering.cpp
		src/engine/specific.cpp
		src/engine/glslUtils.cpp
		src/engine/renderingUtils.cpp
        src/geometry/complexGeo.cpp
        src/geometry/discreteGeometry.cpp
#        src/geometry/hyperbolic.cpp
        src/geometry/planarGeometry.cpp
        src/geometry/smoothParametric.cpp
        src/geometry/smoothImplicit.cpp
        src/physics/solidMeshes.cpp
		src/physics/rigid.cpp
		src/utils/flows.cpp
		src/utils/flows.hpp
		src/utils/modules.hpp
		src/utils/integralTransforms.hpp
		src/utils/integralTransforms.cpp
		src/geometry/pde.hpp
		src/geometry/pde.cpp
		src/geometry/pdeDiscrete.hpp
		src/geometry/pdeDiscrete.cpp
		src/engine/interface.cpp
		src/engine/dynamicalSystems.cpp
		src/utils/elemFunc.cpp
		src/utils/elemFunc.hpp
		src/geometry/pdeDiscrete.cpp
		src/geometry/pdeDiscrete.hpp
)

target_link_libraries(zerodiv
        ${ALL_LIBS}
)
# Xcode and Visual working directories
set_target_properties(zerodiv PROPERTIES XCODE_ATTRIBUTE_CONFIGURATION_BUILD_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src/")
create_target_launcher(zerodiv WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/src/")


add_executable(pallettes
		src/render-projects/palettes.cpp
		src/utils/macros.hpp
		src/utils/filesUtils.cpp
		src/utils/mat.cpp
		src/utils/func.cpp
		src/engine/indexedRendering.cpp
		src/engine/specific.cpp
		src/engine/glslUtils.cpp
		src/engine/renderingUtils.cpp
        src/geometry/complexGeo.cpp
        src/geometry/discreteGeometry.cpp
#        src/geometry/hyperbolic.cpp
        src/geometry/planarGeometry.cpp
        src/geometry/smoothParametric.cpp
        src/geometry/smoothImplicit.cpp
        src/physics/solidMeshes.cpp
		src/physics/rigid.cpp
		src/utils/flows.cpp
		src/utils/flows.hpp
		src/utils/modules.hpp
		src/utils/integralTransforms.hpp
		src/utils/integralTransforms.cpp
		src/geometry/pde.hpp
		src/geometry/pde.cpp
		src/geometry/pdeDiscrete.hpp
		src/geometry/pdeDiscrete.cpp
		src/engine/interface.cpp
		src/engine/dynamicalSystems.cpp
		src/utils/elemFunc.cpp
		src/utils/elemFunc.hpp
		src/geometry/pdeDiscrete.cpp
		src/geometry/pdeDiscrete.hpp
)

target_link_libraries(pallettes
        ${ALL_LIBS}
)
# Xcode and Visual working directories
set_target_properties(pallettes PROPERTIES XCODE_ATTRIBUTE_CONFIGURATION_BUILD_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src/")
create_target_launcher(pallettes WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/src/")

#

add_executable(dom
		src/render-projects/domains.cpp
		src/utils/macros.hpp
		src/utils/filesUtils.cpp
		src/utils/mat.cpp
		src/utils/func.cpp
		src/engine/indexedRendering.cpp
		src/engine/specific.cpp
		src/engine/glslUtils.cpp
		src/engine/renderingUtils.cpp
        src/geometry/complexGeo.cpp
        src/geometry/discreteGeometry.cpp
#        src/geometry/hyperbolic.cpp
        src/geometry/planarGeometry.cpp
        src/geometry/smoothParametric.cpp
        src/geometry/smoothImplicit.cpp
        src/physics/solidMeshes.cpp
		src/physics/rigid.cpp
		src/utils/flows.cpp
		src/utils/flows.hpp
		src/utils/modules.hpp
		src/utils/integralTransforms.hpp
		src/utils/integralTransforms.cpp
		src/geometry/pde.hpp
		src/geometry/pde.cpp
		src/geometry/pdeDiscrete.hpp
		src/geometry/pdeDiscrete.cpp
		src/engine/interface.cpp
		src/engine/dynamicalSystems.cpp
		src/utils/elemFunc.cpp
		src/utils/elemFunc.hpp
		src/geometry/pdeDiscrete.cpp
		src/geometry/pdeDiscrete.hpp
)

target_link_libraries(dom
        ${ALL_LIBS}
)
# Xcode and Visual working directories
set_target_properties(dom PROPERTIES XCODE_ATTRIBUTE_CONFIGURATION_BUILD_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src/")
create_target_launcher(dom WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/src/")
