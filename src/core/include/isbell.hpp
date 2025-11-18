#pragma once

// Utils
#include "utils/abstractNonsense.hpp"
#include "utils/concepts.hpp"
#include "utils/elemFunc.hpp"
#include "utils/formatters.hpp"
#include "utils/integralTransforms.hpp"
#include "utils/interfaces.hpp"
#include "utils/logging.hpp"
#include "utils/modules.hpp"
#include "utils/solvers.hpp"
#include "utils/stringUtils.hpp"
#include "utils/symbolic.hpp"
#include "utils/symbolicAlgebras.hpp"

// Controllers
#include "controllers/cameraControllers.hpp"

// Events
#include "events/event.hpp"
#include "events/eventQueue.hpp"
#include "events/keyCodes.hpp"

// File Management
#include "file-management/configFiles.hpp"
#include "file-management/filesUtils.hpp"
#include "file-management/macroParsing.hpp"
#include "file-management/shaderGenerator.hpp"

// Geometry
#include "geometry/discreteGeometry.hpp"
#include "geometry/hyperbolic.hpp"
#include "geometry/hyperbolic_dep.hpp"
#include "geometry/pde.hpp"
#include "geometry/pdeDiscrete.hpp"
#include "geometry/planarGeometry.hpp"
#include "geometry/smoothParametric.hpp"
#include "geometry/sph.hpp"

// Physics
#include "physics/generic.hpp"
#include "physics/rigid.hpp"
#include "physics/solidMeshes.hpp"

// Engine
#include "engine/buffers.hpp"
#include "engine/clock.hpp"
#include "engine/glCommand.hpp"
#include "engine/indexedMesh.hpp"
#include "engine/planarMesh.hpp"
#include "engine/renderer.hpp"
#include "engine/shaders.hpp"
#include "engine/shaderTypes.hpp"
#include "engine/UILayers.hpp"

// SDF Rendering
#include "engine/sdf-rendering/SDFObjects.hpp"
#include "engine/sdf-rendering/SDFObjects2.hpp"
#include "engine/sdf-rendering/SDFRendering.hpp"
