#pragma once
#include "fluids.hpp"

IndexedHexahedralVolumetricMesh volumetricRing(vec3 center, vec3 normal, float radiusSmall, float radiusBig, float height, int radial_res, int vertical_res, int extrude_res);
IndexedHexahedralVolumetricMesh volumetricCylinder(vec3 center, vec3 normal, float radius, float height, int radial_res, int vertical_res, int extrude_res);
IndexedHexahedralVolumetricMesh primitiveHex(std::array<vec3, 8> vertices);
