#pragma once
#include <functional>

#include <iosfwd>
#include <map>

#include <vector>


#include "fluids.hpp"

HexVolumetricMeshWithBoundary volumetricRing(vec3 center, vec3 normal, float radiusSmall, float radiusBig, float height, int radial_res, int vertical_res, int extrude_res);
HexVolumetricMeshWithBoundary volumetricCylinder(vec3 center, vec3 normal, float radius, float height, int radial_res, int vertical_res, int extrude_res);
HexVolumetricMeshWithBoundary primitiveHex(std::array<vec3, 8> vertices);
HexVolumetricMeshWithBoundary evenGrid(int n, int m, int l, float lenx, float leny, float lenz);



FluidSimulation heatFlowRigidBody(HexVolumetricMeshWithBoundary mesh, std::vector<float> temperatures, float thermalConductivity);

//CellMesh gridCellMesh(glm::ivec3 dims, vec3 len, const std::map<std::string, std::map<glm::ivec3, float>>& attributes);
CellMesh gridCellMesh(glm::ivec3 dims, vec3 len, const std::map<std::string, HOM(glm::ivec3, float)> &attributes);
CellMesh flatGridWithHeight(glm::ivec2 dims, vec2 len, HOM(glm::ivec2, float) h, const std::map<std::string, HOM(glm::ivec2, float)> &attributes);
CellMesh flatCustomGrid(glm::ivec2 dims, const HOM(glm::ivec2, vec3) &position_lower, const HOM(glm::ivec2, vec3) &position_higher, const std::map<std::string, HOM(glm::ivec2, float)> &attributes);

CellMesh flatCustomSubgrid(glm::ivec2 dims, const HOM(glm::ivec2, vec3) &position_lower, const HOM(glm::ivec2, vec3) &position_higher, const std::vector<glm::ivec2>& removed, const std::map<std::string, HOM(glm::ivec2, float)> &attributes);
CellMesh flatSubgrid(glm::ivec2 dims, vec3 len, const std::vector<glm::ivec2>& removed, const std::map<std::string, HOM(glm::ivec2, float)> &attributes);


CellMesh flatCylinder(glm::ivec2 dims, float r, float R, float h, const std::map<std::string, HOM(glm::ivec2, float)> &attributes);





class HeatFlow : public CellMesh {
	float thermalConductivity;
	float heatCapacity;
public:
	HeatFlow(const std::vector<Cell> &cells, float thermalConductivity, float heatCapacity) : CellMesh(cells), thermalConductivity(thermalConductivity), heatCapacity(heatCapacity) {}
	HeatFlow(const CellMesh &cells, float thermalConductivity, float heatCapacity) : CellMesh(cells), thermalConductivity(thermalConductivity), heatCapacity(heatCapacity) {}

	float heat(const Cell &c) const;
	float temperatureFluxPerFace(const Cell &c, FACE f) { return thermalConductivity * gradientOnFace(c, "T", f); }
	float thermalDiffusivity() const { return thermalConductivity/heatCapacity; }
	float heatLaplacian(const Cell &c) { return Laplacian(c, "T"); }
	float changeOverTime(const Cell &c);
	void update(float dt) {
		std::map<PolyGroupID, float> changes;
		for (auto & cell : cells)
			changes[cell.getID()] = changeOverTime(cell);
		for (auto & cell : cells)
			cell.addToAttribute("T",  changes[cell.getID()]*dt);
	}
};

VectorFieldR3Dynamic CouetteFlow(float h, float v0);
