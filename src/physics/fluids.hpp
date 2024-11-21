#pragma once
#include <utility>

#include "rigid.hpp"

enum CELL_ATTRIBUTE {
	VELOCITY="v",
	PRESSURE="p",
	TEMPERATURE="T",
	ANGULAR_VELOCITY="omega",
	ANGULAR_ACCELERATION="alpha",
	FORCE="F",
	TORQUE="tau",
	MASS="m",
	DENSITY="rho",
	HEAT_CAPACITY="k",
	EXTRA_FLOAT="extra1",
	EXTRA_VEC2="extra2",
	EXTRA_VEC3="extra3",
	EXTRA_VEC4="extra4"
};

int cellAttributeSize(CELL_ATTRIBUTE a);
int cellAttributeOffset(CELL_ATTRIBUTE a, std::vector<CELL_ATTRIBUTE> activeAttributes);
int attributeVectorSize(std::vector<CELL_ATTRIBUTE> activeAttributes);



class CellAttributeMatrix : public BigMatrix {
	std::vector<CELL_ATTRIBUTE> activeAttributes;

public:
	CellAttributeMatrix(int n, int m, std::vector<CELL_ATTRIBUTE> activeAttributes) : BigMatrix(n, attributeVectorSize(activeAttributes)), activeAttributes(std::move(activeAttributes)) {}
	CellAttributeMatrix(int n, int m, MATR$X coefs,  std::vector<CELL_ATTRIBUTE> activeAttributes) : BigMatrix(coefs), activeAttributes(std::move(activeAttributes)) {}
};


class CFDHexahedralMesh : public IndexedHexahedralVolumetricMesh {
	std::vector<CELL_ATTRIBUTE> cellAttributes;
	std::vector<CELL_ATTRIBUTE> cellUniforms;
	int parameterSize() const { return attributeVectorSize(cellAttributes); }
	MATR$X buffer;
	std::vector<CellAttributeMatrix> localDeformations;

public:
	CFDHexahedralMesh(const std::vector<Vertex> &vertices, const std::vector<ivec8>& hexahedra, const std::vector<int> &bdFacesInd);
	void updateBuffer(CELL_ATTRIBUTE, int i, float value);
	void updateBuffer(CELL_ATTRIBUTE, int i, vec2 value);
	void updateBuffer(CELL_ATTRIBUTE, int i, vec3 value);
	void updateBuffer(CELL_ATTRIBUTE, int i, vec4 value);
	void updateMatrix(int i, HOM(CellAttributeMatrix*, void) updateOperator);
	void updateMatrices (BIHOM(CellAttributeMatrix*, int, void) updateOperators);
};





class FluidSimulationFVM {
	CFDHexahedralMesh mesh;


	float t=0;
	HOM(float, HOM(HexahedronWithNbhd&, void)) __step;
	void _step(float t) { mesh.applyPerCell(__step(t)); }
public:
	FluidSimulationFVM(CFDHexahedralMesh mesh, HOM(float, HOM(HexahedronWithNbhd&, void)) __step) : mesh(std::move(mesh)), forces(std::move(forces)), __step(std::move(__step)) {}
	void step(float dt) {_step(dt); t+=dt;}



};
