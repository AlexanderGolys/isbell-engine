#pragma once
#include "rigid.hpp"


class FluidSimulation {
protected:
	HexVolumetricMeshWithBoundary mesh;
	std::vector<BigVector> cellAttributes;
	BIHOM(int, float, BigVector) cellAttributeChange;
	BIHOM(BufferedVertex&, std::vector<float>&, void) updateBdWithCellIndex;
public:
	FluidSimulation(const std::vector<ivec6> &hexahedra, const std::vector<glm::ivec4> &faces,
		const std::vector<Vertex> &vertices, const std::vector<ivec2> &faceNbhd,
		const std::vector<int>& bdFaces,
		const MATR$X& initialAttributes, BIHOM(int, float, BigVector) cellAttributeChange,
		BIHOM(BufferedVertex&, std::vector<float>&, void) updateBdWithCellIndex);

	FluidSimulation(HexVolumetricMeshWithBoundary mesh,
		const MATR$X& initialAttributes, BIHOM(int, float, BigVector) cellAttributeChange,
		BIHOM(BufferedVertex&, std::vector<float>&, void) updateBdWithCellIndex);

	BigVector attributesWithNbhrs(int i) const;
	BigVector meanVertexAttributes(int i) const;
	void updateAttributes (float dt);

	void updateBdMesh();
	void update (float dt) { updateAttributes(dt); updateBdMesh(); }
	std::shared_ptr<WeakSuperMesh> getBoundaryMesh() const { return mesh.boundary; }
	BIHOM(float, float, void) updater() { return [this](float t, float dt) { this->update(dt); }; }
};


class FluidSimulationFVM : public FluidSimulation {
	std::vector<BigMatrix> attributeMatrices;
public:
	FluidSimulationFVM(const std::vector<ivec6> &hexahedra, const std::vector<glm::ivec4> &faces,
		const std::vector<Vertex> &vertices, const std::vector<ivec2> &faceNbhd,
		const std::vector<int>& bdFaces,
		const MATR$X& initialAttributes, std::vector<BigMatrix> attributeMatrices,
		BIHOM(BufferedVertex&, std::vector<float>&, void) updateBdWithCellIndex);

	FluidSimulationFVM(HexVolumetricMeshWithBoundary mesh,
		const MATR$X& initialAttributes, std::vector<BigMatrix> attributeMatrices,
		BIHOM(BufferedVertex&, std::vector<float>&, void) updateBdWithCellIndex);

	FluidSimulationFVM(HexVolumetricMeshWithBoundary mesh,
		HOM(HexahedronWithNbhd&, BigVector) initAttr, HOM(HexahedronWithNbhd&, BigMatrix) initMatrices,
		BIHOM(BufferedVertex&, std::vector<float>&, void) updateBdWithCellIndex);

};

class IncompressibleCell {
public:
	bool isValid = true;
	vec2 upleft, upright, downleft, downright;
	bool upBd = false, downBd = false, leftBd = false, rightBd = false;
	float density;
	vec2 velocity;
	float pressure;
	vec2 pressureGradient;
	float temperature;
	float viscosity;
	IncompressibleCell(vec2 upleft, vec2 upright, vec2 downleft, vec2 downright, float density, vec2 velocity, float pressure, float temperature, float viscosity, bool upBd=false, bool downBd=false, bool leftBd=false, bool rightBd=false) : upBd(upBd), downBd(downBd), leftBd(leftBd), rightBd(rightBd),
        density(density), velocity(velocity), pressure(pressure), temperature(temperature), viscosity(viscosity) , upleft(upleft), upright(upright), downleft(downleft), downright(downright) {}

	IncompressibleCell() : isValid(false),  density(0), velocity(vec2(0)), pressure(0), temperature(0), viscosity(0) {}
	vector<float> toVector() const { return {density, velocity.x, velocity.y, pressure, temperature, viscosity}; }
	float rho() const { return density; }
	vec2 v() const { return velocity; }
	float p() const { return pressure; }
	float T() const { return temperature; }
	float mu() const { return viscosity; }
	vec2 corner(int i) const;
	vec2 center() const;
	bool valid() const { return isValid; }

	void updateDensity(float newDensity) { density = newDensity; }
	void updateVelocity(vec2 newVelocity) { velocity = newVelocity; }
	void updatePressure(float newPressure) { pressure = newPressure; }
	void updateTemperature(float newTemperature) { temperature = newTemperature; }
	void updateViscosity(float newViscosity) { viscosity = newViscosity; }
	void updateCell(const IncompressibleCell &newCell) { *this = newCell; }
	void updatePressureGradient(vec2 newPressureGradient) { pressureGradient = newPressureGradient; }
};


class IncompressibleNewtonianFluid2D {
	vector<IncompressibleCell> cells;
	ivec2 gridSize;
	BigMatrix M;


public:
	explicit IncompressibleNewtonianFluid2D(const vector<vector<IncompressibleCell>> &cells) :  cells(flatten(cells)), gridSize(cells.size(), cells[0].size()), M(6*gridSize.x*gridSize.y, 6*gridSize.x*gridSize.y) {}
	void updateMatrix (float dt);
	vec69 velocities_x_vector() const;
	vec69 velocities_y_vector() const;
	void updatePressuresGradients(const vec69 &pressures_x, const vec69 &pressures_y) {for (int i = 0; i < cells.size(); i++) cells[i].updatePressureGradient(vec2(pressures_x[i], pressures_y[i]));}
	void updateVelocities(const vec69 &velocities_x, const vec69 &velocities_y) {for (int i = 0; i < cells.size(); i++) cells[i].updateVelocity(vec2(velocities_x[i], velocities_y[i]));}
	void updateDensities(const vec69 &densities) {for (int i = 0; i < cells.size(); i++) cells[i].updateDensity(densities[i]);}

};
