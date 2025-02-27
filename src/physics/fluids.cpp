#include "fluids.hpp"

//
//FluidSimulation::FluidSimulation(const std::vector<ivec6> &hexahedra, const std::vector<glm::ivec4> &faces,
//		const std::vector<Vertex> &vertices, const std::vector<glm::ivec2> &faceNbhd,
//		const std::vector<int>& bdFaces,
//		const MATR$X& initialAttributes,  BIHOM(int, float, FloatVector) cellAttributeChange,
//		BIHOM(BufferedVertex&, std::vector<float>&, void) updateBdWithCellIndex) :
//			mesh(hexahedra, faces, vertices, bdFaces),
//			 cellAttributeChange(std::move(cellAttributeChange) ),
//			 updateBdWithCellIndex(std::move(updateBdWithCellIndex)) {
//
//	cellAttributes = {};
//	for (const auto & initialAttribute : initialAttributes)
//		cellAttributes.emplace_back(initialAttribute);
//}
//
//FluidSimulation::FluidSimulation(HexVolumetricMeshWithBoundary mesh,
//		const MATR$X& initialAttributes, BIHOM(int, float, FloatVector) cellAttributeChange,
//		BIHOM(BufferedVertex&, std::vector<float>&, void) updateBdWithCellIndex) :
//			mesh(std::move(mesh)),
//			cellAttributeChange(std::move(cellAttributeChange)),
//			updateBdWithCellIndex(std::move(updateBdWithCellIndex)) {
//
//	cellAttributes = {};
//	for (const auto & initialAttribute : initialAttributes)
//		cellAttributes.emplace_back(initialAttribute);
//}
//
//FloatVector FluidSimulation::attributesWithNbhrs(int i) const {
//	FloatVector result = cellAttributes[i];
//	int len = result.size();
//	auto neighbours = mesh.mesh.getCell(i).getNeighboursIndices();
//	for (int nbhr : neighbours) {
//		FloatVector attr = nbhr == -1 ? FloatVector(len) : cellAttributes[nbhr];
//		result.append(attr);
//	}
//	return result;
//}
//
//
//FloatVector FluidSimulation::meanVertexAttributes(int i) const {
//	auto cells = mesh.bdVerticesCells.at(i);
//	if (cells.empty()) throw std::logic_error("No cells for vertex");
//	FloatVector result = attributesWithNbhrs(cells[0])/cells.size();
//	for (int j = 1; j < cells.size(); j++)
//		result += attributesWithNbhrs(cells[j])/cells.size();
//	return result;
//}
//
//void FluidSimulation::updateAttributes(float dt) {
//	for (int i = 0; i < cellAttributes.size(); i++)
//		cellAttributes[i] += cellAttributeChange(i, dt);
//}
//
//
//
//void FluidSimulation::updateBdMesh() {
//	mesh.boundary->deformPerVertex(mesh.boundaryPolygroup, [this](BufferedVertex &v) {
//		vec69 attr = this->meanVertexAttributes(v.getIndex()).getVec();
//		this->updateBdWithCellIndex(v, attr);
//	});
//}
void IncompressibleNewtonianFluid2D::updatePressuresGradients(const vector<float> &pressures_x, const vector<float> &pressures_y) {for (int i = 0; i < cells.size(); i++) cells[i].updatePressureGradient(vec2(pressures_x[i], pressures_y[i]));}
