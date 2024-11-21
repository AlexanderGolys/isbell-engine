#include "fluids.hpp"



HexahedralMesh::HexahedralMesh(const std::vector<Vertex> &vertices, const std::vector<ivec8> &hexahedra, const std::vector<int> &bdFacesInd, std::vector<VolumeElementAttributes> attributes):
	IndexedHexahedralVolumetricMesh(vertices, hexahedra, bdFacesInd),
	attributes(std::move(attributes))
{}
