#pragma once
#include <utility>

#include <array>
#include "src/common/specific.hpp"

enum FACE {DOWN = 0, UP = 1, RIGHT = 2, LEFT = 3, FRONT = 4, BACK = 5};

class Quadrilateral {
	Vertex p0, p1, p2, p3;
public:
	Quadrilateral(vec3 p0, vec3 p1, vec3 p2, vec3 p3) : p0(p0), p1(p1), p2(p2), p3(p3) {}
	std::array<Vertex, 4> getVertices() const;
	static std::array<glm::ivec3, 2> getTriangles() { return {glm::ivec3(0, 1, 2), glm::ivec3(0, 2, 3)}; }
};




inline mat3 toMat3(const std::array<vec3, 3> &m) { return mat3(m[0], m[1], m[2]); }


typedef std::array<vec3, 3> arr33;

class IndexedQuadrilateral {
	glm::ivec4 indices;
	HOM(int, Vertex) bufferManagerInterface;
public:
	IndexedQuadrilateral(glm::ivec4 indices, HOM(int, Vertex) bufferManagerInterface) : indices(indices), bufferManagerInterface(std::move(bufferManagerInterface)) {}
	std::array<Vertex, 4> getVertices() const { return {bufferManagerInterface(indices[0]), bufferManagerInterface(indices[1]), bufferManagerInterface(indices[2]), bufferManagerInterface(indices[3])}; }
	std::array<glm::ivec3, 2> getTriangles() const { return {glm::ivec3(indices[0], indices[1], indices[2]), glm::ivec3(indices[0], indices[2], indices[3])}; }
	std::array<std::array<vec3, 3>, 2> getTriangleVertexPositions() const { return {bufferManagerInterface(indices[0]).getPosition(), bufferManagerInterface(indices[1]).getPosition(), bufferManagerInterface(indices[2]).getPosition()}; }
	float surfaceArea() const { std::array<arr33, 2> v = getTriangleVertexPositions(); return abs(det(toMat3(v[0]))) +  abs(det(toMat3(v[1]))); }
	std::array<vec3, 4> vertexPositions() const {return {bufferManagerInterface(indices[0]).getPosition(), bufferManagerInterface(indices[1]).getPosition(), bufferManagerInterface(indices[2]).getPosition(), bufferManagerInterface(indices[3]).getPosition()}; }
	vec3 vertexPosition(int i) const { return bufferManagerInterface(indices[i]).getPosition(); }
	vec3 normal() const { return normalize(cross(vertexPosition(1) - vertexPosition(0), vertexPosition(2) - vertexPosition(0))); }
};




class IndexedHexahedron {
		glm::ivec4 faceDown, faceUp;
		HOM(int, Vertex) bufferManagerInterface;
public:
		IndexedHexahedron(glm::ivec4 faceDown, glm::ivec4 faceUp, HOM(int, Vertex) bufferManagerInterface) : faceDown(faceDown), faceUp(faceUp), bufferManagerInterface(std::move(bufferManagerInterface)) {}
	IndexedHexahedron(const ivec8 &face, HOM(int, Vertex) bufferManagerInterface) : faceDown(face.xyzw()), faceUp(face.stuv()), bufferManagerInterface(std::move(bufferManagerInterface)) {}
	IndexedQuadrilateral getFaceDown() const { return IndexedQuadrilateral(faceDown, bufferManagerInterface); }
		IndexedQuadrilateral getFaceUp()   const { return IndexedQuadrilateral(faceUp, bufferManagerInterface); }
		IndexedQuadrilateral getFaceRight()const  { return IndexedQuadrilateral(glm::ivec4(faceDown[1], faceDown[2], faceUp[2], faceUp[1]), bufferManagerInterface); }
		IndexedQuadrilateral getFaceLeft() const  { return IndexedQuadrilateral(glm::ivec4(faceDown[0], faceDown[3], faceUp[3], faceUp[0]), bufferManagerInterface); }
		IndexedQuadrilateral getFaceFront()const  { return IndexedQuadrilateral(glm::ivec4(faceDown[0], faceDown[1], faceUp[1], faceUp[0]), bufferManagerInterface); }
		IndexedQuadrilateral getFaceBack() const  { return IndexedQuadrilateral(glm::ivec4(faceDown[3], faceDown[2], faceUp[2], faceUp[3]), bufferManagerInterface); }
		std::array<Vertex, 8> getVertices() const { return {bufferManagerInterface(faceDown[0]), bufferManagerInterface(faceDown[1]), bufferManagerInterface(faceDown[2]), bufferManagerInterface(faceDown[3]), bufferManagerInterface(faceUp[0]), bufferManagerInterface(faceUp[1]), bufferManagerInterface(faceUp[2]), bufferManagerInterface(faceUp[3])}; }
		std::array<vec3, 8> getVerticesPos() const { auto ; return mapMove< Vertex, vec3, 8>(getVertices(), [](const Vertex& v) { return v.getPosition(); }); }
		std::array<IndexedQuadrilateral, 6> getFaces() const { return {getFaceDown(), getFaceUp(), getFaceRight(), getFaceLeft(), getFaceFront(), getFaceBack()}; }
		IndexedQuadrilateral getFace(FACE fc) const { return fc == DOWN ? getFaceDown() : fc == UP ? getFaceUp() : fc == RIGHT ? getFaceRight() : fc == LEFT ? getFaceLeft() : fc == FRONT ? getFaceFront() : getFaceBack(); }
		std::array<glm::ivec3, 12> getTriangles() const;
		std::array<glm::ivec2, 8> getEdges() const;
		ivec8 getIndices() const { return ivec8(faceDown, faceUp); }

		float volume() const;
		float surfaceArea() const;
		vec3 centerOfFace(FACE i) const;
		vec3 centerOfMass() const;
		float radiusOfFaceCenter(FACE i) const;
		mat3 inertiaTensor(float mass) const;
		mat3 inertiaTensor(std::array<float, 8> vertMasses) const;
		mat3 inertiaTensorAroundPoint(vec3 p) const;
		bool convex() const;
		IndexedHexahedron convexHull() const;
		bool contains(vec3 p) const;
		vec3 faceCenterNormalPoint(int i) const;
		vec3 interpolateFromVertices(std::array<float, 8> values) const;
		vec3 interpolateFromVertices(std::array<float, 8> values, vec3 p) const;
		vec3 interpolateFromEdges(std::array<float, 12> values) const;
		vec3 interpolateFromEdges(std::array<float, 12> values, vec3 p) const;
		vec3 interpolateFromFaces(std::array<float, 6> values) const;
		vec3 interpolateFromFaces(std::array<float, 6> values, vec3 p) const;
		HOM(int, SmoothParametricSurface) pencilOfSections() const;
};


class HexahedronWithNbhd : public IndexedHexahedron {
	HOM(int, HexahedronWithNbhd*) _bufferManagerHexa;
	std::map<FACE, int*> nbhd;
public:
	HexahedronWithNbhd(const ivec8 &face, HOM(int, Vertex) bufferManagerInterfacel, HOM(int, HexahedronWithNbhd*) bufferManagerHexa, const std::map<FACE, int*> &neighbours);

	void setNeighbour(FACE fc, int *neighbour) { nbhd[fc] = neighbour; }
	bool hasNeighbour(FACE fc) const { return nbhd.contains(fc) && nbhd.at(fc) != nullptr; }
	int* getNeighbour(FACE fc) const { return hasNeighbour(fc) ? nbhd.at(fc) : nullptr; }
	std::vector<HexahedronWithNbhd*> getNeighbours() const;
};


class IndexedHexahedralVolumetricMesh {
	std::vector<Vertex> vertices;
	std::vector<HexahedronWithNbhd> hexahedra;
	std::vector<IndexedQuadrilateral> faces;
	std::set<int> bdFacesInd;
	HOM(int, Vertex) bufferManagerAccess;
	HOM(int, HexahedronWithNbhd*) bufferManagerHexa;
public:
	IndexedHexahedralVolumetricMesh(std::vector<Vertex> vertices, const std::vector<ivec8>& hexahedra, std::set<int> bdFacesInd);
	IndexedQuadrilateral initFace(glm::ivec4 indices) const { return IndexedQuadrilateral(indices, bufferManagerAccess); }
	HexahedronWithNbhd initHexahedron(const ivec8 &ind) const { return HexahedronWithNbhd(ind, bufferManagerAccess, bufferManagerHexa, std::map<FACE, int*>()); }
	HexahedronWithNbhd initHexahedron(const ivec8 &ind, const std::map<FACE, int*>& nbhd) const;

	void addNeighbour(int nbhdInd, int target, FACE side) { hexahedra[target].setNeighbour(side, &nbhdInd); }
	WeakSuperMesh triangulateBoundary(std::variant<int, std::string> id) const;

	void addVertex(const Vertex& v) { vertices.push_back(v); }
	void addCell(const HexahedronWithNbhd& h, const std::set<FACE>& nbhdFaces={});
	void addFace(const IndexedQuadrilateral &q, bool bd);
	void setAsBoundary(int i) { bdFacesInd.insert(i); }
	int size() const { return hexahedra.size(); }

	void applyPerCell(HOM(HexahedronWithNbhd&, void) f) { for (auto &hex : hexahedra) f(hex); }
	void subdivide();
	void subdivideCell(int i);
};
