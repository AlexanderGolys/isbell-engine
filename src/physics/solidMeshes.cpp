#include "solidMeshes.hpp"

#include <map>
#include <glm/detail/_vectorize.hpp>
#include <glm/detail/_vectorize.hpp>
#include <utility>

using namespace glm;


std::array<ivec3, 12> IndexedHexahedron::getTriangles() const {
	return {
		ivec3(faceDown[0], faceDown[1], faceDown[2]),
		ivec3(faceDown[0], faceDown[2], faceDown[3]),
		ivec3(faceUp[0], faceUp[1], faceUp[2]),
		ivec3(faceUp[0], faceUp[2], faceUp[3]),
		ivec3(faceDown[0], faceDown[1], faceUp[1]),
		ivec3(faceDown[0], faceUp[1], faceUp[0]),
		ivec3(faceDown[1], faceDown[2], faceUp[2]),
		ivec3(faceDown[1], faceUp[2], faceUp[1]),
		ivec3(faceDown[2], faceDown[3], faceUp[3]),
		ivec3(faceDown[2], faceUp[3], faceUp[2]),
		ivec3(faceDown[3], faceDown[0], faceUp[0]),
		ivec3(faceDown[3], faceUp[0], faceUp[3]),
	};
}


std::array<ivec2, 8> IndexedHexahedron::getEdges() const {
	return {
	ivec2(faceDown[0], faceDown[1]),
	ivec2(faceDown[1], faceDown[2]),
	ivec2(faceDown[2], faceDown[3]),
	ivec2(faceDown[3], faceDown[0]),
	ivec2(faceUp[0], faceUp[1]),
	ivec2(faceUp[1], faceUp[2]),
	ivec2(faceUp[2], faceUp[3]),
	ivec2(faceUp[3], faceUp[0]), };
}

float IndexedHexahedron::volume() const {
	return abs(dot(getFaceDown().getVertices()[0].getPosition() - getFaceUp().getVertices()[0].getPosition(), cross(getFaceDown().getVertices()[1].getPosition() - getFaceUp().getVertices()[0].getPosition(), getFaceDown().getVertices()[2].getPosition() - getFaceUp().getVertices()[0].getPosition())))/6.f;
}
float IndexedHexahedron::surfaceArea() const {
	return sum<float, std::array<float, 6>>(mapByVal<IndexedQuadrilateral, float, 6>(getFaces(), [](IndexedQuadrilateral q) { return q.surfaceArea(); }));
}


vec3 IndexedHexahedron::centerOfMass() const {
	return sum<vec3, std::array<vec3, 8>>(mapMove<Vertex, vec3>(getVertices(), [](Vertex v){return v.getPosition();}))/8.f;
}



vec3 IndexedHexahedron::centerOfFace(FACE i) const {
	return sum<vec3, std::array<vec3, 4>>(getFace(i).vertexPositions());
}

float IndexedHexahedron::radiusOfFaceCenter(FACE i) const {
	return norm(projectVectorToPlane(centerOfFace(i) - centerOfMass(), getFace(i).normal()));
}

mat3 IndexedHexahedron::inertiaTensor(float mass) const {
	auto perVert = [](vec3 pos) {
		return mat3( pos.y * pos.y + pos.z * pos.z, -pos.x * pos.y, -pos.x * pos.z,
						-pos.x * pos.y, pos.x * pos.x + pos.z * pos.z, -pos.y * pos.z,
						-pos.x * pos.z, -pos.y * pos.z, pos.x * pos.x + pos.y * pos.y ); };
	return sum<mat3, std::array<mat3, 8>>(mapMove<vec3, mat3>(getVerticesPos(), perVert)) / 8.f * mass;
}

mat3 IndexedHexahedron::inertiaTensor(std::array<float, 8> vertMasses) const {
	auto perVert = [](vec3 pos) {
		return mat3( pos.y * pos.y + pos.z * pos.z, -pos.x * pos.y, -pos.x * pos.z,
						-pos.x * pos.y, pos.x * pos.x + pos.z * pos.z, -pos.y * pos.z,
						-pos.x * pos.z, -pos.y * pos.z, pos.x * pos.x + pos.y * pos.y ); } ;
	auto vert = getVerticesPos();
	return sum<mat3, std::array<mat3, 8>>(arrayComprehension<mat3, 8, HOM(int, mat3)>([&perVert, &vertMasses, &vert](int i){return perVert(vert[i]*vertMasses[i]);}));
}

mat3 IndexedHexahedron::inertiaTensorAroundPoint(vec3 p) const {throw std::logic_error("Not v implemented"); }; // TODO


bool IndexedHexahedron::convex() const {}
IndexedHexahedron IndexedHexahedron::convexHull() const { throw std::logic_error("Not implemented"); }
bool IndexedHexahedron::contains(vec3 p) const { throw std::logic_error("Not implemented"); }
vec3 IndexedHexahedron::faceCenterNormalPoint(int i) const { throw std::logic_error("Not implemented"); }
vec3 IndexedHexahedron::interpolateFromVertices(std::array<float, 8> values) const { throw std::logic_error("Not implemented"); }
vec3 IndexedHexahedron::interpolateFromVertices(std::array<float, 8> values, vec3 p) const { throw std::logic_error("Not implemented"); }
vec3 IndexedHexahedron::interpolateFromEdges(std::array<float, 12> values) const { throw std::logic_error("Not implemented"); }
vec3 IndexedHexahedron::interpolateFromEdges(std::array<float, 12> values, vec3 p) const { throw std::logic_error("Not implemented"); }
vec3 IndexedHexahedron::interpolateFromFaces(std::array<float, 6> values) const { throw std::logic_error("Not implemented"); }
vec3 IndexedHexahedron::interpolateFromFaces(std::array<float, 6> values, vec3 p) const { throw std::logic_error("Not implemented"); }
function<SmoothParametricSurface(int)> IndexedHexahedron::pencilOfSections() const { throw std::logic_error("Not implemented"); }

HexahedronWithNbhd::HexahedronWithNbhd(const ivec8 &face, function<Vertex(int)> bufferManagerInterfacel, function<HexahedronWithNbhd *(int)> bufferManagerHexa, const std::map<FACE, int *> &neighbours)
:	IndexedHexahedron(face, std::move(bufferManagerInterfacel)),
	_bufferManagerHexa(std::move(bufferManagerHexa)),
	nbhd(neighbours) { }


std::vector<HexahedronWithNbhd *> HexahedronWithNbhd::getNeighbours() const {
	std::vector<HexahedronWithNbhd*> nbhd = {};
	for (auto &fc : {DOWN, UP, RIGHT, LEFT, FRONT, BACK})
		if (hasNeighbour(fc))
			nbhd.emplace_back(*getNeighbour(fc), _bufferManagerHexa, std::map<FACE, glm::ivec3*>());
	return nbhd;
}

IndexedHexahedralVolumetricMesh::IndexedHexahedralVolumetricMesh(std::vector<Vertex> vertices, const std::vector<ivec8>& hexahedra, std::set<int> bdFacesInd):
				vertices(std::move(vertices)),
				bdFacesInd(std::move(bdFacesInd))
{
	this->hexahedra = std::vector<HexahedronWithNbhd>();
	bufferManagerAccess = [&vertices](int i) { return vertices[i]; };
	for (auto &hex : hexahedra)
		this->hexahedra.push_back(initHexahedron(hex));

	faces = std::vector<IndexedQuadrilateral>();
	for (auto &hex : this->hexahedra)
		for (auto &face : hex.getFaces())
			faces.push_back(face);
}

HexahedronWithNbhd IndexedHexahedralVolumetricMesh::initHexahedron(const ivec8 &ind, const std::map<FACE, int *> &nbhd) const { return HexahedronWithNbhd(ind, bufferManagerAccess, bufferManagerHexa, nbhd); }

WeakSuperMesh IndexedHexahedralVolumetricMesh::triangulateBoundary(PolyGroupID id) const {
	vector<ivec3> tris = {};
	for (auto &face : bdFacesInd) {
		tris.push_back(hexahedra[face].getTriangles().at(0));
		tris.push_back(hexahedra[face].getTriangles().at(1));
	}
	return WeakSuperMesh(vertices, tris, id);
}

void IndexedHexahedralVolumetricMesh::addCell(const HexahedronWithNbhd &h, const std::set<FACE> &nbhdFaces) { hexahedra.push_back(h); for (auto &fc : nbhdFaces) addNeighbour(hexahedra.size()-1, *h.getNeighbour(fc), fc); }
void IndexedHexahedralVolumetricMesh::addFace(const IndexedQuadrilateral &q, bool bd) { faces.push_back(q); if (bd) bdFacesInd.insert(faces.size()-1); }


void IndexedHexahedralVolumetricMesh::subdivideCell(int i) {
	vec3 p1 = hexahedra[i].faceCenterNormalPoint(DOWN);
	vec3 p2 = hexahedra[i].faceCenterNormalPoint(UP);
	vec3 p3 = hexahedra[i].faceCenterNormalPoint(RIGHT);
	vec3 p4 = hexahedra[i].faceCenterNormalPoint(LEFT);
	vec3 p5 = hexahedra[i].faceCenterNormalPoint(FRONT);
	vec3 p6 = hexahedra[i].faceCenterNormalPoint(BACK);
	vec3 c = hexahedra[i].centerOfMass();

//	addPoint(p1);
//	auto lowRightBack  = IndexedHexahedron()
}


void IndexedHexahedralVolumetricMesh::subdivide() { throw std::logic_error("Not implemented"); }
