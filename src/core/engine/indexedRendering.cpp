// ReSharper disable CppMemberFunctionMayBeConst
#include "indexedRendering.hpp"

#include <stdio.h>
#include <vector>
#include <fstream>
#include <ranges>
#include <sstream>

#include "logging.hpp"
#include "SDFObjects.hpp"

using namespace glm;
using namespace BufferNames;

Vertex::Vertex(vec3 position, vec2 uv, vec3 normal, vec4 color, const std::map<string, vec4>& extraData): position(position), normal(normal), uv(uv), color(color), extraData(extraData) {}

vec4 Vertex::getExtraData(const string& name) const {
	THROW_IF(not extraData.contains(name), IndexError, "No extra data with name '" + name + "' found in vertex.");
	return extraData.at(name);
}

void Vertex::setExtraData(const string& name, vec4 value) {
	extraData[name] = value;
}

bool Vertex::hasExtraData(const string& name) const {
	return extraData.contains(name);
}

Vertex center(const Vertex& v1, const Vertex& v2) {
	std::map<string, vec4> extraData;
	for (const auto& [key, val] : v1.extraData)
		extraData[key] = 0.5f * (val + v2.getExtraData(key));
	return Vertex(
		0.5f * (v1.get_position() + v2.get_position()),
		0.5f * (v1.get_uv() + v2.get_uv()),
		normalise(0.5f * (v1.get_normal() + v2.get_normal())),
		0.5f * (v1.get_color() + v2.get_color()),
		extraData
	);
}

Vertex barycenter(const Vertex& v1, const Vertex& v2, const Vertex& v3) {
		std::map<string, vec4> extraData;
	for (const auto& [key, val] : v1.extraData)
		extraData[key] = (val + v2.getExtraData(key) + v3.getExtraData(key))/3.f;
	return Vertex(
		(v1.get_position() + v2.get_position() + v3.get_position()) / 3.f,
		(v1.get_uv() + v2.get_uv() + v3.get_uv()) / 3.f,
		normalise(v1.get_normal() + v2.get_normal() + v3.get_normal()),
		(v1.get_color() + v2.get_color() + v3.get_color()) / 3.f,
		extraData
	);
}

unique_ptr<vector<vec4>>& BufferManager::extraBufferPtrAt(int slot)  {
	THROW_IF(slot >= extraBufferNames.size(), ValueError, "No extra buffer at slot " + to_string(slot) + ".");
	switch (slot) {
	case 0: return extra0;
	case 1: return extra1;
	case 2: return extra2;
	case 3: return extra3;
	case 4: return extra4;
	}
	THROW(IndexError, "No slot " + to_string(slot) + " for extra buffers available.");
}

vec4 BufferManager::extraElementAt(int index, int slot) const {
	if (slot >= extraBufferNames.size())
		THROW(ValueError, "No extra buffer at slot " + to_string(slot) + ".");
	switch (slot) {
	case 0: return extra0->at(index);
	case 1: return extra1->at(index);
	case 2: return extra2->at(index);
	case 3: return extra3->at(index);
	case 4: return extra4->at(index);
	}
	THROW(IndexError, "No slot " + to_string(slot) + " for extra buffers available.");}

unique_ptr<vector<vec4>>& BufferManager::extraBufferPtrByName(const string& name)  {
	return extraBufferPtrAt(indexOfExtraBuffer(name));
}

void BufferManager::setElementAt(int index, float value, int slot, int component) {
	switch (slot) {
	case 0: (*extra0)[index][component] = value; return;
	case 1: (*extra1)[index][component] = value; return;
	case 2: (*extra2)[index][component] = value; return;
	case 3: (*extra3)[index][component] = value; return;
	case 4: (*extra4)[index][component] = value; return;
	}
	THROW(IndexError, "No slot " + to_string(slot) + " for extra buffers available.");
}

void BufferManager::setElementAt(int index, vec4 value, int slot) {
	switch (slot) {
	case 0: (*extra0)[index] = value; return;
	case 1: (*extra1)[index] = value; return;
	case 2: (*extra2)[index] = value; return;
	case 3: (*extra3)[index] = value; return;
	case 4: (*extra4)[index] = value; return;
	}
	THROW(IndexError, "No slot " + to_string(slot) + " for extra buffers available.");
}

BufferManager::BufferManager(const vector<string>& extra_names)
: extraBufferNames(extra_names) {
	THROW_IF(extra_names.size() > 5, IllegalArgumentError, "BufferManager supports up to 5 extra attribute buffers.");
	THROW_IF(not disjoint<string>(extra_names, RESERVED_BUF), IllegalVariantError, "Extra buffer names cannot collide with reserved buffer names.");

	extra0 = extra_names.size() > 0 ? make_unique<vector<vec4>>() : nullptr;
	extra1 = extra_names.size() > 1 ? make_unique<vector<vec4>>() : nullptr;
	extra2 = extra_names.size() > 2 ? make_unique<vector<vec4>>() : nullptr;
	extra3 = extra_names.size() > 3 ? make_unique<vector<vec4>>() : nullptr;
	extra4 = extra_names.size() > 4 ? make_unique<vector<vec4>>() : nullptr;

	stds = make_unique<Stds>();
	indices = make_unique<vector<ivec3>>();
}

BufferManager::BufferManager(const string& extra_name): BufferManager(vector<string>{extra_name}) {}

BufferManager::BufferManager(const string& extra_name0, const string& extra_name1): BufferManager(vector<string>{extra_name0, extra_name1}) {}

BufferManager::BufferManager(const string& extra_name0, const string& extra_name1, const string& extra_name2): BufferManager(vector<string>{extra_name0, extra_name1, extra_name2}) {}

BufferManager::BufferManager(const string& extra_name0, const string& extra_name1, const string& extra_name2, const string& extra_name3): BufferManager(vector<string>{extra_name0, extra_name1, extra_name2, extra_name3}) {}

BufferManager::BufferManager(const string& extra_name0, const string& extra_name1, const string& extra_name2, const string& extra_name3, const string& extra_name4): BufferManager(vector<string>{extra_name0, extra_name1, extra_name2, extra_name3, extra_name4}) {}

vs_dim BufferManager::getAttributeDimension(const string& attributeName) const {
	THROW_IF(not isActive(attributeName), IllegalVariantError, "Attribute " + attributeName + " is not active in this buffer manager.");
	if (attributeName == POS_BUF or attributeName == NORM_BUF or attributeName == IND_BUF)
		return 3;
	if (attributeName == UV_BUF)
		return 2;
	return 4;
}


array_len BufferManager::attributeDataLength() const {
	return stds->positions.size();
}

array_len BufferManager::indexDataLength() const {
	return indices->size();
}

byte_size BufferManager::attributeDataSize(const string& attributeName) const {
	THROW_IF(attributeName == IND_BUF, IllegalVariantError, "Index buffer size should be queried using indexDataSize().");
	return attributeDataLength() * getAttributeDimension(attributeName) * sizeof(GL_FLOAT);
}

byte_size BufferManager::indexDataSize() const {
	return indexDataLength() * 3 * sizeof(GL_UNSIGNED_INT);
}

byte_size BufferManager::totalAttributeSize() const {
	return attributeDataLength() * (sizeof(vec3)*2 + sizeof(vec2) + sizeof(vec4)*(1+extraBufferNames.size()));
}

raw_data_ptr BufferManager::attributeDataPtr(const string& attributeName) const {
	if (attributeName == POS_BUF)
		return stds->positions.data();
	if (attributeName == NORM_BUF)
		return stds->normals.data();
	if (attributeName == UV_BUF)
		return stds->uvs.data();
	if (attributeName == COL_BUF)
		return stds->colors.data();
	if (contains<string>(extraBufferNames, attributeName)) {
		if (attributeName == extraBufferNames[0])
			return extra0->data();
		if (attributeName == extraBufferNames[1])
			return extra1->data();
		if (attributeName == extraBufferNames[2])
			return extra2->data();
		if (attributeName == extraBufferNames[3])
			return extra3->data();
		if (attributeName == extraBufferNames[4])
			return extra4->data();
	}
	THROW(IllegalVariantError, "Attribute " + attributeName + " is not active in this buffer manager.");
}
raw_data_ptr BufferManager::indexDataPtr() const {
	return indices->data();
}

bool BufferManager::isActive(const string& attributeName) const {
	return contains<string>(extraBufferNames, attributeName) or contains<string>(RESERVED_BUF, attributeName);
}

bool BufferManager::isActive(int i) const {
	return i < extraBufferNames.size();
}

string BufferManager::getExtraBufferName(int slot) const {
	return extraBufferNames[slot];
}

vector<string> BufferManager::getExtraBufferNames() const {
	return extraBufferNames;
}

vs_dim BufferManager::numberOfExtraBuffers() const {
	return static_cast<vs_dim>(extraBufferNames.size());
}

int BufferManager::indexOfExtraBuffer(const string& name) const {
	for (int i = 0; i < extraBufferNames.size(); i++)
		if (extraBufferNames[i] == name)
			return i;
	THROW(IndexError, "No extra buffer named " + name + " found.");
}


int BufferManager::addTriangleVertexIndices(ivec3 ind, int shift) const {
	indices->push_back(ind + ivec3(shift));
	return indices->size() - 1;
}


int BufferManager::addFullVertexData(const Vertex& v) const {
	stds->positions.emplace_back(v.get_position());
	stds->normals.emplace_back(v.get_normal());
	stds->uvs.emplace_back(v.get_uv());
	stds->colors.push_back(v.get_color());
	if (extraBufferNames.size() > 0)
		if (v.hasExtraData(extraBufferNames[0]))
			extra0->push_back(v.getExtraData(extraBufferNames[0]));
		else
			extra0->emplace_back(0, 0, 0, 0);
	if (extraBufferNames.size() > 1)
		if (v.hasExtraData(extraBufferNames[1]))
			extra1->push_back(v.getExtraData(extraBufferNames[1]));
		else
			extra1->emplace_back(0, 0, 0, 0);
	if (extraBufferNames.size() > 2)
		if (v.hasExtraData(extraBufferNames[2]))
			extra2->push_back(v.getExtraData(extraBufferNames[2]));
		else
			extra2->emplace_back(0, 0, 0, 0);
	if (extraBufferNames.size() > 3)
		if (v.hasExtraData(extraBufferNames[3]))
			extra3->push_back(v.getExtraData(extraBufferNames[3]));
		else
			extra3->emplace_back(0, 0, 0, 0);
	if (extraBufferNames.size() > 4)
		if (v.hasExtraData(extraBufferNames[4]))
			extra4->push_back(v.getExtraData(extraBufferNames[4]));
		else
			extra4->emplace_back(0, 0, 0, 0);
	return stds->positions.size() - 1;
}


BufferedVertex::BufferedVertex(BufferManager& bufferBoss, int index)
: bufferBoss(bufferBoss), index(index) {}

BufferedVertex::BufferedVertex(const BufferedVertex& other)
: bufferBoss(other.bufferBoss), index(other.index) {}

BufferedVertex::BufferedVertex(BufferedVertex&& other) noexcept
: bufferBoss(other.bufferBoss), index(other.index) {}

BufferedVertex::BufferedVertex(BufferManager& bufferBoss, const Vertex& v)
: bufferBoss(bufferBoss) {
	index = this->bufferBoss.addFullVertexData(v);
}

BufferedVertex& BufferedVertex::operator=(const BufferedVertex& other) {
	if (this == &other)
		return *this;
	bufferBoss = other.bufferBoss;
	index = other.index;
	return *this;
}

BufferedVertex& BufferedVertex::operator=(BufferedVertex&& other) noexcept {
	if (this == &other)
		return *this;
	bufferBoss = other.bufferBoss;
	index = other.index;
	return *this;
}

int BufferedVertex::getIndex() const {
	return index;
}

vec3 BufferedVertex::getPosition() const {
	return bufferBoss.getPosition(index);
}

vec3 BufferedVertex::getNormal() const {
	return bufferBoss.getNormal(index);
}

vec2 BufferedVertex::getUV() const {
	return bufferBoss.getUV(index);
}

vec4 BufferedVertex::getColor() const {
	return bufferBoss.getColor(index);
}

vec4 BufferedVertex::getExtra(const string& slot) const {
	return bufferBoss.getExtra(index, slot);
}

Vertex BufferedVertex::getVertex() const {
	dict(string, vec4) extraData;
	for (const string& extraName : bufferBoss.getExtraBufferNames()) {
		extraData[extraName] = getExtra(extraName);
	}

	return Vertex(getPosition(), getUV(), getNormal(), getColor(), extraData);
}

void BufferedVertex::setPosition(vec3 value) const {
	bufferBoss.setPosition(index, value);
}

void BufferedVertex::setNormal(vec3 value) const {
	bufferBoss.setNormal(index, value);
}

void BufferedVertex::setUV(vec2 value) const {
	bufferBoss.setUV(index, value);
}

void BufferedVertex::setUV(float value, int i) const {
	bufferBoss.setUV(index, value, i);
}

void BufferedVertex::setColor(vec4 value) const {
	bufferBoss.setColor(index, value);
}

void BufferedVertex::setColor(float value, int i) const {
	bufferBoss.setColor(index, value, i);
}

void BufferedVertex::setExtra(vec4 value, const string& name) {
	bufferBoss.setExtra(index, value, name);
}

void BufferedVertex::setExtra(vec3 value, const string& name) {
	bufferBoss.setExtra(index, value.x, name, 0);
	bufferBoss.setExtra(index, value.y, name, 1);
	bufferBoss.setExtra(index, value.z, name, 2);
}

void BufferedVertex::setExtra(float value, const string& name, int component) {
	bufferBoss.setExtra(index, value, name, component);
}

void BufferedVertex::applyFunction(const SpaceEndomorphism& f) {
	vec3 p = getPosition();
	vec3 n = getNormal();
	setPosition(f(p));
	setNormal(normalise(f.df(p) * n));
}

void BufferedVertex::setVertex(const Vertex& v) {
	setPosition(v.get_position());
	setUV(v.get_uv());
	setNormal(v.get_normal());
	setColor(v.get_color());
	auto extras = bufferBoss.getExtraBufferNames();
	for (int i = 0; i < extras.size(); i++) {
		string extraName = extras[i];
		if (v.hasExtraData(extraName))
			setExtra(v.getExtraData(extraName), extraName);
	}
}

IndexedTriangle::IndexedTriangle(const IndexedTriangle& other)
: index(other.index), bufferBoss(other.bufferBoss) {}

IndexedTriangle::IndexedTriangle(IndexedTriangle&& other) noexcept
: index(other.index), bufferBoss(other.bufferBoss) {}

IndexedTriangle::IndexedTriangle(BufferManager& bufferBoss, ivec3 index, int shift)
: index(bufferBoss.addTriangleVertexIndices(index, shift)), bufferBoss(bufferBoss) {}

IndexedTriangle& IndexedTriangle::operator=(const IndexedTriangle& other) {
	if (this == &other)
		return *this;
	index = other.index;
	bufferBoss = other.bufferBoss;
	return *this;
}

IndexedTriangle& IndexedTriangle::operator=(IndexedTriangle&& other) noexcept {
	if (this == &other)
		return *this;
	index = other.index;
	bufferBoss = other.bufferBoss;
	return *this;
}

ivec3 IndexedTriangle::getVertexIndices() const {
	return bufferBoss.getFaceIndices(index);
}

Vertex IndexedTriangle::getVertex(int i) const {
	return bufferBoss.getVertex(getVertexIndices()[i]);
}


mat3 IndexedTriangle::orthonormalFrame() const {
	vec3 p0 = getVertex(0).get_position();
	vec3 p1 = getVertex(1).get_position();
	vec3 p2 = getVertex(2).get_position();
	return GramSchmidtProcess(mat3(p0 - p2, p1 - p2, cross(p1 - p2, p0 - p2)));
}

vec3 IndexedTriangle::fromPlanar(vec2 v) const {
	mat3 frame = orthonormalFrame();
	return frame[0] * v.x + frame[1] * v.y + getVertex(2).get_position();
}

vec2 IndexedTriangle::toPlanar(vec3 v) const {
	mat3 frame = orthonormalFrame();
	return vec2(dot(frame[0], v - getVertex(2).get_position()), dot(frame[1], v - getVertex(2).get_position()));
}


vec3 IndexedTriangle::fromBars(vec2 v) const {
	return getVertex(0).get_position() * v.x + getVertex(1).get_position() * v.y + getVertex(2).get_position() * (1 - v.x - v.y);
}


std::array<vec3, 3> IndexedTriangle::borderTriangle(float width) const {
	vec3 p0 = getVertex(0).get_position();
	vec3 p1 = getVertex(1).get_position();
	vec3 p2 = getVertex(2).get_position();
	vec3 n = normalize(cross(p1 - p0, p2 - p0));
	vec3 v01 = p1 - p0;
	vec3 v12 = p2 - p1;
	vec3 v20 = p0 - p2;
	vec3 w01 = GramSchmidtProcess(mat3(n, v01, v12))[2];
	vec3 w12 = GramSchmidtProcess(mat3(n, v12, v20))[2];
	vec3 w20 = GramSchmidtProcess(mat3(n, v20, v01))[2];
	w01 = w01 * sign(dot(p2 - p0, w01));
	w12 = w12 * sign(dot(p0 - p1, w12));
	w20 = w20 * sign(dot(p1 - p2, w20));

	vec3 p0_01 = p0 + w01 * width;
	vec3 p1_01 = p1 + w01 * width;
	vec3 p1_12 = p1 + w12 * width;
	vec3 p2_12 = p2 + w12 * width;
	vec3 p0_20 = p0 + w20 * width;
	vec3 p2_20 = p2 + w20 * width;

	vec2 b0_01 = toPlanar(p0_01);
	vec2 b1_01 = toPlanar(p1_01);
	vec2 b1_12 = toPlanar(p1_12);
	vec2 b2_12 = toPlanar(p2_12);
	vec2 b0_20 = toPlanar(p0_20);
	vec2 b2_20 = toPlanar(p2_20);

	vec2 b0 = intersectLines(b0_01, b1_01, b0_20, b2_20);
	vec2 b1 = intersectLines(b1_01, b0_01, b1_12, b2_12);
	vec2 b2 = intersectLines(b2_20, b0_20, b2_12, b1_12);

	return {fromPlanar(b0), fromPlanar(b1), fromPlanar(b2)};
}

vec3 IndexedTriangle::faceNormal() const {
	return normalize(cross(getVertex(1).get_position() - getVertex(0).get_position(), getVertex(2).get_position() - getVertex(0).get_position()));
}

vec3 IndexedTriangle::center() const {
	return (getVertex(0).get_position() + getVertex(1).get_position() + getVertex(2).get_position()) / 3.f;
}

float IndexedTriangle::area() const {
	return 0.5f * length(cross(getVertex(1).get_position() - getVertex(0).get_position(), getVertex(2).get_position() - getVertex(0).get_position()));
}

bool IndexedTriangle::containsEdge(int i, int j) const {
	return contains(getVertexIndices(), i) && contains(getVertexIndices(), j);
}

bool IndexedTriangle::containsEdge(ivec2 edge) const {
	return containsEdge(edge.x, edge.y);
}

void IndexedTriangle::setVertexIndices(const ivec3& in) const {
	bufferBoss.setFaceIndices(index, in);
}

void IndexedTriangle::changeOrientation() const {
	bufferBoss.setFaceIndices(index, ivec3(getVertexIndices().z, getVertexIndices().y, getVertexIndices().x));
}


IndexedMesh3D::IndexedMesh3D(IndexedMesh3D&& other) noexcept
: boss(std::move(other.boss)), vertices(std::move(other.vertices)), triangles(std::move(other.triangles)), polygroupIndexOrder(std::move(other.polygroupIndexOrder)) {}

IndexedMesh3D& IndexedMesh3D::operator=(IndexedMesh3D&& other) noexcept {
	if (this == &other)
		return *this;
	boss = std::move(other.boss);
	polygroupIndexOrder = std::move(other.polygroupIndexOrder);
	vertices = std::move(other.vertices);
	triangles = std::move(other.triangles);
	return *this;
}

IndexedMesh3D::IndexedMesh3D() {
	boss = make_unique<BufferManager>();
	vertices = vector<vector<BufferedVertex>>();
	triangles = vector<vector<IndexedTriangle>>();
	polygroupIndexOrder = unordered_map<PolyGroupID, int>();
}

IndexedMesh3D::IndexedMesh3D(const vector<Vertex>& hardVertices, const vector<ivec3>& faceIndices, const PolyGroupID& id)
: IndexedMesh3D() {
	addNewPolygroup(hardVertices, faceIndices, id);
}

IndexedMesh3D::IndexedMesh3D(const char* filename, const PolyGroupID& id)
: IndexedMesh3D() {
	addNewPolygroup(filename, id);
}

void IndexedMesh3D::addNewPolygroup(const char* filename, const PolyGroupID& id) {
	int shift = boss->attributeDataLength();
	int index = polygroupIndexOrder.size();
	if (polygroupIndexOrder.contains(id))
		throw IllegalVariantError("Polygroup ID already exists in mesh. ", __FILE__, __LINE__);
	polygroupIndexOrder[id] = index;
	vertices.emplace_back(vector<BufferedVertex>());
	triangles.emplace_back(vector<IndexedTriangle>());
	vector<vec3> norms = {};
	vector<vec2> uvs = {};
	vector<vec3> pos = {};

	vec2 UV_DEFAULT = vec2(2137, 0);
	vec3 NORM_DEFAULT = vec3(0, 0, 0);

	std::set<int> vertexIndicesFilled = {};

	// auto estimatedSizes = countEstimatedBufferSizesInOBJFile(filename);
	// vertices[id].reserve(estimatedSizes["positions"]);
	// triangles[id].reserve(estimatedSizes["faces"]);
	// norms.reserve(estimatedSizes["normals"]);
	// uvs.reserve(estimatedSizes["uvs"]);


	std::ifstream in(filename, std::ios::in);
	if (!in)
		throw FileNotFoundError(string(filename), __FILE__, __LINE__);

	string line;
	while (getline(in, line)) {
		if (line.substr(0, 2) == "v ") {
			std::istringstream s(line.substr(2));
			vec3 vert;
			s >> vert.x;
			s >> vert.y;
			s >> vert.z;
			pos.emplace_back(vert);
			vertices[index].emplace_back(*boss, Vertex(vert, UV_DEFAULT, NORM_DEFAULT, vec4(1, 1, 1, 1)));
		}
		if (line.substr(0, 3) == "vt ") {
			std::istringstream s(line.substr(3));
			vec2 uv;
			s >> uv.x;
			s >> uv.y;
			uvs.emplace_back(uv);
		}

		if (line.substr(0, 3) == "vn ") {
			std::istringstream s(line.substr(3));
			vec3 norm;
			s >> norm.x;
			s >> norm.y;
			s >> norm.z;
			norms.emplace_back(norm);
		}

		else if (line.substr(0, 2) == "_f ") {
			ivec3 vertexIndex, uvIndex, normalIndex;
			sscanf(line.c_str(), "_f %d/%d/%d %d/%d/%d %d/%d/%d", &vertexIndex[0], &uvIndex[0], &normalIndex[0], &vertexIndex[1], &uvIndex[1], &normalIndex[1], &vertexIndex[2],
				   &uvIndex[2], &normalIndex[2]);

			vertexIndex -= ivec3(1, 1, 1);
			uvIndex -= ivec3(1, 1, 1);
			normalIndex -= ivec3(1, 1, 1);

			vector uvsVec = {uvs.at(uvIndex[0]), uvs.at(uvIndex[1]), uvs.at(uvIndex[2])};
			vector normsVec = {norms.at(normalIndex[0]), norms.at(normalIndex[1]), norms.at(normalIndex[2])};

			ivec3 face = ivec3(0);
			for (int i = 0; i < 3; i++)
				if (!vertexIndicesFilled.contains(vertexIndex[i])) {
					face[i] = vertexIndex[i];
					vertices[index][vertexIndex[i]].setNormal(normsVec[i]);
					vertices[index][vertexIndex[i]].setUV(uvsVec[i]);
					vertexIndicesFilled.insert(vertexIndex[i]);
				}

				else if (vertices[index][vertexIndex[i]].getNormal() == normsVec[i] && vertices[index][vertexIndex[i]].getUV() == uvsVec[i])
					face[i] = vertexIndex[i];

				else {
					vertices[index].emplace_back(*boss, Vertex(pos[vertexIndex[i]], uvsVec[i], normsVec[i], vec4(1, 1, 1, 1)));
					face[i] = vertices[index].back().getIndex();
				}

			triangles[index].emplace_back(*boss, face, shift);
		}
	}
}


IndexedMesh3D& IndexedMesh3D::operator=(const IndexedMesh3D& other) {
	if (this == &other)
		return *this;
	boss = std::make_unique<BufferManager>(*other.boss);
	vertices = other.vertices;
	triangles = other.triangles;
	polygroupIndexOrder = other.polygroupIndexOrder;
	return *this;
}

IndexedMesh3D::IndexedMesh3D(const SmoothParametricSurface& surf, int tRes, int uRes, const PolyGroupID& id)
: IndexedMesh3D() {
	addUniformSurface(surf, tRes, uRes, id);
}


void IndexedMesh3D::addNewPolygroup(const vector<Vertex>& hardVertices, const vector<ivec3>& faceIndices, const PolyGroupID& id) {
	if (polygroupIndexOrder.contains(id))
		throw IllegalVariantError("Polygroup ID already exists in mesh. ", __FILE__, __LINE__);

	int shift = boss->attributeDataLength();
	int index = polygroupIndexOrder.size();
	polygroupIndexOrder[id] = index;
	vertices.emplace_back(vector<BufferedVertex>());
	triangles.emplace_back(vector<IndexedTriangle>());
	vertices[index].reserve(hardVertices.size());
	triangles[index].reserve(faceIndices.size());

	for (const Vertex& v : hardVertices)
		vertices[index].emplace_back(*boss, v);

	for (const ivec3& ind : faceIndices)
		triangles[index].emplace_back(*boss, ind, shift);
}

IndexedMesh3D::IndexedMesh3D(const IndexedMesh3D& other)
: boss(&*other.boss), vertices(other.vertices), polygroupIndexOrder(other.polygroupIndexOrder), triangles(other.triangles) {}


void IndexedMesh3D::addUniformSurface(const SmoothParametricSurface& surf, int tRes, int uRes, const PolyGroupID& id) {
	vector<Vertex> hardVertices = {};
	vector<ivec3> faceIndices = {};
	hardVertices.reserve(tRes * uRes + uRes + tRes + 1);
	faceIndices.reserve(2 * (tRes) * (uRes));


	for (int i = 0; i < tRes; i++)
		for (int j = 0; j < uRes; j++) {
			float t_ = 1.f * i / (tRes - 1);
			if (surf.isPeriodicT())
				t_ = 1.f * i / (1.f * tRes);

			float u_ = 1.f * j / (uRes - 1);
			if (surf.isPeriodicU())
				u_ = 1.f * j / (1.f * uRes);

			float t = lerp(surf.tMin(), surf.tMax(), t_);
			float u = lerp(surf.uMin(), surf.uMax(), u_);
			hardVertices.emplace_back(surf(t, u), vec2(t_, u_), surf.normal(t, u), vec4(t, u, 0, 1));
		}

	ivec2 size = ivec2(tRes, uRes);

	for (int i = 1; i < tRes; i++)
		for (int j = 1; j < uRes; j++) {
			int i0 = flattened2DVectorIndex(i, j, size);
			int i1 = flattened2DVectorIndex(i - 1, j, size);
			int i2 = flattened2DVectorIndex(i - 1, j - 1, size);
			int i3 = flattened2DVectorIndex(i, j - 1, size);

			faceIndices.emplace_back(i0, i1, i2);
			faceIndices.emplace_back(i0, i3, i2);
		}


	if (surf.isPeriodicT())
		for (int j = 1; j < uRes; j++) {
			int i0 = flattened2DVectorIndex(0, j, size);
			int i1 = flattened2DVectorIndex(-1, j, size);
			int i2 = flattened2DVectorIndex(-1, j - 1, size);
			int i3 = flattened2DVectorIndex(0, j - 1, size);
			faceIndices.emplace_back(i0, i1, i2);
			faceIndices.emplace_back(i0, i3, i2);
		}
	//
	if (surf.isPeriodicU())
		for (int i = 1; i < tRes; i++) {
			int i0 = flattened2DVectorIndex(i, 0, size);
			int i1 = flattened2DVectorIndex(i - 1, 0, size);
			int i2 = flattened2DVectorIndex(i - 1, -1, size);
			int i3 = flattened2DVectorIndex(i, -1, size);
			faceIndices.emplace_back(i0, i1, i2);
			faceIndices.emplace_back(i0, i3, i2);
		}

	if (surf.isPeriodicT() && surf.isPeriodicU()) {
		int i0 = flattened2DVectorIndex(0, 0, size);
		int i1 = flattened2DVectorIndex(-1, 0, size);
		int i2 = flattened2DVectorIndex(-1, -1, size);
		int i3 = flattened2DVectorIndex(0, -1, size);
		faceIndices.emplace_back(i0, i1, i2);
		faceIndices.emplace_back(i0, i3, i2);
	}

	addNewPolygroup(hardVertices, faceIndices, id);
}


void IndexedMesh3D::merge(const IndexedMesh3D& other) {
	for (const auto& id : other.getPolyGroupIDs())
		addNewPolygroup(other.getVertices(id), other.getIndices(id), make_unique_id(id));
}

void IndexedMesh3D::mergeAndKeepID(const IndexedMesh3D& other) {
	for (const auto& id : other.getPolyGroupIDs())
		addNewPolygroup(other.getVertices(id), other.getIndices(id), id);
}

void IndexedMesh3D::copyPolygroup(const IndexedMesh3D& other, const PolyGroupID& id, const PolyGroupID& newId) {
	vector<Vertex> vertices = other.getVertices(id);
	vector<ivec3> indices = other.getIndices(id);
	addNewPolygroup(vertices, indices, newId);
}

void IndexedMesh3D::copyPolygroup(const PolyGroupID& id, const PolyGroupID& newId) {
	vector<Vertex> vertices = getVertices(id);
	vector<ivec3> indices = getIndices(id);
	addNewPolygroup(vertices, indices, newId);
}

vector<string> IndexedMesh3D::getActiveExtraBuffers() const {
	return boss->getExtraBufferNames();
}

raw_data_ptr IndexedMesh3D::bufferIndexLocation() const {
	return boss->indexDataPtr();
}

byte_size IndexedMesh3D::faceIndicesDataSize() const {
	return boss->indexDataSize();
}

array_len IndexedMesh3D::bufferIndexLength() const {
	return boss->indexDataLength();
}

raw_data_ptr IndexedMesh3D::getBufferLocation(const string& name) const {
	return boss->attributeDataPtr(name);
}

array_len IndexedMesh3D::getBufferLength() const {
	return boss->attributeDataLength();
}

byte_size IndexedMesh3D::getAttributeDataSize(const string& name) const {
	return boss->attributeDataSize(name);
}

byte_size IndexedMesh3D::totalAttributeDataSize() const {
	return boss->totalAttributeSize();
}

byte_size IndexedMesh3D::totalByteSize() const {
	return totalAttributeDataSize() + faceIndicesDataSize();
}


vector<PolyGroupID> IndexedMesh3D::getPolyGroupIDs() const {
	vector<PolyGroupID> ids;
	ids.reserve(polygroupIndexOrder.size());
	for (const auto& key : polygroupIndexOrder | std::views::keys)
		ids.emplace_back(key);
	return ids;
}

BufferManager& IndexedMesh3D::getBufferBoss() const {
	return *boss;
}

// bool IndexedMesh::hasGlobalTextures() const {
// 	return !isActive(MATERIAL1) && material->textured();
// }

BufferedVertex& IndexedMesh3D::getAnyVertexFromPolyGroup(const PolyGroupID& id) {
	return vertices.at(polygroupIndexOrder[id]).front();
}

void IndexedMesh3D::deformPerVertex(const PolyGroupID& id, const HOM(BufferedVertex&, void)& deformation) {
	for (auto& v : vertices.at(polygroupIndexOrder[id]))
		deformation(v);
}

void IndexedMesh3D::deformPerVertex(const HOM(BufferedVertex&, void)& deformation) {
	for (auto id : getPolyGroupIDs())
		for (auto& v : vertices.at(polygroupIndexOrder[id]))
			deformation(v);
}

void IndexedMesh3D::deformPerVertex(const PolyGroupID& id, const BIHOM(int, BufferedVertex&, void)& deformation) {
	for (int i = 0; i < vertices.at(polygroupIndexOrder[id]).size(); i++)
		deformation(i, vertices.at(polygroupIndexOrder[id])[i]);
}

void IndexedMesh3D::deformPerId(const BIHOM(BufferedVertex&, PolyGroupID, void)& deformation) {
	for (auto id : getPolyGroupIDs())
		for (auto& v : vertices.at(polygroupIndexOrder[id]))
			deformation(v, id);
}

vec2 IndexedMesh3D::getSurfaceParameters(const BufferedVertex& v) {
	return vec2(v.getColor().x, v.getColor().y);
}


void IndexedMesh3D::encodeSurfacePoint(BufferedVertex& v, const SmoothParametricSurface& surf, vec2 tu) {
	v.setPosition(surf(tu));
	v.setNormal(surf.normal(tu));
	v.setUV(vec2((tu.x - surf.tMin()) / (surf.tMax() - surf.tMin()), (tu.y - surf.uMin()) / (surf.uMax() - surf.uMin())));
	v.setColor(tu.x, 0);
	v.setColor(tu.y, 1);
}

void IndexedMesh3D::adjustToNewSurface(const SmoothParametricSurface& surf, const PolyGroupID& id) {
	deformPerVertex(id, [&](BufferedVertex& v) {
		encodeSurfacePoint(v, surf, getSurfaceParameters(v));
	});
}

void IndexedMesh3D::adjustToNewSurface(const SmoothParametricSurface& surf) {
	for (auto id : getPolyGroupIDs())
		adjustToNewSurface(surf, id);
}

void IndexedMesh3D::moveAlongVectorField(const PolyGroupID& id, const VectorField& X, float delta) {
	for (BufferedVertex& v : vertices.at(polygroupIndexOrder[id]))
		v.setPosition(X.moveAlong(v.getPosition(), delta));
}

void IndexedMesh3D::deformWithAmbientMap(const PolyGroupID& id, const SpaceEndomorphism& f) {
	for (BufferedVertex& v : vertices.at(polygroupIndexOrder[id]))
		v.applyFunction(f);
}

void IndexedMesh3D::deformWithAmbientMap(const SpaceEndomorphism& f) {
	for (auto id : getPolyGroupIDs())
		deformWithAmbientMap(id, f);
}


void IndexedMesh3D::affineTransform(const mat3& M, vec3 v, const PolyGroupID& id) {
	deformWithAmbientMap(id, SpaceEndomorphism::affine(M, v));
}

void IndexedMesh3D::affineTransform(const mat3& M, vec3 v) {
	for (auto& name : getPolyGroupIDs())
		affineTransform(M, v, name);
}

void IndexedMesh3D::shift(vec3 v, const PolyGroupID& id) {
	affineTransform(mat3(1), v, id);
}

void IndexedMesh3D::shift(vec3 v) {
	affineTransform(mat3(1), v);
}

void IndexedMesh3D::scale(float s, const PolyGroupID& id) {
	affineTransform(mat3(s), vec3(0), id);
}

void IndexedMesh3D::scale(float s) {
	affineTransform(mat3(s), vec3(0));
}

vector<Vertex> IndexedMesh3D::getVertices(const PolyGroupID& id) const {
	vector<Vertex> verts = {};
	verts.reserve(vertices.at(polygroupIndexOrder.at(id)).size());
	for (const BufferedVertex& v : vertices.at(polygroupIndexOrder.at(id)))
		verts.push_back(v.getVertex());
	return verts;
}

vector<BufferedVertex> IndexedMesh3D::getBufferedVertices(const PolyGroupID& id) const {
	return vertices.at(polygroupIndexOrder.at(id));
}


vector<ivec3> IndexedMesh3D::getIndices(const PolyGroupID& id) const {
	vector<ivec3> inds = {};
	inds.reserve(triangles.at(polygroupIndexOrder.at(id)).size());
	for (const IndexedTriangle& t : triangles.at(polygroupIndexOrder.at(id)))
		inds.push_back(t.getVertexIndices());
	return inds;
}

vector<IndexedTriangle> IndexedMesh3D::getTriangles(const PolyGroupID& id) const {
	return triangles.at(polygroupIndexOrder.at(id));
}

// vec4 IndexedMesh::getIntencities() const {
// 	return material->compressIntencities();
// }
//
// MaterialPhong IndexedMesh::getMaterial() const {
// 	return *material;
// }
//
//
// void IndexedMesh::addGlobalMaterial(const MaterialPhong &mat) {
// 	material = std::make_shared<MaterialPhong>(mat);
// }

void IndexedMesh3D::flipNormals(const PolyGroupID& id) {
	deformPerVertex(id, [](BufferedVertex& v) {
		v.setNormal(-v.getNormal());
	});
}

void IndexedMesh3D::flipNormals() {
	for (auto& name : getPolyGroupIDs())
		flipNormals(name);
}

void IndexedMesh3D::pointNormalsInDirection(vec3 dir, const PolyGroupID& id) {
	deformPerVertex(id, [dir](BufferedVertex& v) {
		vec3 n = v.getNormal();
		v.setNormal(n * 1.f * sgn(dot(n, dir)));
	});
}

void IndexedMesh3D::pointNormalsInDirection(vec3 dir) {
	for (auto& name : getPolyGroupIDs())
		pointNormalsInDirection(dir, name);
}

IndexedMesh3D IndexedMesh3D::subdivideBarycentric(const PolyGroupID& id) const {
	vector<ivec3> trInds = getIndices(id);
	vector<Vertex> verts = getVertices(id);
	verts.reserve(verts.size() + trInds.size());
	vector<ivec3> newInds = {};
	newInds.reserve(3 * trInds.size());
	for (ivec3 tr : trInds) {
		Vertex b = barycenter(verts[tr.x], verts[tr.y], verts[tr.z]);
		verts.push_back(b);
		int centerInd = verts.size() - 1;
		newInds.push_back(ivec3(tr.x, tr.y, centerInd));
		newInds.push_back(ivec3(tr.y, tr.z, centerInd));
		newInds.emplace_back(tr.z, tr.x, centerInd);
	}
	return IndexedMesh3D(verts, newInds, id);
}


IndexedMesh3D IndexedMesh3D::subdivideEdgecentric(const PolyGroupID& id) const {
	vector<ivec3> trInds = getIndices(id);
	vector<Vertex> verts = getVertices(id);
	verts.reserve(verts.size() + 3 * trInds.size());
	vector<ivec3> newInds = {};
	newInds.reserve(4 * trInds.size());
	for (ivec3 tr : trInds) {
		verts.push_back(center(verts[tr.x], verts[tr.y]));
		int center1 = verts.size() - 1;
		verts.push_back(center(verts[tr.x], verts[tr.z]));
		int center2 = verts.size() - 1;
		verts.push_back(center(verts[tr.y], verts[tr.z]));
		int center3 = verts.size() - 1;
		newInds.emplace_back(tr.x, center1, center2);
		newInds.emplace_back(tr.y, center1, center3);
		newInds.emplace_back(tr.z, center2, center3);
		newInds.emplace_back(center1, center2, center3);
	}
	return IndexedMesh3D(verts, newInds, id);
}

IndexedMesh3D IndexedMesh3D::wireframe(PolyGroupID id, PolyGroupID targetId, float width, float heightCenter, float heightSide) const {
	vector<IndexedTriangle> trs = triangles.at(polygroupIndexOrder.at(id));
	vector<BufferedVertex> verts = vertices.at(polygroupIndexOrder.at(id));
	vector<Vertex> new_verts = getVertices(id);
	vector<ivec3> new_inds = {};
	for (const auto& v : verts) {
		Vertex extr = v.getVertex();
		extr.set_position(extr.get_position() + extr.get_normal() * heightCenter);
		new_verts.push_back(extr);
	}
	for (const auto& tr : trs) {
		array<vec3, 3> border = tr.borderTriangle(width);
		vec3 n = tr.faceNormal() * sign(dot(tr.getVertex(0).get_normal(), tr.faceNormal()));


		for (vec3 p : border)
			new_verts.emplace_back(p + n * heightSide, vec2(0, 0), n, BLACK);
		for (vec3 p : border)
			new_verts.emplace_back(p, vec2(0, 0), n, BLACK);
		new_inds.push_back(ivec3(tr.getVertexIndices().x, tr.getVertexIndices().y, tr.getVertexIndices().x + verts.size()));
		new_inds.push_back(ivec3(tr.getVertexIndices().y, tr.getVertexIndices().y + verts.size(), tr.getVertexIndices().x + verts.size()));
		new_inds.push_back(ivec3(tr.getVertexIndices().x + verts.size(), tr.getVertexIndices().y + verts.size(), new_verts.size() - 6));
		new_inds.push_back(ivec3(tr.getVertexIndices().y + verts.size(), new_verts.size() - 5, new_verts.size() - 6));
		// new_inds.push_back(ivec3(new_verts.size()-6,new_verts.size()-5, new_verts.size()-3));
		// new_inds.push_back(ivec3(new_verts.size()-5,new_verts.size()-2, new_verts.size()-3));

		new_inds.push_back(ivec3(tr.getVertexIndices().x, tr.getVertexIndices().z, tr.getVertexIndices().x + verts.size()));
		new_inds.push_back(ivec3(tr.getVertexIndices().z, tr.getVertexIndices().z + verts.size(), tr.getVertexIndices().x + verts.size()));
		// new_inds.push_back(ivec3(tr.getLocalIndices().x+verts.size(), tr.getLocalIndices().z+verts.size(), new_verts.size()-6));
		// new_inds.push_back(ivec3(tr.getLocalIndices().z+verts.size(), new_verts.size()-4, new_verts.size()-6));
		// new_inds.push_back(ivec3(new_verts.size()-6,new_verts.size()-4, new_verts.size()-3));
		// new_inds.push_back(ivec3(new_verts.size()-4,new_verts.size()-1, new_verts.size()-3));

		new_inds.push_back(ivec3(tr.getVertexIndices().y, tr.getVertexIndices().z, tr.getVertexIndices().y + verts.size()));
		new_inds.push_back(ivec3(tr.getVertexIndices().z, tr.getVertexIndices().z + verts.size(), tr.getVertexIndices().y + verts.size()));
		// new_inds.push_back(ivec3(tr.getLocalIndices().y+verts.size(), tr.getLocalIndices().z+verts.size(), new_verts.size()-5));
		// new_inds.push_back(ivec3(tr.getLocalIndices().z+verts.size(), new_verts.size()-4, new_verts.size()-5));
		// new_inds.push_back(ivec3(new_verts.size()-5,new_verts.size()-4, new_verts.size()-2));
		// new_inds.push_back(ivec3(new_verts.size()-4,new_verts.size()-1, new_verts.size()-2));
	}
	return IndexedMesh3D(new_verts, new_inds, targetId);
}


BufferManager::BufferManager(const BufferManager& other)
: stds(std::make_unique<Stds>(*other.stds)), extra0(std::make_unique<vector<vec4>>(*other.extra0)), extra1(std::make_unique<BUFF4>(*other.extra1)),
  extra2(std::make_unique<BUFF4>(*other.extra2)), extra3(std::make_unique<BUFF4>(*other.extra3)), extra4(std::make_unique<BUFF4>(*other.extra4)),
  indices(std::make_unique<IBUFF3>(*other.indices)), extraBufferNames(other.extraBufferNames) {}

BufferManager& BufferManager::operator=(BufferManager&& other) noexcept {
	if (this == &other)
		return *this;
	stds = std::move(other.stds);
	extra0 = std::move(other.extra0);
	extra1 = std::move(other.extra1);
	extra2 = std::move(other.extra2);
	extra3 = std::move(other.extra3);
	extra4 = std::move(other.extra4);
	indices = std::move(other.indices);
	extraBufferNames = std::move(other.extraBufferNames);
	return *this;
}


BufferManager::BufferManager(BufferManager&& other) noexcept
: stds(std::move(other.stds)), extra0(std::move(other.extra0)), extra1(std::move(other.extra1)), extra2(std::move(other.extra2)), extra3(std::move(other.extra3)),
  extra4(std::move(other.extra4)), indices(std::move(other.indices)),extraBufferNames(std::move(other.extraBufferNames)) {}

BufferManager& BufferManager::operator=(const BufferManager& other) {
	if (this == &other)
		return *this;
	stds = std::make_unique<Stds>(*other.stds);
	extra0 = std::make_unique<BUFF4>(*other.extra0);
	extra1 = std::make_unique<BUFF4>(*other.extra1);
	extra2 = std::make_unique<BUFF4>(*other.extra2);
	extra3 = std::make_unique<BUFF4>(*other.extra3);
	extra4 = std::make_unique<BUFF4>(*other.extra4);
	indices = std::make_unique<IBUFF3>(*other.indices);
	extraBufferNames = other.extraBufferNames;
	return *this;
}

void BufferManager::reserveAdditionalSpaceForIndex(int extraStorage) {
	reserveSpaceForIndex(indexDataLength() + extraStorage);
}

vec3 BufferManager::getPosition(int index) const {
	return stds->positions[index];
}

vec3 BufferManager::getNormal(int index) const {
	return stds->normals[index];
}


vec2 BufferManager::getUV(int index) const {
	return stds->uvs[index];
}

vec4 BufferManager::getColor(int index) const {
	return stds->colors[index];
}

vec4 BufferManager::getExtra(int index, const string& name) const {
	return extraElementAt(index, indexOfExtraBuffer(name));

}

ivec3 BufferManager::getFaceIndices(int index) const {
	return (*indices)[index];
}

Vertex BufferManager::getVertex(int index) const {
	dict(string, vec4) extraData = {};
	for (const string& name : extraBufferNames)
		extraData[name] = getExtra(index, name);
	return Vertex(getPosition(index), getUV(index), getNormal(index), getColor(index), extraData);
}

void BufferManager::setPosition(int index, vec3 value) {
	stds->positions[index] = value;
}

void BufferManager::setNormal(int index, vec3 value) {
	stds->normals[index] = value;
}

void BufferManager::setUV(int index, vec2 value) {
	stds->uvs[index] = value;
}

void BufferManager::setUV(int index, float value, int coord) {
	stds->uvs[index][coord] = value;
}

void BufferManager::setColor(int index, vec4 value) {
	stds->colors[index] = value;
}

void BufferManager::setColor(int index, float value, int component) {
	stds->colors[index][component] = value;
}

void BufferManager::setFaceIndices(int index, const ivec3& in) {
	indices->at(index) = in;
}

void BufferManager::setExtra(int index, vec4 value, const string& name) {
	int slot = indexOfExtraBuffer(name);
	setElementAt(index, value, slot);
}

void BufferManager::setExtra(int index, vec3 value, const string& name) {
	int slot = indexOfExtraBuffer(name);
	setElementAt(index, value[0], slot, 0);
	setElementAt(index, value[1], slot, 1);
	setElementAt(index, value[2], slot, 2);
}

void BufferManager::setExtra(int index, float value, const string& name, int component) {
	int slot = indexOfExtraBuffer(name);
	setElementAt(index, value, slot, component);
}

void BufferManager::reserveSpace(int targetSize) const {
	stds->positions.reserve(targetSize);
	stds->normals.reserve(targetSize);
	stds->uvs.reserve(targetSize);
	stds->colors.reserve(targetSize);

	if (extra0)
		extra0->reserve(targetSize);
	if (extra1)
		extra1->reserve(targetSize);
	if (extra2)
		extra2->reserve(targetSize);
	if (extra3)
		extra3->reserve(targetSize);
	if (extra4)
		extra4->reserve(targetSize);
}

void BufferManager::reserveSpaceForIndex(int targetSize) const {
	indices->reserve(targetSize);
}

void BufferManager::reserveAdditionalSpace(int extraStorage) const {
	reserveSpace(attributeDataLength() + extraStorage);
}

vec3 IndexedMesh3D::centerOfMass(PolyGroupID id) const {
	vec3 sum = vec3(0);
	float totalArea = 0;
	for (const IndexedTriangle& t : triangles.at(polygroupIndexOrder.at(id))) {
		sum += t.center() * t.area();
		totalArea += t.area();
	}
	return sum / totalArea;
}

vec3 IndexedMesh3D::centerOfMass() const {
	vec3 sum = vec3(0);
	float totalArea = 0;
	for (const auto& id : getPolyGroupIDs())
		for (const IndexedTriangle& t : triangles.at(polygroupIndexOrder.at(id))) {
			sum += t.center() * t.area();
			totalArea += t.area();
		}
	return sum / totalArea;
}


Wireframe::Wireframe(const SmoothParametricSurface& surf, float width, int n, int m, int curve_res_rad, int curve_res_hor)
: IndexedMesh3D(), width(width), surf(surf), n(n), m(m), curve_res_rad(curve_res_rad), curve_res_hor(curve_res_hor) {
	for (float t_i : linspace(surf.tMin(), surf.tMax(), n)) {
		SmoothParametricCurve curve = surf.constT(t_i);
		auto id = randomID();
		addUniformSurface(curve.pipe(width), curve_res_rad, curve_res_hor, id);
		deformPerVertex(id, [t_i](BufferedVertex& v) {
			v.setColor(t_i, 2);
			v.setColor(v.getColor().x, 3);
		});
	}

	for (float u_i : linspace(surf.uMin(), surf.uMax(), m)) {
		SmoothParametricCurve curve = surf.constU(u_i);
		auto id = randomID();
		addUniformSurface(curve.pipe(width), curve_res_rad, curve_res_hor, id);
		deformPerVertex(id, [u_i](BufferedVertex& v) {
			v.setColor(u_i, 3);
			v.setColor(v.getColor().x, 2);
		});
	}
}


void Wireframe::changeBaseSurface(const SmoothParametricSurface& newsurf) {
	for (auto id : getPolyGroupIDs())
		for (BufferedVertex& v : vertices.at(polygroupIndexOrder.at(id))) {
			vec2 tu = getSurfaceParameters(v);
			vec3 old_p0 = surf(tu);
			vec3 new_p0 = newsurf(tu);

			vec3 old_n0 = surf.normal(tu);
			vec3 new_n0 = newsurf.normal(tu);
			v.applyFunction(SpaceEndomorphism::affine(mat3(1), new_p0 - old_p0));
			v.applyFunction(SpaceAutomorphism::deltaRotation(old_n0, new_n0, v.getPosition()));
		}
	surf = newsurf;
}

vec2 Wireframe::getSurfaceParameters(const BufferedVertex& v) const {
	return vec2(v.getColor().z, v.getColor().w);
}

PlanarFlowLines::PlanarFlowLines(const VectorFieldR2& X, float dt, int steps, const std::function<float(float, float, float, vec2, vec2)>& width,
								 const std::function<vec4(float, float, float, vec2, vec2)>& color)
: IndexedMesh3D(), X(X), dt(dt), steps(steps), width(width), color(color) {
	boss = make_unique<BufferManager>("extra0");
}


void PlanarFlowLines::generateGrid(vec2 v_min, vec2 v_max, ivec2 res) {
	startPoints = flatten(linspace2D(v_min, vec2((v_max - v_min).x / res.x, 0), vec2(0, (v_max - v_min).y / res.y), res.x + 1, res.y + 1));
	ids = {};
	ids.reserve(startPoints.size());
	for (int i = 0; i < startPoints.size(); i++)
		ids.push_back(randomID());
}


void PlanarFlowLines::generateRandomUniform(vec2 v_min, vec2 v_max, int n) {
	startPoints = {};
	startPoints.reserve(n);
	ids = {};
	ids.reserve(n);
	for (int i = 0; i < n; i++) {
		startPoints.emplace_back(randomFloat(v_min.x, v_max.x), randomFloat(v_min.y, v_max.y));
		ids.push_back(randomID());
	}
}

void PlanarFlowLines::generateStartTimesAll0() {
	startTimes = vector<float>(startPoints.size(), 0);
}

void PlanarFlowLines::generateStartTimesUniform(float t_max) {
	for (int i = 0; i < startPoints.size(); i++)
		startTimes.push_back(randomFloat(0, t_max));
}

void PlanarFlowLines::generateLine(int i) {
	float t0 = startTimes[i];
	vec2 p0 = startPoints[i];
	vec2 p = startPoints[i];
	vector<Vertex> verts = {};
	vector<ivec3> trs = {};
	float len = 0;

	for (int j = 1; j <= steps; j++) {
		len += norm(X(p) * dt);
		p += X(p) * dt;
		vec2 n = X.normal(p);
		float t = t0 + j * dt;
		float speed = X.speed(p);
		float w = width(t, t0, speed, p, p0);
		vec2 q0 = p - n * w;
		vec2 q1 = p + n * w;
		vec4 c = color(t, t0, speed, p, p0);
		c.w = speed;
		verts.emplace_back(vec3(q0, t / (steps * dt * 20)), vec2(t, w), e3, c, dict(string, vec4)({{"extra0", vec4(t0, p0.x, p0.y, len)}}));
		verts.emplace_back(vec3(p, t / (steps * dt * 20)), vec2(t, 0), e3, c, dict(string, vec4)({{"extra0", vec4(t0, p0.x, p0.y, len)}}));
		verts.emplace_back(vec3(q1, t / (steps * dt * 20)), vec2(t, -w), e3, c, dict(string, vec4)({{"extra0", vec4(t0, p0.x, p0.y, len)}}));


		if (j > 1) {
			int last_ind = verts.size() - 1;
			trs.emplace_back(last_ind, last_ind - 3, last_ind - 4);
			trs.emplace_back(last_ind, last_ind - 1, last_ind - 4);
			trs.emplace_back(last_ind - 1, last_ind - 5, last_ind - 4);
			trs.emplace_back(last_ind - 1, last_ind - 2, last_ind - 5);
		}
	}

	addNewPolygroup(verts, trs, ids[i]);
}

void PlanarFlowLines::generateLines() {
	for (int i = 0; i < startPoints.size(); i++)
		generateLine(i);
}

float PlanarFlowLines::getTimeRelative(const BufferedVertex& v) {
	return v.getUV().x;
}

float PlanarFlowLines::getT0(const BufferedVertex& v) {
	return v.getExtra("extra0").x;
}

float PlanarFlowLines::getTimeAbsolute(const BufferedVertex& v) {
	return getT0(v) + getTimeRelative(v);
}

vec2 PlanarFlowLines::getPos(const BufferedVertex& v) {
	return vec2(v.getPosition());
}

vec2 PlanarFlowLines::getStartPoint(const BufferedVertex& v) {
	return vec2(v.getExtra("extra0").y, v.getExtra("extra0").z);
}

float PlanarFlowLines::getSpeed(const BufferedVertex& v) {
	return v.getColor().w;
}

float PlanarFlowLines::getLength(const BufferedVertex& v) {
	return v.getExtra("extra0").w;
}

vec4 PlanarFlowLines::getColor(const BufferedVertex& v) {
	return vec4(vec3(v.getColor()), 1);
}

float PlanarFlowLines::getWidth(const BufferedVertex& v) {
	return v.getUV().y;
}

PlanarDiffusedInterval::PlanarDiffusedInterval(const VectorFieldR2& X, float dt, int steps, vec4 color, const HOM(float, float)& time_dump,
											   const HOM(float, float)& width_ratio_dump, vec2 a, vec2 b, float t0)
: X(X), dt(dt), steps(steps), color(color), width_ratio_dump(width_ratio_dump), time_dump(time_dump), a(a), b(b), t0(t0) {
	vec2 p_a = a;
	vec2 p_b = b;
	vector<Vertex> verts = {};
	vector<ivec3> trs = {};
	float t = t0;
	float eps = 0.1;
	float width = length(b - a);
	for (int j = 0; j <= steps; j++) {
		float width_ratio = width / length(p_b - p_a);
		float t_dump = time_dump(t - t0);
		float w_dump = width_ratio_dump(width_ratio);
		vec4 col = color;
		col.w = t_dump * w_dump;
		if (t_dump * w_dump < .05)
			col.w = 0.f;
		verts.emplace_back(vec3(p_a, (t - t0) * eps), vec2(t, t0), e3, col);
		verts.emplace_back(vec3(p_b, (t - t0) * eps), vec2(t, t0), e3, col);
		p_a += X(p_a) * dt;
		p_b += X(p_b) * dt;
		t += dt;
		if (j > 0 && width_ratio > .15) {
			int last_ind = verts.size() - 1;
			trs.emplace_back(last_ind, last_ind - 3, last_ind - 2);
			trs.emplace_back(last_ind, last_ind - 1, last_ind - 3);
		}
	}
	addNewPolygroup(verts, trs, id);
}

PlanarDiffusedCurve::PlanarDiffusedCurve(const VectorFieldR2& X, float dt, int steps, vec4 color, const HOM(float, float)& time_dump, const HOM(float, float)& width_ratio_dump,
										 const SmoothParametricPlaneCurve& curve, vec2 domain, int resolution, float t0)
: X(X), dt(dt), steps(steps), color(color), width_ratio_dump(width_ratio_dump), time_dump(time_dump), domain(domain), resolution(resolution), t0(t0), curve(curve) {
	for (int i = 0; i < resolution; i++) {
		float a = lerp(domain.x, domain.y, 1.f * i / resolution);
		float b = lerp(domain.x, domain.y, 1.f * (i + 1) / resolution);
		auto line = PlanarDiffusedInterval(X, dt, steps, color, time_dump, width_ratio_dump, curve(a), curve(b), t0);
		merge(line);
	}
}

PlanarDiffusedPatterns::PlanarDiffusedPatterns(const VectorFieldR2& X, float dt, int steps, const vector<vec4>& colors, const vector<vec2>& shifts,
											   const HOM(float, float)& time_dump, const HOM(float, float)& width_ratio_dump, const SmoothParametricPlaneCurve& curve_pattern,
											   vec2 domain, int resolution, float t0) {
	for (int i = 0; i < colors.size(); i++) {
		auto c = curve_pattern + shifts[i];
		auto line = PlanarDiffusedCurve(X, dt, steps, colors[i], time_dump, width_ratio_dump, c, domain, resolution, t0);
		merge(line);
	}
}


PipeCurveVertexShader::PipeCurveVertexShader(const SmoothParametricCurve& curve, float r, int horRes, int radialRes, const PolyGroupID& id, dict(string, vec4) extra_defaults)
: PipeCurveVertexShader(curve, PIPE_SETTINGS(r, horRes, radialRes, extra_defaults), id) {}

PipeCurveVertexShader::PipeCurveVertexShader(const RealFunction& plot, vec2 dom, float r, int horRes, int radialRes, const PolyGroupID& id)
: id(id), settings(r, horRes, radialRes) {
	boss = make_unique<BufferManager>("extra0");
	auto params = linspace(dom.x, dom.y, horRes);

	vector<Vertex> verts = vector<Vertex>();
	vector<ivec3> inds = vector<ivec3>();


	for (int i = 0; i < horRes; i++) {
		float t = params[i];
		vec3 p = vec3(t, 0, plot(t));
		vec3 b = e2;
		vec3 n = normalise(vec3(plot.df(t), 0, -1));
		for (float theta : linspace(0.f, TAU, radialRes))
			verts.emplace_back(p, vec2(t, theta), n, vec4(b.x, b.y, b.z, r), dict(string, vec4)({{"extra0", vec4(0)}}));

		if (i < horRes - 1) {
			for (int j = 0; j < radialRes; j++) {
				inds.emplace_back(i * radialRes + j, (i + 1) * radialRes + j, i * radialRes + (j + 1) % radialRes);
				inds.emplace_back((i + 1) * radialRes + j, (i + 1) * radialRes + (j + 1) % radialRes, i * radialRes + (j + 1) % radialRes);
			}
		}
	}
	addNewPolygroup(verts, inds, this->id);
}

PipeCurveVertexShader::PipeCurveVertexShader(const DiscreteRealFunction& plot, float r, int radialRes, const PolyGroupID& id)
: id(id), settings(r, plot.samples(), radialRes) {
	boss = make_unique<BufferManager>("extra0");


	vec2 dom = plot.getDomain();
	auto horRes = plot.samples();
	auto params = linspace(dom.x, dom.y, horRes);
	auto df = plot.derivative();

	vector<Vertex> verts = vector<Vertex>();
	vector<ivec3> inds = vector<ivec3>();


	for (int i = 0; i < horRes; i++) {
		float t = params[i];
		vec3 p = vec3(t, 0, plot(t));
		vec3 b = e2;
		vec3 n = normalise(vec3(df(t), 0, -1));
		for (float theta : linspace(0.f, TAU, radialRes))
			verts.emplace_back(p, vec2(t, theta), n, vec4(b.x, b.y, b.z, r), dict(string, vec4)({{"extra0", vec4(0)}}));

		if (i < horRes - 1) {
			for (int j = 0; j < radialRes; j++) {
				inds.emplace_back(i * radialRes + j, (i + 1) * radialRes + j, i * radialRes + (j + 1) % radialRes);
				inds.emplace_back((i + 1) * radialRes + j, (i + 1) * radialRes + (j + 1) % radialRes, i * radialRes + (j + 1) % radialRes);
			}
		}
	}
	addNewPolygroup(verts, inds, id);
}

PipeCurveVertexShader::PipeCurveVertexShader(const DiscreteRealFunctionNonUniform& plot, float r, int radialRes, const PolyGroupID& id)
: id(id), settings(r, plot.samples(), radialRes) {
	boss = make_unique<BufferManager>("extra0");

	auto params = plot.args_vector();
	int horRes = params.size();
	auto df = plot.derivative();

	vector<Vertex> verts = vector<Vertex>();
	vector<ivec3> inds = vector<ivec3>();


	for (int i = 0; i < horRes; i++) {
		float t = params[i];
		vec3 p = vec3(t, plot(t), 0);
		vec3 b = e3;
		vec3 n = normalise(vec3(df(t), -1, 0));
		for (float theta : linspace(0.f, TAU, radialRes)) {
			verts.emplace_back(p, vec2(t, theta), n, vec4(b.x, b.y, b.z, r), dict(string, vec4)({{"extra0", vec4(0)}}));
		}
		if (i < horRes - 1) {
			for (int j = 0; j < radialRes; j++) {
				inds.emplace_back(i * radialRes + j, (i + 1) * radialRes + j, i * radialRes + (j + 1) % radialRes);
				inds.emplace_back((i + 1) * radialRes + j, (i + 1) * radialRes + (j + 1) % radialRes, i * radialRes + (j + 1) % radialRes);
			}
		}
	}
	addNewPolygroup(verts, inds, id);
}

PipeCurveVertexShader::PipeCurveVertexShader(const SmoothParametricCurve& curve, const PIPE_SETTINGS& s, const std::variant<int, std::string>& id)
: id(id), settings(s) {

	boss = make_unique<BufferManager>("extra0");

	auto params = linspace(curve.getT0(), curve.getT1(), s.horRes);
	vector<vec3> normals = vector<vec3>();
	vector<vec3> binormals = vector<vec3>();
	vector<vec3> pts = vector<vec3>();
	vector<float> realised_params = vector<float>();

	bool out_of_bounds = false;

	for (float t : params) {
		vec3 p = curve(t);
		if (s.bounding) {
			if (p.x <= s.bound_min.x || p.x >= s.bound_max.x || p.y <= s.bound_min.y || p.y >= s.bound_max.y || p.z <= s.bound_min.z || p.z >= s.bound_max.z) {
				if (out_of_bounds)
					continue;
				p = vec3(min(max(p.x, s.bound_min.x), s.bound_max.x), min(max(p.y, s.bound_min.y), s.bound_max.y), min(max(p.z, s.bound_min.z), s.bound_max.z));
				out_of_bounds = true;
			}
			else
				out_of_bounds = false;
		}
		vec3 b = curve.binormal(t);
		vec3 n = curve.normal(t);

		if (normals.size() > 0) {
			if (dot(n, normals.back()) < 0)
				n = -n;
			if (dot(b, binormals.back()) < 0)
				b = -b;
		}
		normals.push_back(n);
		binormals.push_back(b);
		pts.push_back(p);
		realised_params.push_back(t);
	}

	vector<Vertex> verts = vector<Vertex>();
	vector<ivec3> inds = vector<ivec3>();


	for (int i = 0; i < realised_params.size(); i++) {
		float t = realised_params[i];
		vec3 p = pts[i];
		vec3 b = binormals[i];
		vec3 n = normals[i];
		for (float theta : linspace(0.f, TAU, s.radialRes))
			verts.emplace_back(p, vec2(t, theta), n, vec4(b.x, b.y, b.z, s.radius), s.extra_defaults);


		if (i < realised_params.size() - 1) {
			float t_next = realised_params[i + 1];
			bool skip = false;
			for (float d : s.discontinuities)
				if (t_next >= d && t <= d)
					skip = true;

			if (skip)
				continue;

			for (int j = 0; j < s.radialRes; j++) {
				inds.emplace_back(i * s.radialRes + j, (i + 1) * s.radialRes + j, i * s.radialRes + (j + 1) % s.radialRes);
				inds.emplace_back((i + 1) * s.radialRes + j, (i + 1) * s.radialRes + (j + 1) % s.radialRes, i * s.radialRes + (j + 1) % s.radialRes);
			}
		}
	}
	addNewPolygroup(verts, inds, this->id);
}

PipeCurveVertexShader::PipeCurveVertexShader(const RealFunction& plot, vec2 dom, const PIPE_SETTINGS& s, const std::variant<int, std::string>& id)
: PipeCurveVertexShader(plot, dom, s.radius, s.horRes, s.radialRes, id) {}

PipeCurveVertexShader::PipeCurveVertexShader(const DiscreteRealFunction& plot, const PIPE_SETTINGS& s, const std::variant<int, std::string>& id)
: PipeCurveVertexShader(plot, s.radius, s.radialRes, id) {}

PipeCurveVertexShader::PipeCurveVertexShader(const DiscreteRealFunctionNonUniform& plot, const PIPE_SETTINGS& s, const std::variant<int, std::string>& id)
: PipeCurveVertexShader(plot, s.radius, s.radialRes, id) {}


void PipeCurveVertexShader::updateCurve(const SmoothParametricCurve& curve) {
	deformPerVertex(id, [this, &curve](BufferedVertex& v) {
		float t = v.getUV().x;
		vec3 p = curve(t);
		p = vec3(min(max(p.x, settings.bound_min.x), settings.bound_max.x), min(max(p.y, settings.bound_min.y), settings.bound_max.y),
				 min(max(p.z, settings.bound_min.z), settings.bound_max.z));
		vec3 b = curve.binormal(t);
		vec3 n = curve.normal(t);
		// if (dot(n, v.get_normal()) < 0)
		// 	n = -n;
		// if (dot(b, vec3(v.get_color())) < 0)
		// 	b = -b;
		v.setPosition(p);
		v.setNormal(n);
		v.setColor(b.x, 0);
		v.setColor(b.y, 1);
		v.setColor(b.z, 2);
	});
}

void PipeCurveVertexShader::duplicateCurve(const PolyGroupID& copy_id) {
	copyPolygroup(id, copy_id);
}

void PipeCurveVertexShader::updateCurve(const DiscreteRealFunction& plot) {
	auto df = plot.derivative();
	deformPerVertex(id, [this, &plot, &df](BufferedVertex& v) {
		float t = v.getUV().x;
		vec3 p = vec3(t, 0, plot(t));
		vec3 n = normalize(vec3(df(t), 0, -1));
		v.setPosition(p);
		v.setNormal(n);
	});
}


void PipeCurveVertexShader::updateCurve(const DiscreteRealFunction& plot, const DiscreteRealFunction& df_precomputed) {
	deformPerVertex(id, [this, &plot, &df_precomputed](BufferedVertex& v) {
		float t = v.getUV().x;
		vec3 p = vec3(t, 0, plot(t));
		vec3 n = normalize(vec3(df_precomputed(t), 0, -1));
		v.setPosition(p);
		v.setNormal(n);
	});
}

void PipeCurveVertexShader::updateCurve(const DiscreteRealFunctionNonUniform& plot) {
	auto df = plot.derivative();
	deformPerVertex(id, [this, &plot, &df](BufferedVertex& v) {
		float t = v.getUV().x;
		vec3 p = vec3(t, 0, plot(t));
		vec3 n = normalize(vec3(df(t), 0, -1));
		v.setPosition(p);
		v.setNormal(n);
	});
}

void PipeCurveVertexShader::updateRadius(const HOM(float, float)& r) {
	deformPerVertex(id, [this, &r](BufferedVertex& v) {
		float t = v.getUV().x;
		v.setColor(r(t), 3);
	});
}

void PipeCurveVertexShader::updateRadius(float r) {
	deformPerVertex(id, [this, &r](BufferedVertex& v) {
		v.setColor(r, 3); // set radius in color
	});
	settings.radius = r;
}

vec3 PipeCurveVertexShader::getBinormal(const BufferedVertex& v) {
	return vec3(v.getColor());
}

float PipeCurveVertexShader::getAngle(const BufferedVertex& v) {
	return v.getUV().y;
}

float PipeCurveVertexShader::getRadius(const BufferedVertex& v) {
	return v.getColor().w;
}

float PipeCurveVertexShader::getParameter(const BufferedVertex& v) {
	return v.getUV().x;
}

vec4 PipeCurveVertexShader::getExtra(const BufferedVertex& v) {
	return v.getExtra("extra0");
}

void PipeCurveVertexShader::setBinormal(BufferedVertex& v, vec3 b) {
	v.setColor(b.x, 0);
	v.setColor(b.y, 1);
	v.setColor(b.z, 2);
}

void PipeCurveVertexShader::setAngle(BufferedVertex& v, float theta) {
	v.setUV(theta, 1);
}

void PipeCurveVertexShader::setRadius(BufferedVertex& v, float r) {
	v.setColor(r, 3);
}

void PipeCurveVertexShader::setParameter(BufferedVertex& v, float t) {
	v.setUV(t, 0);
}

void PipeCurveVertexShader::setExtra(BufferedVertex& v, float value, int extra_index) {
	v.setExtra(value, "extra0", extra_index);
}


std::function<void(float, float)> deformationOperator(const std::function<void(BufferedVertex&, float, float)>& deformation, IndexedMesh3D& mesh, const PolyGroupID& id) {
	return [&deformation, &mesh, id](float t, float delta) {
		mesh.deformPerVertex(id, [deformation, t, delta](BufferedVertex& v) {
			deformation(v, t, delta);
		});
	};
}

std::function<void(float)> deformationOperator(const std::function<void(BufferedVertex&, float)>& deformation, IndexedMesh3D& mesh, const PolyGroupID& id) {
	return [&deformation, &mesh, id](float t) {
		mesh.deformPerVertex(id, [deformation, t](BufferedVertex& v) {
			deformation(v, t);
		});
	};
}

std::function<void(float, float)> moveAlongCurve(const SmoothParametricCurve& curve, IndexedMesh3D& mesh, const PolyGroupID& id) {
	return deformationOperator([curve](BufferedVertex& v, float t, float delta) {
		v.setPosition(v.getPosition() + curve(t) - curve(t - delta));
	}, mesh, id);
}

// ReSharper disable once CppPassValueParameterByConstReference
mat3 IndexedMesh3D::inertiaTensorCMAppBd(PolyGroupID id) const {
	return inertiaTensorAppBd(id, centerOfMass(id));
}

vector<int> IndexedMesh3D::findVertexNeighbours(int i, const PolyGroupID& id) const {
	std::set<int> neighbours = {};
	for (const IndexedTriangle& t : triangles.at(polygroupIndexOrder.at(id)))
		if (contains(i, t.getVertexIndices())) {
			neighbours.insert(setMinus(t.getVertexIndices(), i)[0]);
			neighbours.insert(setMinus(t.getVertexIndices(), i)[1]);
		}
	return vector(neighbours.begin(), neighbours.end());
}

vector<int> IndexedMesh3D::findVertexParentTriangles(int i, const PolyGroupID& id) const {
	std::set<int> parentTriangles = {};
	for (int j = 0; j < triangles.at(polygroupIndexOrder.at(id)).size(); j++)
		if (contains(i, triangles.at(polygroupIndexOrder.at(id))[j].getVertexIndices()))
			parentTriangles.insert(j);
	return vector(parentTriangles.begin(), parentTriangles.end());
}

void IndexedMesh3D::recalculateNormal(int i, const PolyGroupID& id) {
	auto trs = findVertexParentTriangles(i, id);
	vec3 n = vec3(0);
	for (int j : trs)
		n += triangles.at(polygroupIndexOrder.at(id))[j].faceNormal() * triangles.at(polygroupIndexOrder.at(id))[j].area();
	vertices.at(polygroupIndexOrder.at(id))[i].setNormal(normalize(n));
}

void IndexedMesh3D::recalculateNormalsNearby(int i, const PolyGroupID& id) {
	for (int j : findVertexNeighbours(i, id))
		recalculateNormal(j, id);
}

void IndexedMesh3D::recalculateNormals(const PolyGroupID& id) {
	for (int i = 0; i < vertices.at(polygroupIndexOrder.at(id)).size(); i++)
		recalculateNormal(i, id);
}

void IndexedMesh3D::recalculateNormals() {
	for (auto id : getPolyGroupIDs())
		recalculateNormals(id);
}

void IndexedMesh3D::orientFaces(const PolyGroupID& id) {
	for (auto& t : triangles.at(polygroupIndexOrder.at(id)))
		if (dot(t.faceNormal(), t.getVertex(0).get_normal()) < 0)
			t.changeOrientation();
}

void IndexedMesh3D::orientFaces() {
	for (auto id : getPolyGroupIDs())
		orientFaces(id);
}

vector<int> IndexedMesh3D::findNeighboursSorted(int i, const PolyGroupID& id) const {
	vector<int> neighbours = findVertexNeighbours(i, id);
	std::map<int, float> angles = {};
	auto t = orthogonalComplementBasis(vertices.at(polygroupIndexOrder.at(id))[i].getNormal());
	for (int j = 0; j < neighbours.size(); j++)
		angles[i] = polarAngle(vertices.at(polygroupIndexOrder.at(id))[neighbours[j]].getPosition() - vertices.at(polygroupIndexOrder.at(id))[i].getPosition(),
							   vertices.at(polygroupIndexOrder.at(id))[i].getNormal());
	std::ranges::sort(neighbours, [&angles](int a, int b) {
		return angles[a] < angles[b];
	});
	return neighbours;
}

bool IndexedMesh3D::checkIfHasCompleteNeighbourhood(int i, const PolyGroupID& id) const {
	vector<int> trs = findVertexParentTriangles(i, id);
	vector<int> neighbours = findVertexNeighbours(i, id);
	for (int p : neighbours) {
		int found = 0;
		for (int t : trs)
			if (triangles.at(polygroupIndexOrder.at(id))[t].containsEdge(i, p))
				found++;
		if (found < 2)
			return false;
	}
	return true;
}

float IndexedMesh3D::meanCurvature(int i, const PolyGroupID& id) const {
	if (!checkIfHasCompleteNeighbourhood(i, id))
		return 0;
	float sum = 0;
	auto nbhd = findNeighboursSorted(i, id);
	vec3 p = vertices.at(polygroupIndexOrder.at(id))[i].getPosition();
	for (int j = 0; j < nbhd.size(); j++) {
		vec3 prev = vertices.at(polygroupIndexOrder.at(id))[nbhd[(j - 1 + nbhd.size()) % nbhd.size()]].getPosition();
		vec3 next = vertices.at(polygroupIndexOrder.at(id))[nbhd[(j + 1) % nbhd.size()]].getPosition();
		vec3 current = vertices.at(polygroupIndexOrder.at(id))[nbhd[j]].getPosition();
		float angle1 = abs(angle(prev - p, prev - current));
		float angle2 = abs(angle(next - p, next - current));
		sum += 0.5f * (cot(angle1) + cot(angle2)) * length(current - p);
	}
	return sum;
}

vec3 IndexedMesh3D::meanCurvatureVector(int i, const PolyGroupID& id) const {
	return meanCurvature(i, id) * vertices.at(polygroupIndexOrder.at(id))[i].getNormal();
}

float IndexedMesh3D::GaussCurvature(int i, const PolyGroupID& id) const {
	auto nbhd = findNeighboursSorted(i, id);
	float sum = 0;
	vec3 p = vertices.at(polygroupIndexOrder.at(id))[i].getPosition();
	for (int j = 0; j < nbhd.size(); j++) {
		vec3 prev = vertices.at(polygroupIndexOrder.at(id))[nbhd[(j - 1 + nbhd.size()) % nbhd.size()]].getPosition();
		vec3 next = vertices.at(polygroupIndexOrder.at(id))[nbhd[(j + 1) % nbhd.size()]].getPosition();
		vec3 current = vertices.at(polygroupIndexOrder.at(id))[nbhd[j]].getPosition();
		sum += angle(prev - current, next - current) / length(prev - next) * length(p - current) / 2;
	}
	return sum;
}


mat3 IndexedMesh3D::inertiaTensorAppBd(PolyGroupID id, vec3 p) const {
	vec3 cm = p;
	return integrateOverTriangles<mat3>([ cm](const IndexedTriangle& t) {
		vec3 p = t.center();
		return mat3(pow((p - cm).z, 2) + pow((p - cm).y, 2), -(p - cm).x * (p - cm).y, -(p - cm).x * (p - cm).z, -(p - cm).x * (p - cm).y, pow((p - cm).x, 2) + pow((p - cm).z, 2),
					-(p - cm).y * (p - cm).z, -(p - cm).x * (p - cm).z, -(p - cm).y * (p - cm).z, pow((p - cm).x, 2) + pow((p - cm).y, 2));
	}, id);
}


void PipeCurveVertexShader::transform(SpaceAutomorphism F) {
	deformPerVertex(id, [this, &F](BufferedVertex& v) {
		auto p = v.getPosition();
		mat3 J = F.df(p);
		auto b = J * getBinormal(v);
		auto n = J * v.getNormal();
		v.setPosition(F(p));
		v.setNormal(n);
		setBinormal(v, b);
	});
}

SurfacePlotDiscretisedMesh::SurfacePlotDiscretisedMesh(const DiscreteRealFunctionR2& plot) {
	vector<Vertex> points = vector<Vertex>();
	vector<ivec3> triangles = vector<ivec3>();

	auto ts = plot.args_t();
	auto xs = plot.args_x();

	for (int i = 0; i < plot.samples_t(); i++) {
		float t = ts[i];
		for (int j = 0; j < plot.samples_x(); j++) {
			float x = xs[j];
			vec3 p = vec3(x, t, plot[i][j]);
			vec3 n = e3;
			if (i > 0 && j > 0) {
				auto p1 = points.back().get_position();
				auto p2 = points.at(points.size() - plot.samples_x()).get_position();
				n = normalize(cross(p1 - p, p2 - p));
				if (dot(n, e3) < 0)
					n = -n;
				triangles.emplace_back((i - 1) * plot.samples_x() + j - 1, (i - 1) * plot.samples_x() + j, i * plot.samples_x() + j);
				triangles.emplace_back((i - 1) * plot.samples_x() + j - 1, i * plot.samples_x() + j, i * plot.samples_x() + j - 1);
			}
			points.emplace_back(p, vec2(x, t), n, vec4(0));
		}
	}
	addNewPolygroup(points, triangles, randomID());
}

void SurfacePlotDiscretisedMesh::transform(SpaceAutomorphism F) {
	deformPerVertex([F](BufferedVertex& v) {
		vec3 p = v.getPosition();
		mat3 J = F.df(p);
		vec3 n = J * v.getNormal();
		v.setPosition(F(p));
		v.setNormal(n);
	});
}

SurfacePolarPlotDiscretisedMesh::SurfacePolarPlotDiscretisedMesh(const DiscreteRealFunctionR2& plot, float r, float rot_speed) {
	DiscreteRealFunctionR2 f = plot;
	f.setDomain_x(vec2(-PI, PI));
	vector<Vertex> points = vector<Vertex>();
	vector<ivec3> triangles = vector<ivec3>();

	auto ts = f.args_t();
	auto xs = f.args_x();
	int nt = ts.size();
	int nx = xs.size();

	for (int i = 0; i < nt; i++) {
		float t = ts[i];
		for (int j = 0; j < nx; j++) {
			float phi = xs[j] + t * rot_speed;
			// float R = r-_f[i][j];
			float R = r - log(1 + f[i][j]);

			vec3 p = vec3(R * sin(phi), t, R * cos(phi));
			vec3 n = -vec3(sin(phi), 0, cos(phi));
			if (i > 0) {
				auto p1 = points.back().get_position();
				auto p2 = points.at(points.size() - nx).get_position();
				if (j == 0)
					p2 = points.at(points.size() - 2).get_position();
				n = normalize(cross(p1 - p, p2 - p));
				if (dot(n, vec3(sin(phi), 0, cos(phi))) < 0)
					n = -n;
				triangles.emplace_back((i - 1) * nx + (j - 1 + nx) % nx, (i - 1) * nx + j, i * nx + j);
				triangles.emplace_back((i - 1) * nx + (j - 1 + nx) % nx, i * nx + j, i * nx + (j - 1 + nx) % nx);
			}
			points.emplace_back(p, vec2(phi, t), n, vec4(0));
		}
	}
	addNewPolygroup(points, triangles, randomID());
}

FoliatedParametricSurfaceMesh::FoliatedParametricSurfaceMesh(const FoliatedParametricSurfaceMesh& other)
: IndexedMesh3D(other), foliation(other.foliation), special_leaves(other.special_leaves), continuous_leaves(other.continuous_leaves), leaf_radial_res(other.leaf_radial_res),
  leaf_hor_res(other.leaf_hor_res), color_map(other.color_map), special_leaf_colors(other.special_leaf_colors), leaf_radius_map(other.leaf_radius_map),
  special_leaf_radii(other.special_leaf_radii), id_prefix(other.id_prefix) {
	// boss = make_unique<BufferManager>(std::set({POSITION, NORMAL, UV, COLOR, INDEX, EXTRA0}));
}

FoliatedParametricSurfaceMesh::FoliatedParametricSurfaceMesh(FoliatedParametricSurfaceMesh&& other) noexcept
: IndexedMesh3D(std::move(other)), foliation(std::move(other.foliation)), special_leaves(other.special_leaves), continuous_leaves(other.continuous_leaves),
  leaf_radial_res(other.leaf_radial_res), leaf_hor_res(other.leaf_hor_res), color_map(std::move(other.color_map)), special_leaf_colors(std::move(other.special_leaf_colors)),
  leaf_radius_map(std::move(other.leaf_radius_map)), special_leaf_radii(std::move(other.special_leaf_radii)), id_prefix(std::move(other.id_prefix)) {
	// boss = make_unique<BufferManager>(*other.boss);
}

FoliatedParametricSurfaceMesh& FoliatedParametricSurfaceMesh::operator=(const FoliatedParametricSurfaceMesh& other) {
	if (this == &other)
		return *this;
	IndexedMesh3D::operator =(other);
	foliation = other.foliation;
	special_leaves = other.special_leaves;
	continuous_leaves = other.continuous_leaves;
	leaf_radial_res = other.leaf_radial_res;
	leaf_hor_res = other.leaf_hor_res;
	color_map = other.color_map;
	special_leaf_colors = other.special_leaf_colors;
	leaf_radius_map = other.leaf_radius_map;
	special_leaf_radii = other.special_leaf_radii;
	id_prefix = other.id_prefix;
	return *this;
}

FoliatedParametricSurfaceMesh& FoliatedParametricSurfaceMesh::operator=(FoliatedParametricSurfaceMesh&& other) noexcept {
	if (this == &other)
		return *this;
	IndexedMesh3D::operator =(std::move(other));
	foliation = std::move(other.foliation);
	special_leaves = other.special_leaves;
	continuous_leaves = other.continuous_leaves;
	leaf_radial_res = other.leaf_radial_res;
	leaf_hor_res = other.leaf_hor_res;
	color_map = std::move(other.color_map);
	special_leaf_colors = std::move(other.special_leaf_colors);
	leaf_radius_map = std::move(other.leaf_radius_map);
	special_leaf_radii = std::move(other.special_leaf_radii);
	id_prefix = std::move(other.id_prefix);
	return *this;
}

FoliatedParametricSurfaceMesh::FoliatedParametricSurfaceMesh(ParametricSurfaceFoliation foliation, int special_leaves, int continuous_leaves, int leaf_radial_res, int leaf_hor_res,
															 std::function<vec4(float)> color_map, const vector<vec4>& special_leaf_colors, HOM(float, float) leaf_radius_map,
															 const vector<float>& special_leaf_radii, std::function<vec4(float)> extra1_map, const vector<vec4>& special_extra1s)
: IndexedMesh3D(), foliation(foliation), special_leaves(special_leaves), continuous_leaves(continuous_leaves), leaf_radial_res(leaf_radial_res), leaf_hor_res(leaf_hor_res),
  color_map(color_map), special_leaf_colors(special_leaf_colors), leaf_radius_map(leaf_radius_map), special_leaf_radii(special_leaf_radii), extra1_map(extra1_map),
  special_extra1s(special_extra1s) {
	boss = make_unique<BufferManager>("extra0", "extra1");

	int i = 0;
	for (float t : linspace(foliation.getDomain()[0], foliation.getDomain()[1], continuous_leaves, false)) {
		float r = leaf_radius_map(t);
		auto c = foliation.getLeaf(t);
		auto leaf = PipeCurveVertexShader(c, r, leaf_hor_res, leaf_radial_res, id_prefix + "_" + std::to_string(i),
										  dict(string, vec4)({{"extra0", color_map(t)}, {"extra1", extra1_map(t)}}));
		mergeAndKeepID(leaf);
		i++;
	}
	for (int j = 0; j < special_leaves; j++) {
		float r = special_leaf_radii[j];
		auto c = foliation.getSpecialLeaf(j);
		auto leaf = PipeCurveVertexShader(c, r, leaf_hor_res, leaf_radial_res, id_prefix + "_s" + std::to_string(j),
										  dict(string, vec4)({{"extra0", special_leaf_colors[j]}, {"extra1", special_extra1s[j]}}));
		mergeAndKeepID(leaf);
	}
}
