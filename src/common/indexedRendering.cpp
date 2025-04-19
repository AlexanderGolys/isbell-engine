#include "indexedRendering.hpp"

#include <stdio.h>
#include <vector>
#include<fstream>
#include<sstream>

using namespace glm;
using std::vector, std::string, std::shared_ptr, std::unique_ptr, std::pair, std::make_unique, std::make_shared, std::array, std::weak_ptr;

BufferManager::BufferManager(const std::set<CommonBufferType> &activeBuffers) {
    stds = make_unique<Stds>();

    extra0 = activeBuffers.contains(EXTRA0) ? make_unique<BUFF4>() : nullptr;
    extra = nullptr;
    if (activeBuffers.contains(EXTRA1) || activeBuffers.contains(EXTRA2) || activeBuffers.contains(EXTRA3) || activeBuffers.contains(EXTRA4))
        extra = make_unique<buff4x4>();
    indices = activeBuffers.contains(INDEX) ? make_unique<IBUFF3>() : nullptr;
    this->activeBuffers = activeBuffers;
}



BufferManager::BufferManager(bool materials, const std::set<CommonBufferType> &extras) : BufferManager({POSITION, NORMAL, UV, COLOR, INDEX}) {
    if (materials) {
        this->activeBuffers = {POSITION, NORMAL, UV, COLOR, MATERIAL1, MATERIAL2, MATERIAL3, MATERIAL4, INDEX};
    }
    extra0 = extras.contains(EXTRA0) ? make_unique<BUFF4>() : nullptr;
    extra = nullptr;
    if (extras.contains(EXTRA1) || extras.contains(EXTRA2) || extras.contains(EXTRA3) || extras.contains(EXTRA4))
        extra = make_unique<buff4x4>();
    for (CommonBufferType type: extras)
        activeBuffers.insert(type);
}

int BufferManager::bufferLength(CommonBufferType type) const {
    if (!isActive(type))
        return 0;
    switch (type) {
        case POSITION:
        case NORMAL:
        case UV:
        case COLOR:
        case MATERIAL1:
        case MATERIAL2:
        case MATERIAL3:
        case MATERIAL4:
        case EXTRA0:
        case EXTRA1:
        case EXTRA2:
        case EXTRA3:
        case EXTRA4:
            return stds->positions.size();
        case INDEX:
            return indices->size();
    }
    throw UnknownVariantError("Buffer not recognised among common types. ");
}

void *BufferManager::firstElementAddress(CommonBufferType type) const {
    switch (type) {
    case POSITION:
        return &stds->positions[0];
    case NORMAL:
        return &stds->normals[0];
    case UV:
        return &stds->uvs[0];
    case COLOR:
        return &stds->colors[0];
    case EXTRA0:
        return extra0->data();
    case EXTRA1:
        return extra->a.data();
    case EXTRA2:
        return extra->b.data();
    case EXTRA3:
        return extra->c.data();
    case EXTRA4:
        return extra->d.data();
    case INDEX:
        return indices->data();
    default:
		throw UnknownVariantError("Buffer not recognised among common types. ");
    }
    throw UnknownVariantError("Buffer not recognised among common types. ");
}


int BufferManager::addTriangleVertexIndices(ivec3 ind, int shift) {
    indices->push_back(ind+ivec3(shift));
    return bufferLength(INDEX) - 1;
}

int BufferManager::addStdAttributesFromVertex(vec3 pos, vec3 norm, vec2 uv, vec4 col) {
    stds->positions.push_back(pos);
    stds->normals.push_back(norm);
    stds->uvs.push_back(uv);
    stds->colors.push_back(col);
	if (activeBuffers.contains(EXTRA0))
		extra0->push_back(vec4(0));
	if (activeBuffers.contains(EXTRA1))
		extra->a.push_back(vec4(0));
	if (activeBuffers.contains(EXTRA2))
		extra->b.push_back(vec4(0));
	if (activeBuffers.contains(EXTRA3))
		extra->c.push_back(vec4(0));
	if (activeBuffers.contains(EXTRA4))
		extra->d.push_back(vec4(0));
    return stds->positions.size() - 1;
}

int BufferManager::addAttributesFromVertex(vec3 pos, vec3 norm, vec2 uv, vec4 col, vec4 extra0Data) {
	stds->positions.push_back(pos);
	stds->normals.push_back(norm);
	stds->uvs.push_back(uv);
	stds->colors.push_back(col);
	extra0->push_back(extra0Data);
	if (activeBuffers.contains(EXTRA1))
		extra->a.push_back(vec4(0));
	if (activeBuffers.contains(EXTRA2))
		extra->b.push_back(vec4(0));
	if (activeBuffers.contains(EXTRA3))
		extra->c.push_back(vec4(0));
	if (activeBuffers.contains(EXTRA4))
		extra->d.push_back(vec4(0));
	return stds->positions.size() - 1;
}

int BufferManager::addFullVertexData(vec3 pos, vec3 norm, vec2 uv, vec4 col) {
    return addStdAttributesFromVertex(pos, norm, uv, col);
}

int BufferManager::addFullVertexData(vec3 pos, vec3 norm, vec2 uv, vec4 col, vec4 extra0Data) {
	return addAttributesFromVertex(pos, norm, uv, col, extra0Data);
}

int BufferManager::addFullVertexData(const Vertex &v) {
    return addStdAttributesFromVertex(v.getPosition(), v.getNormal(), v.getUV(), v.getColor());
}

int BufferManager::addFullVertexData(const Vertex &v, vec4 extra0Data) {
	return addAttributesFromVertex(v.getPosition(), v.getNormal(), v.getUV(), v.getColor(), extra0Data);
}

BufferedVertex::BufferedVertex(BufferManager &bufferBoss, const Vertex &v) : bufferBoss(bufferBoss) {
    index = this->bufferBoss.addFullVertexData(v);
}

void BufferedVertex::applyFunction(const SpaceEndomorphism &f) {
    vec3 p = getPosition();
    vec3 n = getNormal();
    setPosition(f(p));
    setNormal(normalise(f.df(p)*n));
}


mat3 IndexedTriangle::orthonormalFrame() const {
    vec3 p0 = getVertex(0).getPosition();
    vec3 p1 = getVertex(1).getPosition();
    vec3 p2 = getVertex(2).getPosition();
    return  GramSchmidtProcess(mat3(p0-p2, p1-p2, cross(p1-p2, p0-p2)));
}

vec3 IndexedTriangle::fromPlanar(vec2 v) const {
    mat3 frame = orthonormalFrame();
    return frame[0]*v.x + frame[1]*v.y + getVertex(2).getPosition();
}

vec2 IndexedTriangle::toPlanar(vec3 v) const {
    mat3 frame = orthonormalFrame();
    return vec2(dot(frame[0], v-getVertex(2).getPosition()), dot(frame[1], v-getVertex(2).getPosition()));
}


vec3 IndexedTriangle::fromBars(vec2 v) const {
    return getVertex(0).getPosition() * v.x + getVertex(1).getPosition() * v.y +
           getVertex(2).getPosition() * (1 - v.x - v.y);
}


std::array<vec3, 3> IndexedTriangle::borderTriangle(float width) const {
    vec3 p0 = getVertex(0).getPosition();
    vec3 p1 = getVertex(1).getPosition();
    vec3 p2 = getVertex(2).getPosition();
    vec3 n = normalize(cross(p1 - p0, p2 - p0));
    vec3 v01 = p1-p0;
    vec3 v12 = p2-p1;
    vec3 v20 = p0-p2;
    vec3 w01 = GramSchmidtProcess(mat3(n, v01, v12))[2];
    vec3 w12 = GramSchmidtProcess(mat3(n, v12, v20))[2];
    vec3 w20 = GramSchmidtProcess(mat3(n, v20, v01))[2];
    w01 = w01*sign(dot(p2-p0, w01));
    w12 = w12*sign(dot(p0-p1, w12));
    w20 = w20*sign(dot(p1-p2, w20));

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

bool IndexedTriangle::containsEdge(int i, int j) const {
	return contains(getVertexIndices(), i) && contains(getVertexIndices(), j);
}


WeakSuperMesh::WeakSuperMesh(WeakSuperMesh &&other) noexcept: boss(std::move(other.boss)),
															  vertices(std::move(other.vertices)),
															  triangles(std::move(other.triangles)),
															  material(std::move(other.material)) {}

WeakSuperMesh & WeakSuperMesh::operator=(WeakSuperMesh &&other) noexcept {
    if (this == &other)
        return *this;
    boss = std::move(other.boss);
    vertices = std::move(other.vertices);
    triangles = std::move(other.triangles);
    material = std::move(other.material);
    return *this;
}

WeakSuperMesh::WeakSuperMesh() {boss = make_unique<BufferManager>();}

WeakSuperMesh::WeakSuperMesh(const vector<Vertex> &hardVertices, const vector<ivec3> &faceIndices,const PolyGroupID &id) : WeakSuperMesh() {
    addNewPolygroup(hardVertices, faceIndices, id);
}

WeakSuperMesh::WeakSuperMesh(const char *filename, const PolyGroupID &id) : WeakSuperMesh() {
    addNewPolygroup(filename, id);
}

void WeakSuperMesh::addNewPolygroup(const char *filename, const PolyGroupID &id) {
	int shift = boss->bufferLength(POSITION);
    vertices[id] = vector<BufferedVertex>();
    triangles[id] = vector<IndexedTriangle>();
    vector<vec3> norms = {};
    vector<vec2> uvs = {};
    vector<vec3> pos = {};

    // auto estimatedSizes = countEstimatedBufferSizesInOBJFile(filename);
    // vertices[id].reserve(estimatedSizes["positions"]);
    // triangles[id].reserve(estimatedSizes["faces"]);
    // norms.reserve(estimatedSizes["normals"]);
    // uvs.reserve(estimatedSizes["uvs"]);



    std::ifstream in(filename, std::ios::in);
    if (!in)
    {
        std::cerr << "Cannot open " << filename << std::endl;
        exit(1);
    }

    string line;
    while (getline(in, line))
    {
        if (line.substr(0,2) == "v ")
        {
            std::istringstream s(line.substr(2));
            vec3 vert;
            s >> vert.x; s >> vert.y; s >> vert.z;
            vertices[id].emplace_back(*boss, Vertex(vert, vec2(2137, 0), vec3(0, 0, 0), BLACK));
            pos.push_back(vert);
        }
        if (line.substr(0,3) == "vt ")
        {
            std::istringstream s(line.substr(3));
            vec2 uv;
            s >> uv.x; s >> uv.y;
            uvs.push_back(uv);
        }

        if (line.substr(0,3) == "vn ")
        {
            std::istringstream s(line.substr(3));
            vec3 norm;
            s >> norm.x; s >> norm.y; s >> norm.z;
            norms.push_back(norm);
        }

        else if (line.substr(0,2) == "f ")
        {
            ivec3 vertexIndex, uvIndex, normalIndex;
            sscanf(line.c_str(),"f %d/%d/%d %d/%d/%d %d/%d/%d", &vertexIndex[0], &uvIndex[0], &normalIndex[0], &vertexIndex[1], &uvIndex[1], &normalIndex[1], &vertexIndex[2], &uvIndex[2], &normalIndex[2]);

            vertexIndex -= ivec3(1, 1, 1);
            uvIndex -= ivec3(1, 1, 1);
            normalIndex -= ivec3(1, 1, 1);


            vector uvsVec = { uvs.at(uvIndex[0]), uvs.at(uvIndex[1]), uvs.at(uvIndex[2]) };
            vector normsVec = { norms.at(normalIndex[0]), norms.at(normalIndex[1]), norms.at(normalIndex[2]) };

            for (int i = 0; i < 3; i++) {
                if (vertices[id][vertexIndex[i]].getNormal() == vec3(0, 0, 0)) {
                    vertices[id][vertexIndex[i]].setNormal(normsVec[i]);
                }
                if (vertices[id][vertexIndex[i]].getUV() == vec2(2137, 0)) {
                    vertices[id][vertexIndex[i]].setUV(uvsVec[i]);
                }
            }
            ivec3 face = ivec3(0);
            for (int i = 0; i < 3; i++)
				if (vertices[id][vertexIndex[i]].getNormal() == normsVec[i] && vertices[id][vertexIndex[i]].getUV() == uvsVec[i])
					face[i] = vertexIndex[i];
				else {
					vertices[id].emplace_back(*boss, Vertex(vertices[id][vertexIndex[i]].getPosition(), uvsVec[i], normsVec[i], BLACK));
					face[i] = vertices[id].size() - 1;
				}

			triangles[id].emplace_back(*boss, face, shift);
        }
    }
}


WeakSuperMesh & WeakSuperMesh::operator=(const WeakSuperMesh &other) {
    if (this == &other)
        return *this;
    boss = std::make_unique<BufferManager>(*other.boss);
    vertices = other.vertices;
    triangles = other.triangles;
    material = other.material;
    return *this;
}

WeakSuperMesh::WeakSuperMesh(const SmoothParametricSurface &surf, int tRes, int uRes, const PolyGroupID &id)
: WeakSuperMesh() {
	addUniformSurface(surf, tRes, uRes, id);
}


void WeakSuperMesh::addNewPolygroup(const vector<Vertex> &hardVertices, const vector<ivec3> &faceIndices, const PolyGroupID &id) {

    if (vertices.contains(id) || triangles.contains(id))
        throw IllegalVariantError("Polygroup ID already exists in mesh. ");

	int shift = boss->bufferLength(POSITION);
    vertices[id] = vector<BufferedVertex>();
    triangles[id] = vector<IndexedTriangle>();
    vertices[id].reserve(hardVertices.size());
    triangles[id].reserve(faceIndices.size());

    for (const Vertex &v: hardVertices)
		vertices[id].emplace_back(*boss, v);

	for (const ivec3 &ind: faceIndices)
        triangles[id].emplace_back(*boss, ind, shift);
}

void WeakSuperMesh::addNewPolygroup(const std::vector<Vertex> &hardVertices, const std::vector<ivec3> &faceIndices, const PolyGroupID &id, const std::vector<vec4> &extra0) {
	addNewPolygroup(hardVertices, faceIndices, id);
	for (int i = 0; i < vertices[id].size(); i++)
		boss->setExtra(vertices[id][i].getIndex(), extra0[i], 0);
}

void WeakSuperMesh::addNewPolygroup(const vector<Vertex> &hardVertices, const vector<ivec3> &faceIndices, const PolyGroupID &id, const vector<vec4> &extra0, const vector<mat4> &extra) {
	addNewPolygroup(hardVertices, faceIndices, id);
	for (int i = 0; i < vertices[id].size(); i++)
		boss->setExtra(vertices[id][i].getIndex(), extra0[i], 0);
	for (int e = 0; e < extra.size(); e++)
		for (int i = 0; i < vertices[id].size(); i++)
			boss->setExtra(vertices[id][i].getIndex(), extra[e][i], e+1);
}

WeakSuperMesh::WeakSuperMesh(const WeakSuperMesh &other): boss(&*other.boss),
														  vertices(other.vertices),
														  triangles(other.triangles),
														  material(other.material) {}


void WeakSuperMesh::addUniformSurface(const SmoothParametricSurface &surf, int tRes, int uRes, const PolyGroupID &id) {
    vector<Vertex> hardVertices = {};
    vector<ivec3> faceIndices = {};
    hardVertices.reserve(tRes * uRes + uRes + tRes + 1);
    faceIndices.reserve(2 * (tRes) * (uRes));



    for (int i = 0; i < tRes; i++)
        for (int j = 0; j < uRes; j++) {

            float t_ = 1.f * i / (1.f * tRes - 1);
        	if (surf.isPeriodicT()) t_ = 1.f * i / (1.f * tRes);

            float u_ = 1.f * j / (1.f * uRes - 1);
        	if (surf.isPeriodicU()) u_ = 1.f * j / (1.f * uRes);

            float t = lerp(surf.tMin(), surf.tMax(), t_);
            float u = lerp(surf.uMin(), surf.uMax(), u_);
            hardVertices.emplace_back(surf(t, u), vec2(t_, u_), surf.normal(t, u), vec4(t, u, 0, 1));
        }

	ivec2 size = ivec2(tRes, uRes);

    for (int i = 1; i < tRes ; i++)
        for (int j = 1; j < uRes; j++) {

        	int i0 = flattened2DVectorIndex(i, j, size);
            int i1 = flattened2DVectorIndex(i-1, j, size);
            int i2 = flattened2DVectorIndex(i-1, j-1, size);
        	int i3 = flattened2DVectorIndex(i, j-1, size);

			faceIndices.emplace_back(i0, i1, i2);
			faceIndices.emplace_back(i0, i3, i2);
        }


	if (surf.isPeriodicT())
		for (int j = 1; j < uRes; j++) {
			int i0 = flattened2DVectorIndex(0, j, size);
			int i1 = flattened2DVectorIndex(-1, j, size);
			int i2 = flattened2DVectorIndex(-1, j-1, size);
			int i3 = flattened2DVectorIndex(0, j-1, size);
			faceIndices.emplace_back(i0, i1, i2);
			faceIndices.emplace_back(i0, i3, i2);
		}
//
	if (surf.isPeriodicU())
		for (int i = 1; i < tRes; i++) {
			int i0 = flattened2DVectorIndex(i, 0, size);
			int i1 = flattened2DVectorIndex(i-1, 0, size);
			int i2 = flattened2DVectorIndex(i-1, -1, size);
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


void WeakSuperMesh::merge(const WeakSuperMesh &other) {
	for (const auto& id: other.getPolyGroupIDs())
		addNewPolygroup(other.getVertices(id), other.getIndices(id), make_unique_id(id));
}
void WeakSuperMesh::mergeAndKeepID(const WeakSuperMesh &other) {
	for (const auto& id: other.getPolyGroupIDs())
		addNewPolygroup(other.getVertices(id), other.getIndices(id), id);
}

void WeakSuperMesh::copyPolygroup(const WeakSuperMesh &other, const PolyGroupID &id, const PolyGroupID &newId) {
	vector<Vertex> vertices = other.getVertices(id);
	vector<ivec3> indices = other.getIndices(id);
	addNewPolygroup(vertices, indices, newId);}

void WeakSuperMesh::copyPolygroup(const PolyGroupID &id, PolyGroupID newId) {
	vector<Vertex> vertices = getVertices(id);
	vector<ivec3> indices = getIndices(id);
	addNewPolygroup(vertices, indices, newId);
}

bool WeakSuperMesh::hasExtra0() const {
	return boss->hasExtra0();
}

bool WeakSuperMesh::hasExtra1() const {
	return boss->hasExtra1();
}

bool WeakSuperMesh::hasExtra2() const {
	return boss->hasExtra2();
}

bool WeakSuperMesh::hasExtra3() const {
	return boss->hasExtra3();
}

bool WeakSuperMesh::hasExtra4() const {
	return boss->hasExtra4();
}

bool WeakSuperMesh::hasExtra(int slot) const {
	if (slot == 0) return hasExtra0();
	if (slot == 1) return hasExtra1();
	if (slot == 2) return hasExtra2();
	if (slot == 3) return hasExtra3();
	if (slot == 4) return hasExtra4();
	throw std::invalid_argument("slot must be between 0 and 4");
}

void * WeakSuperMesh::bufferIndexLocation() const {
	return boss->firstElementAddress(INDEX);
}

size_t WeakSuperMesh::bufferIndexSize() const {
	return boss->bufferSize(INDEX);
}

int WeakSuperMesh::bufferIndexLength() const {
	return boss->bufferLength(INDEX);
}

const void * WeakSuperMesh::getBufferLocation(CommonBufferType type) const {
	return boss->firstElementAddress(type);
}

unsigned int WeakSuperMesh::getBufferLength(CommonBufferType type) const {
	return boss->bufferLength(type);
}

size_t WeakSuperMesh::getBufferSize(CommonBufferType type) const {
	return boss->bufferSize(type);
}


std::vector<PolyGroupID> WeakSuperMesh::getPolyGroupIDs() const {
    vector<PolyGroupID> ids = {};
    for (const pair<PolyGroupID, vector<BufferedVertex>> &p: vertices)
        ids.push_back(p.first);
    return ids;
}

BufferManager & WeakSuperMesh::getBufferBoss() const {
	return *boss;
}

bool WeakSuperMesh::isActive(CommonBufferType type) const {
	return boss->isActive(type);
}

bool WeakSuperMesh::hasGlobalTextures() const {
	return !isActive(MATERIAL1) && material->textured();
}

BufferedVertex & WeakSuperMesh::getAnyVertexFromPolyGroup(const PolyGroupID &id) {
	return vertices.at(id).front();
}

void WeakSuperMesh::deformPerVertex(const PolyGroupID &id, const HOM(BufferedVertex&, void) &deformation) {
	for (auto &v : vertices.at(id))
		deformation(v);
}

void WeakSuperMesh::deformPerVertex(const HOM(BufferedVertex&, void) &deformation) {
	for (auto id: getPolyGroupIDs())
		for (auto &v : vertices.at(id))
			deformation(v);
}

void WeakSuperMesh::deformPerVertex(const PolyGroupID &id, const BIHOM(int, BufferedVertex&, void) &deformation) {
	for (int i=0; i<vertices.at(id).size(); i++)
		deformation(i, vertices.at(id)[i]);
}

void WeakSuperMesh::deformPerId(const BIHOM(BufferedVertex&, PolyGroupID, void) &deformation) {
	for (auto id: getPolyGroupIDs())
		for (auto &v : vertices.at(id))
			deformation(v, id);
}

vec2 WeakSuperMesh::getSurfaceParameters(const BufferedVertex &v) {
	return vec2(v.getColor().x, v.getColor().y);
}


void WeakSuperMesh::encodeSurfacePoint(BufferedVertex &v, const SmoothParametricSurface &surf, vec2 tu) {
	v.setPosition(surf(tu));
	v.setNormal(surf.normal(tu));
	v.setUV(vec2((tu.x-surf.tMin())/(surf.tMax()-surf.tMin()), (tu.y-surf.uMin())/(surf.uMax()-surf.uMin())));
	v.setColor(tu.x, 0);
	v.setColor(tu.y, 1);
}

void WeakSuperMesh::adjustToNewSurface(const SmoothParametricSurface &surf, const PolyGroupID &id) {
	deformPerVertex(id, [&](BufferedVertex &v) {
		encodeSurfacePoint(v, surf, getSurfaceParameters(v));
	});
}

void WeakSuperMesh::adjustToNewSurface(const SmoothParametricSurface &surf) {
	for (auto id: getPolyGroupIDs())
		adjustToNewSurface(surf, id);
}

void WeakSuperMesh::moveAlongVectorField(const PolyGroupID &id, VectorField X, float delta) {
    for (BufferedVertex &v: vertices.at(id))
        v.setPosition(X.moveAlong(v.getPosition(), delta));
}

void WeakSuperMesh::deformWithAmbientMap(const PolyGroupID &id, SpaceEndomorphism f) {
    for (BufferedVertex &v: vertices.at(id))
        v.applyFunction(f);
}

void WeakSuperMesh::deformWithAmbientMap(const SpaceEndomorphism &f) {
	for (auto id: getPolyGroupIDs())
		deformWithAmbientMap(id, f);
}

void WeakSuperMesh::initGlobalTextures() const {
	if (hasGlobalTextures())
		material->initTextures();
}

void WeakSuperMesh::affineTransform(const mat3 &M, vec3 v, const PolyGroupID &id) {
	deformWithAmbientMap(id, SpaceEndomorphism::affine(M, v));
}

void WeakSuperMesh::affineTransform(const mat3 &M, vec3 v) {
	for (auto& name: getPolyGroupIDs())
		affineTransform(M, v, name);
}

void WeakSuperMesh::shift(vec3 v, const PolyGroupID &id) {
	affineTransform(mat3(1), v, id);
}

void WeakSuperMesh::shift(vec3 v) {
	affineTransform(mat3(1), v);
}

void WeakSuperMesh::scale(float s, const PolyGroupID &id) {
	affineTransform(mat3(s), vec3(0), id);
}

void WeakSuperMesh::scale(float s) {
	affineTransform(mat3(s), vec3(0));
}

vector<Vertex> WeakSuperMesh::getVertices(const PolyGroupID &id) const {
    vector<Vertex> verts = {};
    verts.reserve(vertices.at(id).size());
    for (const BufferedVertex &v: vertices.at(id))
        verts.push_back(v.getVertex());
    return verts;
}

vector<BufferedVertex> WeakSuperMesh::getBufferedVertices(const PolyGroupID &id) const {
	return vertices.at(id);
}


vector<ivec3> WeakSuperMesh::getIndices(const PolyGroupID &id) const {
    vector<ivec3> inds = {};
    inds.reserve(triangles.at(id).size());
    for (const IndexedTriangle &t: triangles.at(id))
        inds.push_back(t.getVertexIndices());
    return inds;
}

vector<IndexedTriangle> WeakSuperMesh::getTriangles(const PolyGroupID &id) const {
	return triangles.at(id);
}

vec4 WeakSuperMesh::getIntencities() const {
	return material->compressIntencities();
}

MaterialPhong WeakSuperMesh::getMaterial() const {
	return *material;
}


void WeakSuperMesh::addGlobalMaterial(const MaterialPhong &mat) {
	material = std::make_shared<MaterialPhong>(mat);
}

void WeakSuperMesh::flipNormals(const PolyGroupID &id) {
	deformPerVertex(id, [](BufferedVertex &v) {
		v.setNormal(-v.getNormal());
	});
}

void WeakSuperMesh::pointNormalsInDirection(vec3 dir, const PolyGroupID &id) {
	deformPerVertex(id, [dir](BufferedVertex &v) {
	vec3 n = v.getNormal();
	v.setNormal(n*1.f*sgn(dot(n, dir)));
}); }

void WeakSuperMesh::pointNormalsInDirection(vec3 dir) {
	for (auto &name: getPolyGroupIDs())
		pointNormalsInDirection(dir, name);
}

WeakSuperMesh WeakSuperMesh::subdivideBarycentric(const PolyGroupID &id) const {
    vector<ivec3> trInds = getIndices(id);
    vector<Vertex> verts = getVertices(id);
    verts.reserve(verts.size() + trInds.size());
    vector<ivec3> newInds = {};
    newInds.reserve(3 * trInds.size());
    for (ivec3 tr: trInds) {
        Vertex b = barycenter(verts[tr.x], verts[tr.y], verts[tr.z]);
        verts.push_back(b);
        int centerInd = verts.size() - 1;
        newInds.push_back(ivec3(tr.x, tr.y, centerInd));
        newInds.push_back(ivec3(tr.y, tr.z, centerInd));
        newInds.emplace_back(tr.z, tr.x, centerInd);
    }
    return WeakSuperMesh(verts, newInds, id);
}


WeakSuperMesh WeakSuperMesh::subdivideEdgecentric(const PolyGroupID &id) const {
    vector<ivec3> trInds = getIndices(id);
    vector<Vertex> verts = getVertices(id);
    verts.reserve(verts.size() + 3 * trInds.size());
    vector<ivec3> newInds = {};
    newInds.reserve(4 * trInds.size());
    for (ivec3 tr: trInds) {
        verts.push_back(center(verts[tr.x], verts[tr.y]));
        int center1 = verts.size() - 1;
        verts.push_back(center(verts[tr.x], verts[tr.z]));
        int center2 = verts.size() - 1;
        verts.push_back(center(verts[tr.y], verts[tr.z]));
        int center3 = verts.size() - 1;
        newInds.push_back(ivec3(tr.x, center1, center2));
        newInds.push_back(ivec3(tr.y, center1, center3));
        newInds.push_back(ivec3(tr.z, center2, center3));
        newInds.push_back(ivec3(center1, center2, center3));
    }
    return WeakSuperMesh(verts, newInds, id);
}

WeakSuperMesh WeakSuperMesh::wireframe(PolyGroupID id, PolyGroupID targetId, float width, float heightCenter, float heightSide) const {
    vector<IndexedTriangle> trs = triangles.at(id);
    vector<BufferedVertex> verts = vertices.at(id);
    vector<Vertex> new_verts = getVertices(id);
    vector<ivec3> new_inds = {};
    for (const auto &v: verts) {
        Vertex extr = v.getVertex();
        extr.setPosition(extr.getPosition() + extr.getNormal() * heightCenter);
        new_verts.push_back(extr);
    }
    for (const auto &tr : trs) {
        array<vec3, 3> border = tr.borderTriangle(width);
        vec3 n = tr.faceNormal()*sign(dot(tr.getVertex(0).getNormal(), tr.faceNormal()));


        for (vec3 p : border)
            new_verts.emplace_back(p+n*heightSide, vec2(0, 0), n, BLACK, std::nullopt);
        for (vec3 p : border)
            new_verts.emplace_back(p, vec2(0, 0), n, BLACK, std::nullopt);
        new_inds.push_back(ivec3(tr.getVertexIndices().x, tr.getVertexIndices().y, tr.getVertexIndices().x + verts.size()));
        new_inds.push_back(ivec3(tr.getVertexIndices().y, tr.getVertexIndices().y + verts.size(), tr.getVertexIndices().x + verts.size()));
        new_inds.push_back(ivec3(tr.getVertexIndices().x+verts.size(), tr.getVertexIndices().y+verts.size(), new_verts.size()-6));
        new_inds.push_back(ivec3(tr.getVertexIndices().y+verts.size(), new_verts.size()-5, new_verts.size()-6));
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
    return WeakSuperMesh(new_verts, new_inds, targetId);
}


void BufferManager::insertValueToSingleBuffer(CommonBufferType type, void *valueAddress) {
    if (type == INDEX) {
        ivec3* value = static_cast<ivec3*>(valueAddress);
        indices->push_back(*value);
    }
   if (bufferElementLength(type) == 2) {
        vec2* v = static_cast<vec2*>(valueAddress);
        BUFF2* b = static_cast<BUFF2*>(firstElementAddress(type));
        b->push_back(*v);
        return;
    }
   if (bufferElementLength(type) == 3) {
        vec3* v = static_cast<vec3*>(valueAddress);
        BUFF3* b = static_cast<BUFF3*>(firstElementAddress(type));
        b->push_back(*v);
        return;
    }
   if (bufferElementLength(type) == 4) {
        vec4* v = static_cast<vec4*>(valueAddress);
        BUFF4* b = static_cast<BUFF4*>(firstElementAddress(type));
        b->push_back(*v);
        return;
    }
    throw IllegalVariantError("Buffer element has illegal length.");
}

void BufferManager::insertDefaultValueToSingleBuffer(CommonBufferType type) {
    if (type == INDEX) {
        indices->emplace_back(0, 0, 0);
        return;
    }
    if (bufferElementLength(type) == 2) {
        auto *b = static_cast<BUFF2 *>(firstElementAddress(type));
        b->emplace_back(0, 0);
        return;
    }
    if (bufferElementLength(type) == 3) {
        BUFF3 *b = static_cast<BUFF3 *>(firstElementAddress(type));
        b->emplace_back(0, 0, 0);
        return;
    }
    if (bufferElementLength(type) == 4) {
        BUFF4 *b = static_cast<BUFF4 *>(firstElementAddress(type));
        b->emplace_back(0, 0, 0, 0);
        return;
    }
    throw IllegalVariantError("Buffer element has illegal length.");
}

BufferManager::BufferManager(const BufferManager &other) :
    stds(std::make_unique<Stds>(*other.stds)),
    extra0(std::make_unique<BUFF4>(*other.extra0)),
    extra(std::make_unique<buff4x4>(*other.extra)),
    indices(std::make_unique<IBUFF3>(*other.indices)),
    activeBuffers(other.activeBuffers) {}

BufferManager & BufferManager::operator=(BufferManager &&other) noexcept {
    if (this == &other)
        return *this;
    stds = std::move(other.stds);
    extra0 = std::move(other.extra0);
    extra = std::move(other.extra);
    indices = std::move(other.indices);
    activeBuffers = std::move(other.activeBuffers);
    return *this;
}


BufferManager::BufferManager(BufferManager &&other) noexcept: stds(std::move(other.stds)),
                                                              extra0(std::move(other.extra0)),
                                                              extra(std::move(other.extra)),
                                                              indices(std::move(other.indices)),
                                                              activeBuffers(std::move(other.activeBuffers)) {}

BufferManager & BufferManager::operator=(const BufferManager &other) {
    if (this == &other)
        return *this;
    stds = std::make_unique<Stds>(*other.stds);
    extra0 = std::make_unique<BUFF4>(*other.extra0);
    extra = std::make_unique<buff4x4>(*other.extra);
    indices = std::make_unique<IBUFF3>(*other.indices);
    activeBuffers = other.activeBuffers;
    return *this;
}

void BufferManager::initialiseExtraBufferSlot(int slot) {
    switch (slot) {
        case 0:
            if (extra0 == nullptr)
                extra0 = make_unique<BUFF4>();
            activeBuffers.insert(EXTRA0);
            return;
    case 1:
        if (extra == nullptr)
            extra = make_unique<buff4x4>();
        activeBuffers.insert(EXTRA1);
        return;
    case 2:
        if (extra == nullptr)
            extra = make_unique<buff4x4>();
        activeBuffers.insert(EXTRA2);
        return;
        case 3:
            if (extra == nullptr)
                extra = make_unique<buff4x4>();
        activeBuffers.insert(EXTRA3);
        return;
        case 4:
            if (extra == nullptr)
                extra = make_unique<buff4x4>();
        activeBuffers.insert(EXTRA4);
		default: ;
	}
    }

vec4 BufferManager::getExtra(int index, int slot) const {
    switch (slot) {
    case 0:
        return (*extra0)[index];
    case 1:
        return  extra->a[index];
    case 2:
        return  extra->b[index];
    case 3:
        return  extra->c[index];
    case 4:
        return  extra->d[index];
    }
    throw UnknownVariantError("Extra slot not recognised. ");
}


void BufferManager::setExtra(int index, vec4 value, int slot) {
    switch (slot) {
        case 0:
            (*extra0)[index] = value;
            return;
        case 1:
            extra->a[index] = value;
            return;
        case 2:
            extra->b[index] = value;
            return;
        case 3:
            extra->c[index] = value;
            return;
        case 4:
            extra->d[index] = value;
            return;
    }
    throw UnknownVariantError("Extra slot not recognised. ");
}
void BufferManager::setExtra(int index, vec3 value, int slot) {
    setExtra(index, value.x, slot, 0);
    setExtra(index, value.y, slot, 1);
    setExtra(index, value.z, slot, 2);
}
void BufferManager::setExtra(int index, float value, int slot, int component) {
    switch (slot) {
        case 0:
            (*extra0)[index][component] = value;
            return;
        case 1:
            extra->a[index][component] = value;
            return;
        case 2:
            extra->b[index][component] = value;
            return;
        case 3:
            extra->c[index][component] = value;
            return;
        case 4:
            extra->d[index][component] = value;
            return;
        default:
            throw UnknownVariantError("Extra slot not recognised. ");

    }
}

void BufferManager::reserveSpace(int targetSize) {
    stds->positions.reserve(targetSize);
    stds->normals.reserve(targetSize);
    stds->uvs.reserve(targetSize);
    stds->colors.reserve(targetSize);

    if (isActive(EXTRA0))
        extra0->reserve(targetSize);
    if (isActive(EXTRA1) || isActive(EXTRA2) || isActive(EXTRA3) || isActive(EXTRA4)) {
        extra->a.reserve(targetSize);
        extra->b.reserve(targetSize);
        extra->c.reserve(targetSize);
        extra->d.reserve(targetSize);
    }

}

vec3 WeakSuperMesh::centerOfMass(PolyGroupID id) const {
	vec3 sum        = vec3(0);
	float totalArea = 0;
	for (const IndexedTriangle &t: triangles.at(id)) {
		sum += t.center()*t.area();
		totalArea += t.area();
	}
	return sum/totalArea;
}

vec3 WeakSuperMesh::centerOfMass() const {
	vec3 sum        = vec3(0);
	float totalArea = 0;
	for (const auto &id: getPolyGroupIDs())
		for (const IndexedTriangle &t: triangles.at(id)) {
			sum += t.center()*t.area();
			totalArea += t.area();
		}
	return sum/totalArea;
}


Wireframe::Wireframe(const SmoothParametricSurface &surf, float width, int n, int m, int curve_res_rad, int curve_res_hor)
	:	WeakSuperMesh(), width(width), surf(surf) , n(n), m(m), curve_res_rad(curve_res_rad), curve_res_hor(curve_res_hor)
{
	for (float t_i: linspace(surf.tMin(), surf.tMax(), n)) {
		SmoothParametricCurve curve = surf.constT(t_i);
		auto id = randomID();
			addUniformSurface(curve.pipe(width), curve_res_rad, curve_res_hor, id);
			deformPerVertex(id, [t_i](BufferedVertex &v) {
				v.setColor(t_i, 2);
				v.setColor(v.getColor().x, 3);
		});
	}

	for (float u_i: linspace(surf.uMin(), surf.uMax(), m)) {
		SmoothParametricCurve curve = surf.constU(u_i);
		auto id                     = randomID();
			addUniformSurface(curve.pipe(width), curve_res_rad, curve_res_hor, id);
			deformPerVertex(id, [u_i](BufferedVertex &v) {
				v.setColor(u_i, 3);
				v.setColor(v.getColor().x, 2);
		});
	}
}



void Wireframe::changeBaseSurface(const SmoothParametricSurface &newsurf) {
	for (auto id: getPolyGroupIDs())
		for (BufferedVertex &v: vertices.at(id)) {
			vec2 tu     = getSurfaceParameters(v);
			vec3 old_p0 = surf(tu);
			vec3 new_p0 = newsurf(tu);

			vec3 old_n0 = surf.normal(tu);
			vec3 new_n0 = newsurf.normal(tu);
			v.applyFunction(SpaceEndomorphism::affine(mat3(1), new_p0 - old_p0));
			v.applyFunction(SpaceAutomorphism::deltaRotation(old_n0, new_n0, v.getPosition()));
		}
	surf = newsurf;
}
vec2 Wireframe::getSurfaceParameters(const BufferedVertex &v) const {
	return vec2(v.getColor().z, v.getColor().w);
}

PlanarFlowLines::PlanarFlowLines(const VectorFieldR2 &X, float dt, int steps, const std::function<float(float, float, float, vec2, vec2)> &width, const std::function<vec4(float, float, float, vec2, vec2)> &color)
	: WeakSuperMesh(), X(X), dt(dt), steps(steps), width(width), color(color) {
	boss = make_unique<BufferManager>(std::set({POSITION, NORMAL, UV, COLOR, INDEX, EXTRA0}));
}


void PlanarFlowLines::generateGrid(vec2 v_min, vec2 v_max, ivec2 res) {
	startPoints = flatten(linspace2D(v_min, vec2((v_max-v_min).x/res.x, 0), vec2(0, (v_max-v_min).y/res.y), res.x+1, res.y+1));
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
	vector<vec4> extra0 = {};
	float len = 0;

	for (int j = 1; j <= steps; j++) {
		len += norm(X(p)*dt);
		p += X(p)*dt;
		vec2 n = X.normal(p);
		float t = t0 + j*dt;
		float speed = X.speed(p);
		float w = width(t, t0, speed, p, p0);
		vec2 q0 = p - n*w;
		vec2 q1 = p + n*w;
		vec4 c = color(t, t0, speed, p, p0);
		c.w = speed;
		verts.emplace_back(vec3(q0, t/(steps*dt*20)), vec2(t, w), e3, c);
		verts.emplace_back(vec3(p, t/(steps*dt*20)), vec2(t, 0), e3, c);
		verts.emplace_back(vec3(q1, t/(steps*dt*20)), vec2(t, -w), e3, c);

		extra0.emplace_back(t0, p0.x, p0.y, len);
		extra0.emplace_back(t0, p0.x, p0.y, len);
		extra0.emplace_back(t0, p0.x, p0.y, len);


		if (j > 1) {
			int last_ind = verts.size() - 1;
			trs.emplace_back(last_ind, last_ind-3, last_ind-4);
			trs.emplace_back(last_ind, last_ind-1, last_ind-4);
			trs.emplace_back(last_ind-1, last_ind-5, last_ind-4);
			trs.emplace_back(last_ind-1, last_ind-2, last_ind-5);
		}
	}

	addNewPolygroup(verts, trs, ids[i], extra0);
}

void PlanarFlowLines::generateLines() {
	for (int i = 0; i < startPoints.size(); i++)
		generateLine(i);
}

float PlanarFlowLines::getTimeRelative(const BufferedVertex &v) {
	return v.getUV().x;
}

float PlanarFlowLines::getT0(const BufferedVertex &v) {
	return v.getExtra(0).x;
}

float PlanarFlowLines::getTimeAbsolute(const BufferedVertex &v) {
	return getT0(v) + getTimeRelative(v);
}

vec2 PlanarFlowLines::getPos(const BufferedVertex &v) {
	return vec2(v.getPosition());
}

vec2 PlanarFlowLines::getStartPoint(const BufferedVertex &v) {
	return vec2(v.getExtra(0).y, v.getExtra(0).z);
}

float PlanarFlowLines::getSpeed(const BufferedVertex &v) {
	return v.getColor().w;
}

float PlanarFlowLines::getLength(const BufferedVertex &v) {
	return v.getExtra(0).w;
}

vec4 PlanarFlowLines::getColor(const BufferedVertex &v) {
	return vec4(vec3(v.getColor()), 1);
}

float PlanarFlowLines::getWidth(const BufferedVertex &v) {
	return v.getUV().y;
}

PlanarDiffusedInterval::	PlanarDiffusedInterval(const VectorFieldR2 &X, float dt, int steps, vec4 color, const HOM(float, float) &time_dump, const HOM(float, float) &width_ratio_dump, vec2 a, vec2 b, float t0)
: X(X),
  dt(dt),
  steps(steps),
  color(color),
  width_ratio_dump(width_ratio_dump),
  time_dump(time_dump),
  a(a),
  b(b),
  t0(t0)
{
	vec2 p_a = a;
	vec2 p_b = b;
	vector<Vertex> verts = {};
	vector<ivec3> trs = {};
	float t = t0;
	float eps = 0.1;
	float width = length(b-a);
	for (int j = 0; j <= steps; j++) {
		float width_ratio = width/length(p_b-p_a);
		float t_dump = time_dump(t-t0);
		float w_dump = width_ratio_dump(width_ratio);
		vec4 col = color;
		col.w = t_dump*w_dump;
		if (t_dump*w_dump < .05) col.w=0.f;
		verts.emplace_back(vec3(p_a, (t-t0)*eps), vec2(t, t0), e3, col);
		verts.emplace_back(vec3(p_b, (t-t0)*eps), vec2(t, t0), e3, col);
		p_a += X(p_a)*dt;
		p_b += X(p_b)*dt;
		t += dt;
		if (j > 0 && width_ratio > .15) {
			int last_ind = verts.size() - 1;
			trs.emplace_back(last_ind, last_ind-3, last_ind-2);
			trs.emplace_back(last_ind, last_ind-1, last_ind-3);
		}
	}
	addNewPolygroup(verts, trs, id);
}

PlanarDiffusedCurve::PlanarDiffusedCurve(const VectorFieldR2 &X, float dt, int steps, vec4 color, const HOM(float, float) &time_dump, const HOM(float, float) &width_ratio_dump, const SmoothParametricPlaneCurve &curve, vec2 domain, int resolution, float t0)
: X(X),
  dt(dt),
  steps(steps),
  color(color),
  width_ratio_dump(width_ratio_dump),
  time_dump(time_dump),
  domain(domain),
  resolution(resolution),
  t0(t0),
  curve(curve)
{
	for (int i=0; i<resolution; i++) {
		float a = lerp(domain.x, domain.y, 1.f*i/resolution);
		float b = lerp(domain.x, domain.y, 1.f*(i+1)/resolution);
		auto line = PlanarDiffusedInterval(X, dt, steps, color, time_dump, width_ratio_dump, curve(a), curve(b), t0);
		merge(line);
	}
}

PlanarDiffusedPatterns::PlanarDiffusedPatterns(const VectorFieldR2 &X, float dt, int steps, const vector<vec4> &colors, const vector<vec2> &shifts, const std::function<float(float)> &time_dump, const std::function<float(float)> &width_ratio_dump, const SmoothParametricPlaneCurve &curve_pattern, vec2 domain, int resolution, float t0) {
	for (int i=0; i<colors.size(); i++) {
		auto c = curve_pattern + shifts[i];
		auto line = PlanarDiffusedCurve(X, dt, steps, colors[i], time_dump, width_ratio_dump, c, domain, resolution, t0);
		merge(line);
	}
}

PipeCurveVertexShader::PipeCurveVertexShader(const SmoothParametricCurve &curve, float r, int horRes, int radialRes)
: r(r), id(randomID())
{
	boss = make_unique<BufferManager>(std::set({POSITION, NORMAL, UV, COLOR, INDEX, EXTRA0}));

	auto params = linspace(curve.getT0(), curve.getT1(), horRes);
	vector<vec3> normals = vector<vec3>();
	vector<vec3> binormals = vector<vec3>();

	for (float t: params) {
		vec3 b = curve.binormal(t);
		vec3 n = curve.normal(t);
		normals.push_back(n);
		binormals.push_back(b);
	}

	vector<Vertex> verts = vector<Vertex>();
	vector<ivec3> inds = vector<ivec3>();


	for (int i = 0; i < horRes; i++) {
		float t = params[i];
		vec3 p = curve(t);
		vec3 b = binormals[i];
		vec3 n = normals[i];
		for (float theta: linspace(0.f, TAU, radialRes)) {
			verts.emplace_back(p, vec2(t, theta), n, vec4(b.x, b.y, b.z, r));
			verts.back().addExtraData("extra0", vec4(0));
		}
		if (i < horRes - 1) {
			for (int j = 0; j < radialRes; j++) {
				inds.emplace_back(i*radialRes+j, (i+1)*radialRes+j, i*radialRes+(j+1)%radialRes);
				inds.emplace_back((i+1)*radialRes+j, (i+1)*radialRes+(j+1)%radialRes, i*radialRes+(j+1)%radialRes);
			}
		}
	}
	addNewPolygroup(verts, inds, id);

}

PipeCurveVertexShader::PipeCurveVertexShader(const RealFunction &plot, vec2 dom, float r, int horRes, int radialRes) : id(randomID()) {
	boss = make_unique<BufferManager>(std::set({POSITION, NORMAL, UV, COLOR, INDEX, EXTRA0}));


	auto params = linspace(dom.x, dom.y, horRes);

	vector<Vertex> verts = vector<Vertex>();
	vector<ivec3> inds = vector<ivec3>();


	for (int i = 0; i < horRes; i++) {
		float t = params[i];
		vec3 p = vec3(t, 0, plot(t));
		vec3 b = e2;
		vec3 n = normalise(vec3(plot.df(t), 0, -1));
		for (float theta: linspace(0.f, TAU, radialRes)) {
			verts.emplace_back(p, vec2(t, theta), n, vec4(b.x, b.y, b.z, r));
			verts.back().addExtraData("extra0", vec4(0));
		}
		if (i < horRes - 1) {
			for (int j = 0; j < radialRes; j++) {
				inds.emplace_back(i*radialRes+j, (i+1)*radialRes+j, i*radialRes+(j+1)%radialRes);
				inds.emplace_back((i+1)*radialRes+j, (i+1)*radialRes+(j+1)%radialRes, i*radialRes+(j+1)%radialRes);
			}
		}
	}
	addNewPolygroup(verts, inds, id);
}

PipeCurveVertexShader::PipeCurveVertexShader(const DiscreteRealFunction &plot, float r, int radialRes) : id(randomID()) {
	boss = make_unique<BufferManager>(std::set({POSITION, NORMAL, UV, COLOR, INDEX, EXTRA0}));


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
		for (float theta: linspace(0.f, TAU, radialRes)) {
			verts.emplace_back(p, vec2(t, theta), n, vec4(b.x, b.y, b.z, r));
			verts.back().addExtraData("extra0", vec4(0));
		}
		if (i < horRes - 1) {
			for (int j = 0; j < radialRes; j++) {
				inds.emplace_back(i*radialRes+j, (i+1)*radialRes+j, i*radialRes+(j+1)%radialRes);
				inds.emplace_back((i+1)*radialRes+j, (i+1)*radialRes+(j+1)%radialRes, i*radialRes+(j+1)%radialRes);
			}
		}
	}
	addNewPolygroup(verts, inds, id);
}

void PipeCurveVertexShader::updateCurve(const SmoothParametricCurve &curve) {
	deformPerVertex(id, [this, &curve](BufferedVertex &v) {
		float t = v.getUV().x;
		vec3 p = curve(t);
		vec3 b = curve.binormal(t);
		vec3 n = curve.normal(t);
		if (dot(n, v.getNormal()) < 0)
			n = -n;
		if (dot(b, vec3(v.getColor())) < 0)
			b = -b;
		v.setPosition(p);
		v.setNormal(n);
		v.setColor(b.x, 0);
		v.setColor(b.y, 1);
		v.setColor(b.z, 2);
	});
}

void PipeCurveVertexShader::duplicateCurve(const PolyGroupID &copy_id) {
	copyPolygroup(id, copy_id);
}

void PipeCurveVertexShader::updateCurve(const DiscreteRealFunction &plot) {
	auto df = plot.derivative();
	deformPerVertex(id, [this, &plot, &df](BufferedVertex &v) {
		float t = v.getUV().x;
		vec3 p = vec3(t, 0, plot(t));
		vec3 n = normalize(vec3(df(t), 0, -1));
		v.setPosition(p);
		v.setNormal(n);
	});}


void PipeCurveVertexShader::updateCurve(const DiscreteRealFunction &plot, const DiscreteRealFunction &df_precomputed) {
	deformPerVertex(id, [this, &plot, &df_precomputed](BufferedVertex &v) {
		float t = v.getUV().x;
		vec3 p = vec3(t, 0, plot(t));
		vec3 n = normalize(vec3(df_precomputed(t), 0, -1));
		v.setPosition(p);
		v.setNormal(n);
	});}

void PipeCurveVertexShader::updateRadius(const HOM(float, float) &r) {
	deformPerVertex(id, [this, &r](BufferedVertex &v) {
		float t = v.getUV().x;
		v.setColor(r(t), 3);
	});}



std::function<void(float, float)> deformationOperator(const std::function<void(BufferedVertex &, float, float)> &deformation, WeakSuperMesh &mesh, const PolyGroupID &id) {
    return [&deformation, &mesh, id](float t, float delta) {
        mesh.deformPerVertex(id, [deformation, t, delta](BufferedVertex &v) {
	        deformation(v, t, delta);
        });
    };
}
std::function<void(float)> deformationOperator(const std::function<void(BufferedVertex &, float)> &deformation, WeakSuperMesh &mesh, const PolyGroupID &id) {
    return [&deformation, &mesh, id](float t) {
        mesh.deformPerVertex(id, [deformation, t](BufferedVertex &v) {
	        deformation(v, t);
        });
    };
}

std::function<void(float, float)> moveAlongCurve(const SmoothParametricCurve &curve, WeakSuperMesh &mesh, const PolyGroupID &id) {
    return deformationOperator([curve](BufferedVertex &v, float t, float delta) {
        v.setPosition(v.getPosition() + curve(t) - curve(t - delta));
    }, mesh, id);
}

// ReSharper disable once CppPassValueParameterByConstReference
mat3 WeakSuperMesh::inertiaTensorCMAppBd(PolyGroupID id) const {
	return inertiaTensorAppBd(id, centerOfMass(id));
}

vector<int> WeakSuperMesh::findVertexNeighbours(int i, const PolyGroupID &id) const {
	std::set<int> neighbours = {};
	for (const IndexedTriangle &t: triangles.at(id))
		if (contains(i, t.getVertexIndices())) {
			neighbours.insert(setMinus(t.getVertexIndices(), i)[0]);
			neighbours.insert(setMinus(t.getVertexIndices(), i)[1]);
		}
	return vector(neighbours.begin(), neighbours.end());
}

vector<int> WeakSuperMesh::findVertexParentTriangles(int i, const PolyGroupID &id) const {
	std::set<int> parentTriangles = {};
	for (int j = 0; j < triangles.at(id).size(); j++)
		if (contains(i, triangles.at(id)[j].getVertexIndices()))
			parentTriangles.insert(j);
	return vector(parentTriangles.begin(), parentTriangles.end());
}

void WeakSuperMesh::recalculateNormal(int i, const PolyGroupID &id) {
	auto trs = findVertexParentTriangles(i, id);
	vec3 n   = vec3(0);
	for (int j: trs)
		n += triangles.at(id)[j].faceNormal()*triangles.at(id)[j].area();
	vertices.at(id)[i].setNormal(normalize(n));
}

void WeakSuperMesh::recalculateNormalsNearby(int i, const PolyGroupID &id) {
	for (int j: findVertexNeighbours(i, id))
		recalculateNormal(j, id);
}

void WeakSuperMesh::recalculateNormals(const PolyGroupID &id) {
	for (int i = 0; i < vertices.at(id).size(); i++)
		recalculateNormal(i, id);
}

void WeakSuperMesh::recalculateNormals() {
	for (auto id: getPolyGroupIDs()) recalculateNormals(id);
}

void WeakSuperMesh::orientFaces(const PolyGroupID &id) {
	for (auto &t: triangles.at(id))
		if (dot(t.faceNormal(), t.getVertex(0).getNormal()) < 0)
			t.changeOrientation();
}

void WeakSuperMesh::orientFaces() { for (auto id: getPolyGroupIDs()) orientFaces(id); }

vector<int> WeakSuperMesh::findNeighboursSorted(int i, const PolyGroupID &id) const {
	vector<int> neighbours = findVertexNeighbours(i, id);
	std::map<int, float> angles = {};
	auto t = orthogonalComplementBasis(vertices.at(id)[i].getNormal());
	for (int j = 0; j < neighbours.size(); j++)
		angles[i] = polarAngle(vertices.at(id)[neighbours[j]].getPosition() - vertices.at(id)[i].getPosition(), vertices.at(id)[i].getNormal());
	std::ranges::sort(neighbours, [&angles](int a, int b) { return angles[a] < angles[b]; });
	return neighbours;
}

bool WeakSuperMesh::checkIfHasCompleteNeighbourhood(int i, const PolyGroupID &id) const {
	vector<int> trs = findVertexParentTriangles(i, id);
	vector<int> neighbours = findVertexNeighbours(i, id);
	for (int p: neighbours) {
		int found = 0;
		for (int t: trs)
			if (triangles.at(id)[t].containsEdge(i, p))
				found++;
		if (found < 2)
			return false;
	}
	return true;
}

float WeakSuperMesh::meanCurvature(int i, const PolyGroupID &id) const {
	if (!checkIfHasCompleteNeighbourhood(i, id))
		return 0;
	float sum = 0;
	auto nbhd = findNeighboursSorted(i, id);
	vec3 p = vertices.at(id)[i].getPosition();
	for (int j=0; j<nbhd.size(); j++) {
		vec3 prev = vertices.at(id)[nbhd[(j-1+nbhd.size())%nbhd.size()]].getPosition();
		vec3 next = vertices.at(id)[nbhd[(j+1)%nbhd.size()]].getPosition();
		vec3 current = vertices.at(id)[nbhd[j]].getPosition();
		float angle1 = abs(angle(prev-p, prev-current));
		float angle2 = abs(angle(next-p, next-current));
		sum += 0.5f*(cot(angle1) + cot(angle2))*length(current-p);
	}
//	return sum/
	return sum;
}

vec3 WeakSuperMesh::meanCurvatureVector(int i, const PolyGroupID &id) const {
	return meanCurvature(i, id)*vertices.at(id)[i].getNormal();
}

void WeakSuperMesh::meanCurvatureFlowDeform(float dt, const PolyGroupID &id) {
	for (int i=0; i<vertices.at(id).size(); i++) {
		vertices.at(id)[i].setPosition(vertices.at(id)[i].getPosition() - dt*meanCurvatureVector(i, id));
		recalculateNormalsNearby(i, id);
	}
	recalculateNormals(id);
}

float WeakSuperMesh::GaussCurvature(int i, const PolyGroupID &id) const {
	auto nbhd = findNeighboursSorted(i, id);
	float sum = 0;
	vec3 p = vertices.at(id)[i].getPosition();
	for (int j=0; j<nbhd.size(); j++) {
		vec3 prev = vertices.at(id)[nbhd[(j-1+nbhd.size())%nbhd.size()]].getPosition();
		vec3 next = vertices.at(id)[nbhd[(j+1)%nbhd.size()]].getPosition();
		vec3 current = vertices.at(id)[nbhd[j]].getPosition();
		sum += angle(prev-current, next-current)/length(prev-next)*length(p-current)/2;
	}
	return sum;
}

void WeakSuperMesh::paintMeanCurvature(const PolyGroupID &id) {
	deformPerVertex(id, [this, id](BufferedVertex &v) {
		v.setColor(meanCurvature(v.getIndex(), id), 2);
	});
}

void WeakSuperMesh::paintMeanCurvature() {
	for (auto id: getPolyGroupIDs())
		paintMeanCurvature(id);
}



mat3 WeakSuperMesh::inertiaTensorAppBd(PolyGroupID id, vec3 p) const {
	vec3 cm = p;
	return integrateOverTriangles<mat3>([ cm](const IndexedTriangle &t) {
		vec3 p = t.center();
		return mat3(pow((p-cm).z, 2) + pow((p-cm).y, 2), -(p-cm).x*(p-cm).y, -(p-cm).x*(p-cm).z,
					-(p-cm).x*(p-cm).y, pow((p-cm).x, 2) + pow((p-cm).z, 2), -(p-cm).y*(p-cm).z,
					-(p-cm).x*(p-cm).z, -(p-cm).y*(p-cm).z, pow((p-cm).x, 2) + pow((p-cm).y, 2));
	}, id);
}


SurfacePlotDiscretisedMesh::SurfacePlotDiscretisedMesh(const DiscreteRealFunctionR2 &plot) {
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
				auto p1 = points.back().getPosition();
				auto p2 = points.at(points.size()-plot.samples_x()).getPosition();
				n = normalize(cross(p1-p, p2-p));
				if (dot(n, e3) < 0) n = -n;
				triangles.emplace_back((i-1)*plot.samples_x()+j-1, (i-1)*plot.samples_x()+j, i*plot.samples_x()+j);
				triangles.emplace_back((i-1)*plot.samples_x()+j-1, i*plot.samples_x()+j, i*plot.samples_x()+j-1);
			}
			points.emplace_back(p, vec2(x, t), n, vec4(0));
		}
	}
	addNewPolygroup(points, triangles, randomID());

}

SurfacePolarPlotDiscretisedMesh::SurfacePolarPlotDiscretisedMesh(const DiscreteRealFunctionR2 &plot, float r, float rot_speed) {
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
			float phi = xs[j] + t*rot_speed;
			// float R = r-f[i][j];
			float R = r-log(1+f[i][j]);

			vec3 p = vec3(R*sin(phi), t, R*cos(phi));
			vec3 n = -vec3(sin(phi), 0, cos(phi));
			if (i > 0) {
				auto p1 = points.back().getPosition();
				auto p2 = points.at(points.size()-nx).getPosition();
				if (j == 0)
					p2 = points.at(points.size()-2).getPosition();
				n = normalize(cross(p1-p, p2-p));
				if (dot(n, vec3(sin(phi), 0, cos(phi))) < 0) n = -n;
				triangles.emplace_back((i-1)*nx+(j-1+nx)%nx, (i-1)*nx+j, i*nx+j);
				triangles.emplace_back((i-1)*nx+(j-1+nx)%nx, i*nx+j, i*nx+(j-1+nx)%nx);
			}
			points.emplace_back(p, vec2(phi, t), n, vec4(0));
		}
	}
	addNewPolygroup(points, triangles, randomID());}
