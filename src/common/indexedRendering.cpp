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


int BufferManager::addTriangleVertexIndices(glm::ivec3 ind, int shift) {
    indices->push_back(ind+ivec3(shift));
    return bufferLength(INDEX) - 1;
}

int BufferManager::addStdAttributesFromVertex(vec3 pos, vec3 norm, vec2 uv, vec4 col) {
    stds->positions.push_back(pos);
    stds->normals.push_back(norm);
    stds->uvs.push_back(uv);
    stds->colors.push_back(col);
    return stds->positions.size() - 1;
}

int BufferManager::addFullVertexData(vec3 pos, vec3 norm, vec2 uv, vec4 col) {
    return addStdAttributesFromVertex(pos, norm, uv, col);
}
int BufferManager::addFullVertexData(const Vertex &v) {
    return addStdAttributesFromVertex(v.getPosition(), v.getNormal(), v.getUV(), v.getColor());
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


glm::mat3 IndexedTriangle::orthonormalFrame() const {
    vec3 p0 = getVertex(0).getPosition();
    vec3 p1 = getVertex(1).getPosition();
    vec3 p2 = getVertex(2).getPosition();
    return  GramSchmidtProcess(mat3(p0-p2, p1-p2, cross(p1-p2, p0-p2)));
}

glm::vec3 IndexedTriangle::fromPlanar(glm::vec2 v) const {
    mat3 frame = orthonormalFrame();
    return frame[0]*v.x + frame[1]*v.y + getVertex(2).getPosition();
}

glm::vec2 IndexedTriangle::toPlanar(glm::vec3 v) const {
    mat3 frame = orthonormalFrame();
    return vec2(dot(frame[0], v-getVertex(2).getPosition()), dot(frame[1], v-getVertex(2).getPosition()));
}


glm::vec3 IndexedTriangle::fromBars(glm::vec2 v) const {
    return getVertex(0).getPosition() * v.x + getVertex(1).getPosition() * v.y +
           getVertex(2).getPosition() * (1 - v.x - v.y);
}


std::array<glm::vec3, 3> IndexedTriangle::borderTriangle(float width) const {
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

WeakSuperMesh::WeakSuperMesh(const char *filename, const std::variant<int, std::string> &id) : WeakSuperMesh() {
    addNewPolygroup(filename, id);
}

void WeakSuperMesh::addNewPolygroup(const char *filename, const std::variant<int, std::string> &id) {
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
            glm::vec3 vert;
            s >> vert.x; s >> vert.y; s >> vert.z;
            std::cout << vert.x << " " << vert.y << " " << vert.z << std::endl;
            vertices[id].emplace_back(*boss, Vertex(vert, vec2(2137, 0), vec3(0, 0, 0), BLACK));
            pos.push_back(vert);
        }
        if (line.substr(0,3) == "vt ")
        {
            std::istringstream s(line.substr(3));
            glm::vec2 uv;
            s >> uv.x; s >> uv.y;
            uvs.push_back(uv);
        }

        if (line.substr(0,3) == "vn ")
        {
            std::istringstream s(line.substr(3));
            glm::vec3 norm;
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
            for (int i = 0; i < 3; i++) {
                if (vertices[id][vertexIndex[i]].getNormal() == normsVec[i] && vertices[id][vertexIndex[i]].getUV() == uvsVec[i])
                    face[i] = vertexIndex[i];
                else {
                    vertices[id].emplace_back(*boss, Vertex(vertices[id][vertexIndex[i]].getPosition(), uvsVec[i], normsVec[i], BLACK));
                    face[i] = vertices[id].size() - 1;
                }
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

WeakSuperMesh::WeakSuperMesh(const SmoothParametricSurface &surf, int tRes, int uRes, const PolyGroupID &id) : WeakSuperMesh() {
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
            float u_ = 1.f * j / (1.f * uRes - 1);
            float t = lerp(surf.tMin(), surf.tMax(), t_);
            float u = lerp(surf.uMin(), surf.uMax(), u_);
            hardVertices.emplace_back(surf(t, u), vec2(t_, u_), surf.normal(t, u), vec4(t, u, 0, 1));
        }

    for (int i = 1; i < tRes ; i++)
        for (int j = 1; j < uRes; j++) {
            int i0 =   i * uRes + j;
            int i1 =	(i-1) * uRes + j;
            int i2 =	(i-1) * uRes + j-1;
            int i3 =	i * uRes + j-1;
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


std::vector<PolyGroupID> WeakSuperMesh::getPolyGroupIDs() const {
    vector<PolyGroupID> ids = {};
    for (const pair<PolyGroupID, vector<BufferedVertex>> &p: vertices)
        ids.push_back(p.first);
    return ids;
}

void WeakSuperMesh::deformPerVertex(const std::variant<int, std::string> &id, const std::function<void(int, BufferedVertex &)> &deformation) {
	for (int i=0; i<vertices.at(id).size(); i++)
		deformation(i, vertices.at(id)[i]);
}

vec2 WeakSuperMesh::getSurfaceParameters(const BufferedVertex &v) const { return vec2(v.getColor().x, v.getColor().y); }


void WeakSuperMesh::encodeSurfacePoint(BufferedVertex &v, const SmoothParametricSurface &surf, vec2 tu) {
	v.setPosition(surf(tu));
	v.setNormal(surf.normal(tu));
	v.setUV(vec2((tu.x-surf.tMin())/(surf.tMax()-surf.tMin()), (tu.y-surf.uMin())/(surf.uMax()-surf.uMin())));
	v.setColor(tu.x, 0);
	v.setColor(tu.y, 1);
}

void WeakSuperMesh::adjustToNewSurface(const SmoothParametricSurface &surf, const std::variant<int, std::string> &id) {
	deformPerVertex(id, [&](BufferedVertex &v) {
	encodeSurfacePoint(v, surf, getSurfaceParameters(v));
}); }

void WeakSuperMesh::adjustToNewSurface(const SmoothParametricSurface &surf) { for (auto id: getPolyGroupIDs()) adjustToNewSurface(surf, id); }

void WeakSuperMesh::moveAlongVectorField(const std::variant<int, std::string> &id, VectorFieldR3 X, float delta) {
    for (BufferedVertex &v: vertices.at(id))
        v.setPosition(X.moveAlong(v.getPosition(), delta));
}

void WeakSuperMesh::deformWithAmbientMap(const std::variant<int, std::string> &id, SpaceEndomorphism f) {
    for (BufferedVertex &v: vertices.at(id))
        v.applyFunction(f);
}

vector<Vertex> WeakSuperMesh::getVertices(const std::variant<int, std::string> &id) const {
    vector<Vertex> verts = {};
    verts.reserve(vertices.at(id).size());
    for (const BufferedVertex &v: vertices.at(id))
        verts.push_back(v.getVertex());
    return verts;
}


vector<ivec3> WeakSuperMesh::getIndices(const std::variant<int, std::string> &id) const {
    vector<ivec3> inds = {};
    inds.reserve(triangles.at(id).size());
    for (const IndexedTriangle &t: triangles.at(id))
        inds.push_back(t.getVertexIndices());
    return inds;
}




void WeakSuperMesh::pointNormalsInDirection(vec3 dir, const std::variant<int, std::string> &id) { deformPerVertex(id, [dir](BufferedVertex &v) {
	vec3 n = v.getNormal();
	v.setNormal(n*1.f*sgn(dot(n, dir)));
}); }

WeakSuperMesh WeakSuperMesh::subdivideBarycentric(const std::variant<int, std::string> &id) const {
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


WeakSuperMesh WeakSuperMesh::subdivideEdgecentric(const std::variant<int, std::string> &id) const {
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

glm::vec4 BufferManager::getExtra(int index, int slot) const {
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


void BufferManager::setExtra(int index, glm::vec4 value, int slot) {
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
void BufferManager::setExtra(int index, glm::vec3 value, int slot) {
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

vec3 WeakSuperMesh::centerOfMass(std::variant<int, std::string> id) const {
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


std::function<void(float, float)> deformationOperator(const std::function<void(BufferedVertex &, float, float)> &deformation, WeakSuperMesh &mesh, const std::variant<int, std::string> &id) {
    return [&deformation, &mesh, id](float t, float delta) {
        mesh.deformPerVertex(id, [deformation, t, delta](BufferedVertex &v) { deformation(v, t, delta); });
    };
}
std::function<void(float)> deformationOperator(const std::function<void(BufferedVertex &, float)> &deformation, WeakSuperMesh &mesh, const std::variant<int, std::string> &id) {
    return [&deformation, &mesh, id](float t) {
        mesh.deformPerVertex(id, [deformation, t](BufferedVertex &v) { deformation(v, t); });
    };
}

std::function<void(float, float)> moveAlongCurve(const SmoothParametricCurve &curve, WeakSuperMesh &mesh, const std::variant<int, std::string> &id) {
    return deformationOperator([curve](BufferedVertex &v, float t, float delta) {
        v.setPosition(v.getPosition() + curve(t) - curve(t - delta));
    }, mesh, id);
}
