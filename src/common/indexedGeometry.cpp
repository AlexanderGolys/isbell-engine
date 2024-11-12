#include "indexedGeometry.hpp"

using namespace glm;
using std::vector, std::string, std::map, std::shared_ptr, std::unique_ptr, std::pair, std::make_unique, std::make_shared, std::array, std::weak_ptr;

BufferManager::BufferManager(std::set<CommonBufferType> activeBuffers) {
    stds = make_unique<buffStds>();
    mater =  activeBuffers.contains(MATERIAL1) ? make_unique<buff4x4>() : nullptr;

    extra0 = activeBuffers.contains(EXTRA0) ? make_unique<buff4>() : nullptr;
    extra = nullptr;
    if (activeBuffers.contains(EXTRA1) || activeBuffers.contains(EXTRA2) || activeBuffers.contains(EXTRA3) || activeBuffers.contains(EXTRA4))
        extra = make_unique<buff4x4>();
    indices = activeBuffers.contains(INDEX) ? make_unique<ibuff3>() : nullptr;
    this->activeBuffers = activeBuffers;
}
BufferManager::BufferManager(bool materials, std::set<CommonBufferType> extras) : BufferManager({POSITION, NORMAL, UV, COLOR, INDEX}) {
    if (materials) {
        mater = make_unique<buff4x4>();
        this->activeBuffers = {POSITION, NORMAL, UV, COLOR, MATERIAL1, MATERIAL2, MATERIAL3, MATERIAL4, INDEX};
    }
    extra0 = extras.contains(EXTRA0) ? make_unique<buff4>() : nullptr;
    extra = nullptr;
    if (extras.contains(EXTRA1) || extras.contains(EXTRA2) || extras.contains(EXTRA3) || extras.contains(EXTRA4))
        extra = make_unique<buff4x4>();
    for (CommonBufferType type: extras)
        activeBuffers.insert(type);
}

int BufferManager::bufferLength(CommonBufferType type) const {
    switch (type) {
        case POSITION:
        case NORMAL:
        case UV:
        case COLOR:
            return stds->size();
        case MATERIAL1:
        case MATERIAL2:
        case MATERIAL3:
        case MATERIAL4:
            return mater->size();
        case EXTRA0:
            return extra0->size();
        case EXTRA1:
        case EXTRA2:
        case EXTRA3:
        case EXTRA4:
            return extra->size();
        case INDEX:
            return indices->size();
    }
    throw UnknownVariantError("Buffer not recognised among common types. ");
}

int BufferManager::offset(CommonBufferType type) {
    switch (type) {
        case POSITION:
        case MATERIAL1:
        case EXTRA0:
        case EXTRA1:
        case INDEX:
            return 0;
        case NORMAL:
            return 3*sizeof(FLOAT);
        case UV:
            return 6*sizeof(FLOAT);
        case COLOR:
        case MATERIAL3:
        case EXTRA3:
            return 8*sizeof(FLOAT);
        case MATERIAL2:
        case EXTRA2:
            return 4*sizeof(FLOAT);

        case MATERIAL4:
        case EXTRA4:
            return 12*sizeof(FLOAT);
        default:
            throw UnknownVariantError("Buffer not recognised among common types. ");
    }
    return -1;
}

void *BufferManager::firstElementAddress(CommonBufferType type) const {
    switch (type) {
    case POSITION:
    case NORMAL:
    case UV:
    case COLOR:
        return &stds->at(0);
    case MATERIAL1:
    case MATERIAL2:
    case MATERIAL3:
    case MATERIAL4:
        return &mater->at(0);
    case EXTRA0:
        return &extra0->at(0);
    case EXTRA1:
    case EXTRA2:
    case EXTRA3:
    case EXTRA4:
        return &extra->at(0);
    case INDEX:
        return &indices->at(0);
    }
    throw UnknownVariantError("Buffer not recognised among common types. ");
}



int BufferManager::addMaterialBufferData(mat4 mat) {
    return addMaterialBufferData(mat[0], mat[1], mat[2], mat[3]);
}
int BufferManager::addMaterialBufferData(const MaterialPhong &mat) {
    return addMaterialBufferData(mat.compressToMatrix());
}

int BufferManager::addStdAttributesFromVertex(vec3 pos, vec3 norm, vec2 uv, vec4 col) {
    stds->push_back({pos, norm, uv, col});
    return stds->size() - 1;
}
int BufferManager::addMaterialBufferData(vec4 ambient, vec4 diffuse, vec4 specular, vec4 intencity) {
    mater->push_back({ambient, diffuse, specular, intencity});
    return mater->size() - 1;
}
int BufferManager::addFullVertexData(vec3 pos, vec3 norm, vec2 uv, vec4 col, mat4 mat, bool material) {
    if (material)
        addMaterialBufferData(mat);
    return addStdAttributesFromVertex(pos, norm, uv, col);
}
int BufferManager::addFullVertexData(const Vertex &v, bool material) {
    if (material)
        return addFullVertexData(v.getPosition(), v.getNormal(), v.getUV(), v.getColor(), v.getMaterial().compressToMatrix(), material);
    return addFullVertexData(v.getPosition(), v.getNormal(), v.getUV(), v.getColor(), mat4(0), material);
}


BufferedVertex::BufferedVertex(BufferManager &bufferBoss, const Vertex &v, bool material) : bufferBoss(bufferBoss) {
    index = this->bufferBoss.addFullVertexData(v, material);
}

void BufferedVertex::applyFunction(SpaceEndomorphism f) {
    vec3 p = getPosition();
    vec3 n = getNormal();
    setPosition(f(p));
    setNormal(f.df(p)*n);
}


IndexedTriangle::IndexedTriangle(ivec3 indices) {
    metaIndex = -1;
    indicesWithinPolygroup = indices;
}

IndexedTriangle::IndexedTriangle(ivec3 indices, const vector<BufferedVertex> &arrayWithinPoly, BufferManager &bufferBoss) :
    IndexedTriangle(indices) {
    addToBuffer(arrayWithinPoly, bufferBoss);
}


ivec3 IndexedTriangle::bufferIndices(const std::vector<BufferedVertex> &arrayWithinPoly) const {
    return ivec3(arrayWithinPoly[indicesWithinPolygroup.x].getIndex(), arrayWithinPoly[indicesWithinPolygroup.y].getIndex(), arrayWithinPoly[indicesWithinPolygroup.z].getIndex());
}

void IndexedTriangle::addToBuffer(const std::vector<BufferedVertex> &arrayWithinPoly, BufferManager &bufferBoss) {
    metaIndex = bufferBoss.addTriangleVertexIndices(bufferIndices(arrayWithinPoly));
}


glm::mat3 IndexedTriangle::orthonormalFrame(const std::vector<BufferedVertex> &arrayWithinPoly) const {
    vec3 p0 = getVertex(0, arrayWithinPoly).getPosition();
    vec3 p1 = getVertex(1, arrayWithinPoly).getPosition();
    vec3 p2 = getVertex(2, arrayWithinPoly).getPosition();
    return  GramSchmidtProcess(mat3(p0-p2, p1-p2, cross(p1-p2, p0-p2)));
}

glm::vec3 IndexedTriangle::fromPlanar(glm::vec2 v, const std::vector<BufferedVertex> &arrayWithinPoly) const {
    mat3 frame = orthonormalFrame(arrayWithinPoly);
    return frame[0]*v.x + frame[1]*v.y + getVertex(2, arrayWithinPoly).getPosition();
}

glm::vec2 IndexedTriangle::toPlanar(glm::vec3 v, const std::vector<BufferedVertex> &arrayWithinPoly) const {
    mat3 frame = orthonormalFrame(arrayWithinPoly);
    return vec2(dot(frame[0], v-getVertex(2, arrayWithinPoly).getPosition()), dot(frame[1], v-getVertex(2, arrayWithinPoly).getPosition()));
}


glm::vec3 IndexedTriangle::fromBars(glm::vec2 v, const std::vector<BufferedVertex> &arrayWithinPoly) const {
    return getVertex(0, arrayWithinPoly).getPosition() * v.x + getVertex(1, arrayWithinPoly).getPosition() * v.y +
           getVertex(2, arrayWithinPoly).getPosition() * (1 - v.x - v.y);
}


std::array<glm::vec3, 3> IndexedTriangle::borderTriangle(float width, const std::vector<BufferedVertex> &arrayWithinPoly) const {
    vec3 p0 = getVertex(0, arrayWithinPoly).getPosition();
    vec3 p1 = getVertex(1, arrayWithinPoly).getPosition();
    vec3 p2 = getVertex(2, arrayWithinPoly).getPosition();
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

    vec2 b0_01 = toPlanar(p0_01, arrayWithinPoly);
    vec2 b1_01 = toPlanar(p1_01, arrayWithinPoly);
    vec2 b1_12 = toPlanar(p1_12, arrayWithinPoly);
    vec2 b2_12 = toPlanar(p2_12, arrayWithinPoly);
    vec2 b0_20 = toPlanar(p0_20, arrayWithinPoly);
    vec2 b2_20 = toPlanar(p2_20, arrayWithinPoly);

    vec2 b0 = intersectLines(b0_01, b1_01, b0_20, b2_20);
    vec2 b1 = intersectLines(b1_01, b0_01, b1_12, b2_12);
    vec2 b2 = intersectLines(b2_20, b0_20, b2_12, b1_12);

    return {fromPlanar(b0, arrayWithinPoly), fromPlanar(b1, arrayWithinPoly), fromPlanar(b2, arrayWithinPoly)};
}


WeakSuperMesh::WeakSuperMesh() {boss = make_unique<BufferManager>();}

WeakSuperMesh::WeakSuperMesh(const vector<Vertex> &hardVertices, const vector<ivec3> &faceIndices, PolyGroupID id, MaterialPhong *material) : WeakSuperMesh() {
    addNewPolygroup(hardVertices, faceIndices, id, material);
}

WeakSuperMesh::WeakSuperMesh(SmoothParametricSurface surf, int tRes, int uRes, std::variant<int, std::string> id, const MaterialPhong &material) : WeakSuperMesh() {
    addUniformSurface(surf, tRes, uRes, id, material);
}


void WeakSuperMesh::addNewPolygroup(const vector<Vertex> &hardVertices, const vector<ivec3> &faceIndices, PolyGroupID id, MaterialPhong *material) {
    if (vertices.contains(id) || triangles.contains(id))
        throw IllegalVariantError("Polygroup ID already exists in mesh. ");

    vertices[id] = vector<BufferedVertex>();
    triangles[id] = vector<IndexedTriangle>();
    vertices[id].reserve(hardVertices.size());
    triangles[id].reserve(faceIndices.size());

    if (material)
        materials[id] = *material;

    for (const Vertex &v: hardVertices)
        vertices[id].emplace_back(*boss, v, material != nullptr);

    for (const ivec3 &ind: faceIndices)
        triangles[id].emplace_back(ind, vertices[id], *boss);
}



void WeakSuperMesh::addUniformSurface(SmoothParametricSurface surf, int tRes, int uRes, std::variant<int, std::string> id,
                                      const MaterialPhong &material) {
    vector<Vertex> hardVertices = {};
    vector<ivec3> faceIndices = {};
    hardVertices.reserve(tRes * uRes + uRes + tRes + 1);
    faceIndices.reserve(2 * (tRes) * (uRes));

    for (int i = 0; i < tRes; i++) {
        for (int j = 0; j < uRes; j++) {
            float t_ = 1.f * i / (1.f * tRes - 1);
            float u_ = 1.f * j / (1.f * uRes - 1);
            float t = lerp(surf.tMin(), surf.tMax(), t_);
            float u = lerp(surf.uMin(), surf.uMax(), u_);
            hardVertices.emplace_back(surf(t, u), vec2(t_, u_), surf.normal(t, u), BLACK, material);
        }
    }

    for (int i = 0; i < tRes - 1; i++) {
        for (int j = 0; j < uRes - 1; j++) {
            int i0 = j + i * uRes;
            int i1 = j + (i + 1) * uRes;
            int i2 = j + 1 + i * uRes;
            int i3 = j + 1 + (i + 1) * uRes;
            faceIndices.emplace_back(i0, i1, i3);
            faceIndices.emplace_back(i0, i2, i3);
        }
    }
    // if (surf.isPeriodicT()) {
    //     for (int j = 1; j < uRes; j++) {
    //         int i3 = j;
    //         int i2 = j-1;
    //         int i1 = j + (tRes-1)*uRes;
    //         int i0 = j  + (tRes-1)*uRes - 1;
    //         faceIndices.emplace_back(i0, i1, i3);
    //         faceIndices.emplace_back(i0, i2, i3);
    //     }
    // }

    // if (surf.isPeriodicU()) {
    //     for (int i = 0; i < tRes; i++) {
    //         int i3 = (i+1) * uRes;
    //         int i2 = i * uRes;
    //         int i1 = (i+1) * (uRes-1);
    //         int i0 =  i* (uRes-1);
    //         faceIndices.emplace_back(i0, i1, i3);
    //         faceIndices.emplace_back(i0, i2, i3);
    //     }


    addNewPolygroup(hardVertices, faceIndices, id);
}

std::vector<std::variant<int, std::string>> WeakSuperMesh::getPolyGroupIDs() const {
    vector<PolyGroupID> ids = {};
    for (const pair<PolyGroupID, vector<BufferedVertex>> &p: vertices)
        ids.push_back(p.first);
    return ids;
}


void WeakSuperMesh::moveAlongVectorField(PolyGroupID id, VectorFieldR3 X, float delta) {
    for (BufferedVertex &v: vertices.at(id))
        v.setPosition(X.moveAlong(v.getPosition(), delta));
}

void WeakSuperMesh::deformWithAmbientMap(PolyGroupID id, SpaceEndomorphism f) {
    for (BufferedVertex &v: vertices.at(id))
        v.applyFunction(f);
}

vector<Vertex> WeakSuperMesh::getVertices(std::variant<int, std::string> id) const {
    vector<Vertex> verts = {};
    verts.reserve(vertices.at(id).size());
    for (const BufferedVertex &v: vertices.at(id))
        verts.push_back(v.getVertex());
    return verts;
}


vector<ivec3> WeakSuperMesh::getIndices(std::variant<int, std::string> id) const {
    vector<ivec3> inds = {};
    inds.reserve(triangles.at(id).size());
    for (const IndexedTriangle &t: triangles.at(id))
        inds.push_back(t.getLocalIndices());
    return inds;
}


WeakSuperMesh WeakSuperMesh::subdivideBarycentric(std::variant<int, std::string> id) const {
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
        newInds.push_back(ivec3(tr.z, tr.x, centerInd));
    }
    return WeakSuperMesh(verts, newInds, id);
}
WeakSuperMesh WeakSuperMesh::subdivideEdgecentric(std::variant<int, std::string> id) const {
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
        array<vec3, 3> border = tr.borderTriangle(width, verts);
        vec3 n = tr.faceNormal(verts)*sign(dot(tr.getVertex(0, verts).getNormal(), tr.faceNormal(verts)));


        for (vec3 p : border) {
            new_verts.emplace_back(p+n*heightSide, vec2(0, 0), n, BLACK, MaterialPhong());
        }
        for (vec3 p : border) {
            new_verts.emplace_back(p, vec2(0, 0), n, BLACK, MaterialPhong());
        }
        new_inds.push_back(ivec3(tr.getLocalIndices().x, tr.getLocalIndices().y, tr.getLocalIndices().x + verts.size()));
        new_inds.push_back(ivec3(tr.getLocalIndices().y, tr.getLocalIndices().y + verts.size(), tr.getLocalIndices().x + verts.size()));
        new_inds.push_back(ivec3(tr.getLocalIndices().x+verts.size(), tr.getLocalIndices().y+verts.size(), new_verts.size()-6));
        new_inds.push_back(ivec3(tr.getLocalIndices().y+verts.size(), new_verts.size()-5, new_verts.size()-6));
        // new_inds.push_back(ivec3(new_verts.size()-6,new_verts.size()-5, new_verts.size()-3));
        // new_inds.push_back(ivec3(new_verts.size()-5,new_verts.size()-2, new_verts.size()-3));

        new_inds.push_back(ivec3(tr.getLocalIndices().x, tr.getLocalIndices().z, tr.getLocalIndices().x + verts.size()));
        new_inds.push_back(ivec3(tr.getLocalIndices().z, tr.getLocalIndices().z + verts.size(), tr.getLocalIndices().x + verts.size()));
        // new_inds.push_back(ivec3(tr.getLocalIndices().x+verts.size(), tr.getLocalIndices().z+verts.size(), new_verts.size()-6));
        // new_inds.push_back(ivec3(tr.getLocalIndices().z+verts.size(), new_verts.size()-4, new_verts.size()-6));
        // new_inds.push_back(ivec3(new_verts.size()-6,new_verts.size()-4, new_verts.size()-3));
        // new_inds.push_back(ivec3(new_verts.size()-4,new_verts.size()-1, new_verts.size()-3));

        new_inds.push_back(ivec3(tr.getLocalIndices().y, tr.getLocalIndices().z, tr.getLocalIndices().y + verts.size()));
        new_inds.push_back(ivec3(tr.getLocalIndices().z, tr.getLocalIndices().z + verts.size(), tr.getLocalIndices().y + verts.size()));
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
        buff2* b = static_cast<buff2*>(firstElementAddress(type));
        b->push_back(*v);
        return;
    }
   if (bufferElementLength(type) == 3) {
        vec3* v = static_cast<vec3*>(valueAddress);
        buff3* b = static_cast<buff3*>(firstElementAddress(type));
        b->push_back(*v);
        return;
    }
   if (bufferElementLength(type) == 4) {
        vec4* v = static_cast<vec4*>(valueAddress);
        buff4* b = static_cast<buff4*>(firstElementAddress(type));
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
        buff2 *b = static_cast<buff2 *>(firstElementAddress(type));
        b->emplace_back(0, 0);
        return;
    }
    if (bufferElementLength(type) == 3) {
        buff3 *b = static_cast<buff3 *>(firstElementAddress(type));
        b->emplace_back(0, 0, 0);
        return;
    }
    if (bufferElementLength(type) == 4) {
        buff4 *b = static_cast<buff4 *>(firstElementAddress(type));
        b->emplace_back(0, 0, 0, 0);
        return;
    }
    throw IllegalVariantError("Buffer element has illegal length.");
}

void BufferManager::initialiseExtraBufferSlot(int slot) {
    switch (slot) {
        case 0:
            if (extra0 == nullptr)
                extra0 = make_unique<buff4>();
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

        }
    }

glm::vec4 BufferManager::getExtra(int index, int slot) const {
    switch (slot) {
    case 0:
        return (*extra0)[index];
    case 1:
        return (*extra)[index].x;
    case 2:
        return  (*extra)[index].y;
    case 3:
        return (*extra)[index].z;
    case 4:
        return (*extra)[index].w;
    }
    throw UnknownVariantError("Extra slot not recognised. ");
}
void BufferManager::setExtra(int index, glm::vec4 value, int slot) {
    switch (slot) {
        case 0:
            (*extra0)[index] = value;
            return;
        case 1:
            (*extra)[index].x = value;
            return;
        case 2:
            (*extra)[index].y = value;
            return;
        case 3:
            (*extra)[index].z = value;
            return;
        case 4:
            (*extra)[index].w = value;
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
            (*extra)[index].x[component] = value;
            return;
        case 2:
            (*extra)[index].y[component] = value;
            return;
        case 3:
            (*extra)[index].z[component] = value;
            return;
        case 4:
            (*extra)[index].w[component] = value;
            return;
        default:
            throw UnknownVariantError("Extra slot not recognised. ");

    }
}

void BufferManager::reserveSpace(int targetSize) {
    stds->reserve(targetSize);

    if (isActive(MATERIAL1) || isActive(MATERIAL2) || isActive(MATERIAL3) || isActive(MATERIAL4))
        mater->reserve(targetSize);
    if (isActive(EXTRA0))
        extra0->reserve(targetSize);
    if (isActive(EXTRA1) || isActive(EXTRA2) || isActive(EXTRA3) || isActive(EXTRA4))
        extra->reserve(targetSize);

}





