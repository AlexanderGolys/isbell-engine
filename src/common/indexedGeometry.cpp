#include "indexedGeometry.hpp"

using namespace glm;
using std::vector, std::string, std::map, std::shared_ptr, std::unique_ptr, std::pair, std::make_unique, std::make_shared;

BufferManager::BufferManager(std::set<CommonBufferType> activeBuffers) {
    positions = activeBuffers.contains(POSITION) ? make_unique<buff3>() : nullptr;
    normals = activeBuffers.contains(NORMAL) ? make_unique<buff3>() : nullptr;
    uvs = activeBuffers.contains(UV) ? make_unique<buff2>() : nullptr;
    colors = activeBuffers.contains(COLOR) ? make_unique<buff4>() : nullptr;
    mat1 = activeBuffers.contains(MATERIAL1) ? make_unique<buff4>() : nullptr;
    mat2 = activeBuffers.contains(MATERIAL2) ? make_unique<buff4>() : nullptr;
    mat3 = activeBuffers.contains(MATERIAL3) ? make_unique<buff4>() : nullptr;
    mat4 = activeBuffers.contains(MATERIAL4) ? make_unique<buff4>() : nullptr;
    extra1 = activeBuffers.contains(EXTRA1) ? make_unique<buff4>() : nullptr;
    extra2 = activeBuffers.contains(EXTRA2) ? make_unique<buff4>() : nullptr;
    extra3 = activeBuffers.contains(EXTRA3) ? make_unique<buff4>() : nullptr;
    extra4 = activeBuffers.contains(EXTRA4) ? make_unique<buff4>() : nullptr;
    indices = activeBuffers.contains(INDEX) ? make_unique<ibuff3>() : nullptr;
    this->activeBuffers = activeBuffers;
}

int BufferManager::bufferLength(CommonBufferType type) const {
    switch (type) {
    case POSITION:
        return positions->size();
    case NORMAL:
        return normals->size();
    case UV:
        return uvs->size();
    case COLOR:
        return colors->size();
    case MATERIAL1:
        return mat1->size();
    case MATERIAL2:
        return mat2->size();
    case MATERIAL3:
        return mat3->size();
    case MATERIAL4:
        return mat4->size();
    case EXTRA1:
        return extra1->size();
    case EXTRA2:
        return extra2->size();
    case EXTRA3:
        return extra3->size();
    case EXTRA4:
        return extra4->size();
    case INDEX:
        return indices->size();
    }
    throw UnknownVariantError("Buffer not recognised among common types. ");
}

void *BufferManager::firstElementAddress(CommonBufferType type) const {
    switch (type) {
    case POSITION:
        return &positions->at(0);
    case NORMAL:
        return &normals->at(0);
    case UV:
        return &uvs->at(0);
    case COLOR:
        return &colors->at(0);
    case MATERIAL1:
        return &mat1->at(0);
    case MATERIAL2:
        return &mat2->at(0);
    case MATERIAL3:
        return &mat3->at(0);
    case MATERIAL4:
        return &mat4->at(0);
    case EXTRA1:
        return &extra1->at(0);
    case EXTRA2:
        return &extra2->at(0);
    case EXTRA3:
        return &extra3->at(0);
    case EXTRA4:
        return &extra4->at(0);
    case INDEX:
        return &indices->at(0);
    }
    throw UnknownVariantError("Buffer not recognised among common types. ");
}



int BufferManager::addMaterialBufferData(glm::mat4 mat) {
    return addMaterialBufferData(mat[0], mat[1], mat[2], mat[3]);
}
int BufferManager::addMaterialBufferData(const MaterialPhong &mat) {
    return addMaterialBufferData(mat.compressToMatrix());
}

int BufferManager::addStdAttributesFromVertex(glm::vec3 pos, glm::vec3 norm, glm::vec2 uv, glm::vec4 col) {
    positions->push_back(pos);
    normals->push_back(norm);
    uvs->push_back(uv);
    colors->push_back(col);
    return positions->size() - 1;
}
int BufferManager::addMaterialBufferData(glm::vec4 ambient, glm::vec4 diffuse, glm::vec4 specular, glm::vec4 intencity) {
    mat1->push_back(ambient);
    mat2->push_back(diffuse);
    mat3->push_back(specular);
    mat4->push_back(intencity);
    return mat1->size() - 1;
}
int BufferManager::addFullVertexData(glm::vec3 pos, glm::vec3 norm, glm::vec2 uv, glm::vec4 col, glm::mat4 mat) {
    addMaterialBufferData(mat);
    return addStdAttributesFromVertex(pos, norm, uv, col);
}
int BufferManager::addFullVertexData(const Vertex &v) {
    return addFullVertexData(v.getPosition(), v.getNormal(), v.getUV(), v.getColor(), v.getMaterial().compressToMatrix());
}


BufferedVertex::BufferedVertex(BufferManager &bufferBoss, const Vertex &v) : bufferBoss(bufferBoss) {
    index = this->bufferBoss.addFullVertexData(v);
}

IndexedTriangle::IndexedTriangle(ivec3 indices) {
    metaIndex = -1;
    indicesWithinPolygroup = indices;
}

IndexedTriangle::IndexedTriangle(ivec3 indices, const vector<BufferedVertex> &arrayWithinPoly, BufferManager &bufferBoss)
: IndexedTriangle(indices) {
    addToBuffer(arrayWithinPoly, bufferBoss);
}



ivec3 IndexedTriangle::bufferIndices(const std::vector<BufferedVertex> &arrayWithinPoly) const {
    return ivec3(arrayWithinPoly[indicesWithinPolygroup.x].getIndex(), arrayWithinPoly[indicesWithinPolygroup.y].getIndex(), arrayWithinPoly[indicesWithinPolygroup.z].getIndex());
}

void IndexedTriangle::addToBuffer(const std::vector<BufferedVertex> &arrayWithinPoly, BufferManager &bufferBoss) {
    metaIndex = bufferBoss.addTriangleVertexIndices(bufferIndices(arrayWithinPoly));
}

WeakSuperMesh::WeakSuperMesh() {boss = make_unique<BufferManager>();}

WeakSuperMesh::WeakSuperMesh(const vector<Vertex> &hardVertices, const vector<ivec3> &faceIndices, PolyGroupID id) : WeakSuperMesh() {
    addNewPolygroup(hardVertices, faceIndices, id);
}

void WeakSuperMesh::addNewPolygroup(const vector<Vertex> &hardVertices, const vector<ivec3> &faceIndices, PolyGroupID id) {
    vertices[id] = vector<BufferedVertex>();
    triangles[id] = vector<IndexedTriangle>();
    vertices[id].reserve(hardVertices.size());
    triangles[id].reserve(faceIndices.size());

    for (const Vertex &v: hardVertices)
        vertices[id].emplace_back(*boss, v);

    for (const ivec3 &ind: faceIndices)
        triangles[id].emplace_back(ind, vertices[id], *boss);
}


void WeakSuperMesh::addUniformSurface(SmoothParametricSurface surf, int tRes, int uRes, std::variant<int, std::string> id, const MaterialPhong &material) {
    vector<Vertex> hardVertices = {};
    vector<ivec3> faceIndices = {};
    hardVertices.reserve(tRes * uRes + uRes + tRes + 1);
    faceIndices.reserve(2 * (tRes) * (uRes));
    for (int i = 0; i < tRes; i++) {
        for (int j = 0; j < uRes; j++) {
            float t_ = 1.f * i / (tRes-1);
            float u_ = 1.f * j / (uRes-1);
            float t = lerp(surf.tMin(), surf.tMax(), t_);
            float u = lerp(surf.uMin(), surf.uMax(), u_);
            hardVertices.emplace_back(surf(t, u), vec2(t_, u_), surf.normal(t, u), BLACK,  material);
        }
    }

    for (int i = 0; i < tRes-1; i++) {
        for (int j = 0; j < uRes-1; j++) {
            int i0 = j + i*uRes;
            int i1 = j + (i+1)*uRes;
            int i2 = j+1 + i*uRes;
            int i3 = j+1 + (i+1)*uRes;
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
void BufferManager::reserveSpace(int targetSize) {
    positions->reserve(targetSize);
    normals->reserve(targetSize);
    uvs->reserve(targetSize);
    colors->reserve(targetSize);
    extra1->reserve(targetSize);
    extra2->reserve(targetSize);
    extra3->reserve(targetSize);
    extra4->reserve(targetSize);

    if (isActive(MATERIAL1))
        mat1->reserve(targetSize);
    if (isActive(MATERIAL2))
        mat2->reserve(targetSize);
    if (isActive(MATERIAL3))
        mat3->reserve(targetSize);
    if (isActive(MATERIAL4))
        mat4->reserve(targetSize);
    if (isActive(EXTRA1))
        extra1->reserve(targetSize);
    if (isActive(EXTRA2))
        extra2->reserve(targetSize);
    if (isActive(EXTRA3))
        extra3->reserve(targetSize);
    if (isActive(EXTRA4))
        extra4->reserve(targetSize);
}

void BufferManager::initialiseExtraBufferSlot(int slot) {
    CommonBufferType typ = static_cast<CommonBufferType>(EXTRA1 + slot);
    if (isActive(typ))
        return;
    activeBuffers.insert(typ);
    switch (typ) {
        case EXTRA1: {
            extra1 = make_unique<buff4>();
            extra1->reserve(positions->size());
            for (int i = 0; i < positions->size(); i++)
                extra1->emplace_back(0, 0, 0, 0);
            return;
        }

        case EXTRA2: {
            extra2 = make_unique<buff4>();
            extra2->reserve(positions->size());
            for (int i = 0; i < positions->size(); i++)
                extra2->emplace_back(0, 0, 0, 0);
            return;
        }

        case EXTRA3: {
            extra3 = make_unique<buff4>();
            extra3->reserve(positions->size());
            for (int i = 0; i < positions->size(); i++)
                extra3->emplace_back(0, 0, 0, 0);
            return;
        }

        case EXTRA4: {
            extra4 = make_unique<buff4>();
            extra4->reserve(positions->size());
            for (int i = 0; i < positions->size(); i++)
                extra4->emplace_back(0, 0, 0, 0);
            return;
        }

        default:
            throw UnknownVariantError("Buffer not recognised among 4 extra slots. ");

    }

}


BufferManager::BufferManager(BufferManager &&other) noexcept :
    uvs(std::move(other.uvs)), positions(std::move(other.positions)), normals(std::move(other.normals)), colors(std::move(other.colors)),
    mat1(std::move(other.mat1)), mat2(std::move(other.mat2)), mat3(std::move(other.mat3)), mat4(std::move(other.mat4)),
    extra1(std::move(other.extra1)), extra2(std::move(other.extra2)), extra3(std::move(other.extra3)), extra4(std::move(other.extra4)),
    indices(std::move(other.indices)), activeBuffers(std::move(other.activeBuffers)) {}


BufferManager &BufferManager::operator=(BufferManager &&other) noexcept {
    if (this == &other)
        return *this;
    uvs = std::move(other.uvs);
    positions = std::move(other.positions);
    normals = std::move(other.normals);
    colors = std::move(other.colors);
    mat1 = std::move(other.mat1);
    mat2 = std::move(other.mat2);
    mat3 = std::move(other.mat3);
    mat4 = std::move(other.mat4);
    extra1 = std::move(other.extra1);
    extra2 = std::move(other.extra2);
    extra3 = std::move(other.extra3);
    extra4 = std::move(other.extra4);
    indices = std::move(other.indices);
    activeBuffers = std::move(other.activeBuffers);
    return *this;
}

