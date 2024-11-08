#ifndef INDEXED_HPP
#define INDEXED_HPP

#include "geometry.hpp"

#include <set>

#define buff2 std::vector<glm::vec2>
#define buff3 std::vector<glm::vec3>
#define buff4 std::vector<glm::vec4>
#define ibuff3 std::vector<glm::ivec3>

enum CommonBufferType {
  POSITION, NORMAL, UV, COLOR, MATERIAL1, MATERIAL2, MATERIAL3, MATERIAL4, INDEX, EXTRA1, EXTRA2, EXTRA3, EXTRA4
};
const std::map<CommonBufferType, int> bufferTypeLength = {{POSITION, 3},  {NORMAL, 3},    {UV, 2},        {COLOR, 4},
                                                          {MATERIAL1, 4}, {MATERIAL2, 4}, {MATERIAL3, 4}, {MATERIAL4, 4},
                                                          {INDEX, 3},     {EXTRA1, 4},    {EXTRA2, 4}, {EXTRA3, 4}, {EXTRA4, 4}};

inline int bufferElementLength(CommonBufferType type) { return bufferTypeLength.at(type); }
inline size_t bufferElementSize(CommonBufferType type) { return type == INDEX ? bufferElementLength(type) * sizeof(FLOAT) : bufferElementLength(type) * sizeof(unsigned int); }
class IndexedTriangle;

class BufferManager {
    std::unique_ptr<buff2> uvs;
    std::unique_ptr<buff3> positions, normals;
    std::unique_ptr<buff4> colors, mat1, mat2, mat3, mat4, extra1, extra2, extra3, extra4;
    std::unique_ptr<ibuff3> indices;
    std::set<CommonBufferType> activeBuffers;
    void insertValueToSingleBuffer(CommonBufferType type, void *value);
    void insertDefaultValueToSingleBuffer(CommonBufferType type);


public:
    BufferManager(BufferManager &&other) noexcept;
    BufferManager &operator=(BufferManager &&other) noexcept;
    BufferManager(std::set<CommonBufferType> activeBuffers = {POSITION, NORMAL, UV, COLOR, MATERIAL1, MATERIAL2, MATERIAL3, MATERIAL4, INDEX});
    int bufferLength(CommonBufferType type) const;
    size_t bufferSize(CommonBufferType type) const { return bufferLength(type) * bufferElementSize(type); }
    void *firstElementAddress(CommonBufferType type) const;
    bool isActive(CommonBufferType type) const { return activeBuffers.contains(type); }

    int addTriangleVertexIndices(glm::ivec3 ind) {indices->push_back(ind); return bufferLength(INDEX) - 1; }
    int addStdAttributesFromVertex(glm::vec3 pos, glm::vec3 norm, glm::vec2 uv, glm::vec4 col);
    int addMaterialBufferData(const MaterialPhong &mat);
    int addMaterialBufferData(glm::mat4 mat);
    int addMaterialBufferData(glm::vec4 ambient, glm::vec4 diffuse, glm::vec4 specular, glm::vec4 intencity);
    int addFullVertexData(glm::vec3 pos, glm::vec3 norm, glm::vec2 uv, glm::vec4 col, glm::mat4 mat);
    int addFullVertexData(const Vertex &v);

    void reserveSpace(int targetSize);
    void reserveSpaceForIndex(int targetSize) { indices->reserve(targetSize); }
    void reserveAdditionalSpace(int extraStorage) { reserveSpace(bufferLength(POSITION) + extraStorage); }
    void reserveAdditionalSpaceForIndex(int extraStorage) { reserveSpaceForIndex(bufferLength(INDEX) + extraStorage); }
    void initialiseExtraBufferSlot(int slot);

};



class BufferedVertex {
  BufferManager &bufferBoss;
  int index;
public:
  BufferedVertex(BufferManager &bufferBoss, int index) : bufferBoss(bufferBoss), index(index) {}
  BufferedVertex(const BufferedVertex &other) : bufferBoss(other.bufferBoss), index(other.index) {}
  BufferedVertex(BufferedVertex &&other) noexcept : bufferBoss(other.bufferBoss), index(other.index) {}
  BufferedVertex(BufferManager &bufferBoss, const Vertex &v);

  int getIndex() const { return index; }
};

class IndexedTriangle {
  int metaIndex;
  glm::ivec3 indicesWithinPolygroup;
public:
  IndexedTriangle(glm::ivec3 indices);
  IndexedTriangle(glm::ivec3 indices, const std::vector<BufferedVertex> &arrayWithinPoly, BufferManager &bufferBoss);

  glm::ivec3 bufferIndices (const std::vector<BufferedVertex> &arrayWithinPoly) const;
  void addToBuffer(const std::vector<BufferedVertex> &arrayWithinPoly, BufferManager &bufferBoss);
};


class WeakSuperMesh {
  std::unique_ptr<BufferManager> boss;
  std::map<PolyGroupID, std::vector<BufferedVertex>> vertices = {};
  std::map<PolyGroupID, std::vector<IndexedTriangle>> triangles = {};

public:
  WeakSuperMesh();
  WeakSuperMesh(const std::vector<Vertex> &hardVertices, const std::vector<glm::ivec3> &faceIndices, PolyGroupID id);
  void addNewPolygroup(const std::vector<Vertex> &hardVertices, const std::vector<glm::ivec3> &faceIndices, PolyGroupID id);

  void* bufferIndexLocation() const { return boss->firstElementAddress(INDEX); }
  size_t bufferIndexSize() const { return boss->bufferSize(INDEX); }
  int bufferIndexLength() const { return boss->bufferLength(INDEX); }
  const void* getBufferLocation(CommonBufferType type) const { return boss->firstElementAddress(type); }
  unsigned int getBufferLength(CommonBufferType type) const { return boss->bufferLength(type); }
};

#endif