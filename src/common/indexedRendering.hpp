#pragma once

#include "renderingUtils.hpp"
// #include "src/geometry/smoothParametric.hpp"

#include <set>

struct buff4x4 {
    std::vector<glm::vec4> a, b, c, d;
};

inline glm::mat4 getMat(const buff4x4 &buff, int index) {
  return glm::mat4(buff.a[index], buff.b[index], buff.c[index], buff.d[index]);
}



#define buff2 std::vector<glm::vec2>
#define buff3 std::vector<glm::vec3>
#define buff4 std::vector<glm::vec4>
#define ibuff3 std::vector<glm::ivec3>

struct Stds {
    buff3 positions;
    buff3 normals;
    buff2 uvs;
    buff4 colors;
};

enum CommonBufferType {
  POSITION, NORMAL, UV, COLOR, MATERIAL1, MATERIAL2, MATERIAL3, MATERIAL4, INDEX, EXTRA0, EXTRA1, EXTRA2, EXTRA3, EXTRA4
};
const std::map<CommonBufferType, int> bufferTypeLength = {{POSITION, 3},  {NORMAL, 3},    {UV, 2},        {COLOR, 4},
                                                          {MATERIAL1, 4}, {MATERIAL2, 4}, {MATERIAL3, 4}, {MATERIAL4, 4},
                                                          {INDEX, 3},  {EXTRA0, 4},   {EXTRA1, 4},    {EXTRA2, 4}, {EXTRA3, 4}, {EXTRA4, 4}};

inline int bufferElementLength(CommonBufferType type) { return bufferTypeLength.at(type); }
inline size_t bufferElementSize(CommonBufferType type) { return type != INDEX ? bufferElementLength(type) * sizeof(FLOAT) : bufferElementLength(type) * sizeof(GLuint); }
class IndexedTriangle;


class BufferManager {
    std::unique_ptr<Stds> stds;
    std::unique_ptr<buff4> extra0;
    std::unique_ptr<buff4x4> mater, extra;
    std::unique_ptr<ibuff3> indices;
    std::set<CommonBufferType> activeBuffers;
    void insertValueToSingleBuffer(CommonBufferType type, void *valueAddress);
    void insertDefaultValueToSingleBuffer(CommonBufferType type);


public:
    BufferManager(const BufferManager &other);
    BufferManager & operator=(const BufferManager &other);


    BufferManager(BufferManager &&other) noexcept;
    BufferManager &operator=(BufferManager &&other) noexcept;
    explicit BufferManager(std::set<CommonBufferType> activeBuffers = {POSITION, NORMAL, UV, COLOR, INDEX});
    explicit BufferManager(bool materials,std::set<CommonBufferType> extras = {} );

    int bufferLength(CommonBufferType type) const;
    size_t bufferSize(CommonBufferType type) const { return bufferLength(type) * bufferElementSize(type); }
    void *firstElementAddress(CommonBufferType type) const;
    bool isActive(CommonBufferType type) const { return activeBuffers.contains(type); }
    bool hasMaterial() const { return isActive(MATERIAL1); }

    int addTriangleVertexIndices(glm::ivec3 ind);

    int addStdAttributesFromVertex(glm::vec3 pos, glm::vec3 norm, glm::vec2 uv, glm::vec4 col);
    int addMaterialBufferData(const MaterialPhong &mat);
    int addMaterialBufferData(glm::mat4 mat);
    int addMaterialBufferData(glm::vec4 ambient, glm::vec4 diffuse, glm::vec4 specular, glm::vec4 intencity);
    int addFullVertexData(glm::vec3 pos, glm::vec3 norm, glm::vec2 uv, glm::vec4 col, glm::mat4 mat, bool material=false);
    int addFullVertexData(const Vertex &v, bool material=false);

    void reserveSpace(int targetSize);
    void reserveSpaceForIndex(int targetSize) { indices->reserve(targetSize); }
    void reserveAdditionalSpace(int extraStorage) { reserveSpace(bufferLength(POSITION) + extraStorage); }
    void reserveAdditionalSpaceForIndex(int extraStorage) { reserveSpaceForIndex(bufferLength(INDEX) + extraStorage); }
    void initialiseExtraBufferSlot(int slot);

    glm::vec3 getPosition(int index) const { return stds->positions[index]; }
    glm::vec3 getNormal(int index) const { return stds->normals[index]; }
    glm::vec2 getUV(int index) const { return stds->uvs[index]; }
    glm::vec4 getColor(int index) const { return stds->colors[index]; }
    glm::mat4 getMaterial(int index) const { return getMat(*mater, index); }
    glm::vec4 getExtra(int index, int slot = 1) const;
    float getExtraSlot(int index, int slot = 1, int component = 3) const { return getExtra(index, slot)[component]; }

    void setPosition(int index, glm::vec3 value) { stds->positions[index] = value; }
    void setNormal(int index, glm::vec3 value) { stds->normals[index] = value; }
    void setUV(int index, glm::vec2 value) { stds->uvs[index] = value; }
    void setColor(int index, glm::vec4 value) { stds->colors[index] = value; }
    void setMaterial(int index, glm::mat4 value);

    void setExtra(int index, glm::vec4 value, int slot = 1);
    void setExtra(int index, glm::vec3 value, int slot = 1);
    void setExtra(int index, float value, int slot = 1, int component = 3);
};





class BufferedVertex {
  BufferManager &bufferBoss;
  int index;
public:
  BufferedVertex(BufferManager &bufferBoss, int index) : bufferBoss(bufferBoss), index(index) {}
  BufferedVertex(const BufferedVertex &other) : bufferBoss(other.bufferBoss), index(other.index) {}
  BufferedVertex(BufferedVertex &&other) noexcept : bufferBoss(other.bufferBoss), index(other.index) {}
  BufferedVertex(BufferManager &bufferBoss, const Vertex &v, bool material);

  int getIndex() const { return index; }
  glm::vec3 getPosition() const { return bufferBoss.getPosition(index); }
  glm::vec3 getNormal() const { return bufferBoss.getNormal(index); }
  glm::vec2 getUV() const { return bufferBoss.getUV(index); }
  glm::vec4 getColor() const { return bufferBoss.getColor(index); }
  glm::mat4 getMaterial() const { return bufferBoss.getMaterial(index); }
  glm::vec4 getExtra(int slot = 1) const { return bufferBoss.getExtra(index, slot); }
  Vertex getVertex() const { return Vertex(getPosition(), getUV(), getNormal(), getColor()); }

  void setPosition(glm::vec3 value) { bufferBoss.setPosition(index, value); }
  void setNormal(glm::vec3 value) { bufferBoss.setNormal(index, value); }
  void setUV(glm::vec2 value) { bufferBoss.setUV(index, value); }
  void setColor(glm::vec4 value) { bufferBoss.setColor(index, value); }
  void setMaterial(glm::mat4 value) { bufferBoss.setMaterial(index, value); }
  void setExtra(glm::vec4 value, int slot = 1) { bufferBoss.setExtra(index, value, slot); }
  void setExtra(glm::vec3 value, int slot = 1) { bufferBoss.setExtra(index, value, slot); }
  void setExtra(float value, int slot = 1, int component = 3) { bufferBoss.setExtra(index, value, slot, component); }
  void applyFunction(const SpaceEndomorphism &f);


};

class WeakSuperMesh;





class IndexedTriangle {
  int metaIndex;
  glm::ivec3 indicesWithinPolygroup;

public:
  explicit IndexedTriangle(glm::ivec3 indices);
  IndexedTriangle(glm::ivec3 indices, const std::vector<BufferedVertex> &arrayWithinPoly, BufferManager &bufferBoss);


  glm::ivec3 bufferIndices (const std::vector<BufferedVertex> &arrayWithinPoly) const;
  void addToBuffer(const std::vector<BufferedVertex> &arrayWithinPoly, BufferManager &bufferBoss);
  glm::ivec3 getLocalIndices() const { return indicesWithinPolygroup; }
  glm::mat2 barMatrix(const std::vector<BufferedVertex> &arrayWithinPoly) const;
  glm::mat3 orthonormalFrame(const std::vector<BufferedVertex> &arrayWithinPoly) const;
  glm::vec3 fromPlanar(glm::vec2, const std::vector<BufferedVertex> &arrayWithinPoly) const;
  glm::vec2 toPlanar(glm::vec3, const std::vector<BufferedVertex> &arrayWithinPoly) const;
  glm::vec3 fromBars(glm::vec2, const std::vector<BufferedVertex> &arrayWithinPoly) const;
  glm::vec2 toBars(glm::vec3, const std::vector<BufferedVertex> &arrayWithinPoly) const;
  std::array<glm::vec3, 3> borderTriangle(float width, const std::vector<BufferedVertex> &arrayWithinPoly) const;
  Vertex getVertex(int i, const std::vector<BufferedVertex> &arrayWithinPoly) const { return arrayWithinPoly[indicesWithinPolygroup[i]].getVertex(); }
  glm::vec3 faceNormal(const std::vector<BufferedVertex> &arrayWithinPoly) const { return glm::normalize(glm::cross(getVertex(1, arrayWithinPoly).getPosition() - getVertex(0, arrayWithinPoly).getPosition(), getVertex(2, arrayWithinPoly).getPosition() - getVertex(0, arrayWithinPoly).getPosition())); }
};



class SmoothParametricSurface;

class WeakSuperMesh {
  std::unique_ptr<BufferManager> boss;
  std::map<PolyGroupID, std::vector<BufferedVertex>> vertices = {};
  std::map<PolyGroupID, std::vector<IndexedTriangle>> triangles = {};
  std::shared_ptr<MaterialPhong> material = nullptr;

public:

  WeakSuperMesh();
  WeakSuperMesh(const std::vector<Vertex> &hardVertices, const std::vector<glm::ivec3> &faceIndices, const PolyGroupID &id);
  WeakSuperMesh(const char* filename, const PolyGroupID &id);

  void addNewPolygroup(const std::vector<Vertex> &hardVertices, const std::vector<glm::ivec3> &faceIndices, const PolyGroupID &id);
  void addNewPolygroup(const char* filename, const PolyGroupID &id);



  WeakSuperMesh & operator=(const WeakSuperMesh &other);
  WeakSuperMesh(WeakSuperMesh &&other) noexcept;
  WeakSuperMesh& operator=(WeakSuperMesh &&other) noexcept;
  WeakSuperMesh(const WeakSuperMesh &other);


  WeakSuperMesh(const SmoothParametricSurface &surf, int tRes, int uRes, const PolyGroupID &id);
  void addUniformSurface(const SmoothParametricSurface &surf, int tRes, int uRes, const PolyGroupID &id);

  void* bufferIndexLocation() const { return boss->firstElementAddress(INDEX); }
  size_t bufferIndexSize() const { return boss->bufferSize(INDEX); }
  int bufferIndexLength() const { return boss->bufferLength(INDEX); }
  const void* getBufferLocation(CommonBufferType type) const { return boss->firstElementAddress(type); }
  unsigned int getBufferLength(CommonBufferType type) const { return boss->bufferLength(type); }
  size_t getBufferSize(CommonBufferType type) const { return boss->bufferSize(type); }
  std::vector<PolyGroupID> getPolyGroupIDs() const;
  BufferManager& getBufferBoss() const { return *boss; }
  bool isActive(CommonBufferType type) const { return boss->isActive(type); }
  bool hasGlobalTextures() const { return !isActive(MATERIAL1) && material->textured(); }

  void deformPerVertex(const PolyGroupID &id, const std::function<void(BufferedVertex&)> &deformation) { for (auto &v : vertices.at(id)) deformation(v);  }
  void moveAlongVectorField(const PolyGroupID &id, VectorFieldR3 X, float delta=1);
  void deformWithAmbientMap(const PolyGroupID &id, SpaceEndomorphism f);
  void initGlobalTextures() {if (hasGlobalTextures()) material->initTextures();}

  void affineTransform(const glm::mat3 &M, glm::vec3 v, const PolyGroupID &id) {deformWithAmbientMap(id, SpaceEndomorphism::affine(M, v));}
  void affineTransform(const glm::mat3 &M, glm::vec3 v) {for (auto& name: getPolyGroupIDs()) affineTransform(M, v, name);}
  void shift(glm::vec3 v, const PolyGroupID &id) { affineTransform(glm::mat3(1), v, id); }
  void shift(glm::vec3 v) { affineTransform(glm::mat3(1), v); }
  void scale(float s, const PolyGroupID &id) { affineTransform(glm::mat3(s), glm::vec3(0), id); }
  void scale(float s) { affineTransform(glm::mat3(s), glm::vec3(0)); }
  void addGlobalMaterial(const MaterialPhong &mat) { material = std::make_shared<MaterialPhong>(mat); }



  WeakSuperMesh subdivideBarycentric(const PolyGroupID &id) const;
  WeakSuperMesh subdivideEdgecentric(const PolyGroupID &id) const;
  WeakSuperMesh wireframe(PolyGroupID id, PolyGroupID targetId, float width, float heightCenter, float heightSide) const;

  std::vector<Vertex> getVertices(const PolyGroupID &id) const;
  std::vector<BufferedVertex> getBufferedVertices(const PolyGroupID &id) const { return vertices.at(id); }
  std::vector<glm::ivec3> getIndices(const PolyGroupID &id) const;
  std::vector<IndexedTriangle> getTriangles(const PolyGroupID &id) const { return triangles.at(id); }
  glm::vec4 getIntencities() const { return material->compressIntencities(); }
  MaterialPhong getMaterial() const { return *material; }
};
std::function<void(float, float)> deformationOperator (const std::function<void(BufferedVertex&, float, float)> &deformation, WeakSuperMesh &mesh, const PolyGroupID &id);
std::function<void(float)> deformationOperator (const std::function<void(BufferedVertex&, float)> &deformation, WeakSuperMesh &mesh, const PolyGroupID &id);

std::function<void(float, float)> moveAlongCurve(const SmoothParametricCurve &curve, WeakSuperMesh &mesh, const PolyGroupID &id);
