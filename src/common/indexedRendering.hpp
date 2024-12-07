#pragma once

#include "renderingUtils.hpp"
// #include "src/geometry/smoothParametric.hpp"

#include <set>




struct buff4x4 {
    std::vector<vec4> a, b, c, d;
};

inline mat4 getMat(const buff4x4 &buff, int index) {
  return mat4(buff.a[index], buff.b[index], buff.c[index], buff.d[index]);
}





struct Stds {
    BUFF3 positions;
    BUFF3 normals;
    BUFF2 uvs;
    BUFF4 colors;
};

enum CommonBufferType {
  POSITION, NORMAL, UV, COLOR, MATERIAL1, MATERIAL2, MATERIAL3, MATERIAL4, INDEX, EXTRA0, EXTRA1, EXTRA2, EXTRA3, EXTRA4
};
const std::map<CommonBufferType, int> bufferTypeLength = {{POSITION, 3},  {NORMAL, 3},    {UV, 2},        {COLOR, 4},
                                                          {MATERIAL1, 4}, {MATERIAL2, 4}, {MATERIAL3, 4}, {MATERIAL4, 4},
                                                          {INDEX, 3},  {EXTRA0, 4},   {EXTRA1, 4},    {EXTRA2, 4}, {EXTRA3, 4}, {EXTRA4, 4}};

inline int bufferElementLength(CommonBufferType type) { return bufferTypeLength.at(type); }
inline size_t bufferElementSize(CommonBufferType type) { return type != INDEX ? bufferElementLength(type) * sizeof(FLOAT) : bufferElementLength(type) * sizeof(GLuint); }


class BufferManager {
    std::unique_ptr<Stds> stds;
    std::unique_ptr<BUFF4> extra0;
    std::unique_ptr<buff4x4> extra;
    std::unique_ptr<IBUFF3> indices;
    std::set<CommonBufferType> activeBuffers;
    void insertValueToSingleBuffer(CommonBufferType type, void *valueAddress);
    void insertDefaultValueToSingleBuffer(CommonBufferType type);


public:
    BufferManager(const BufferManager &other);
    BufferManager & operator=(const BufferManager &other);


    BufferManager(BufferManager &&other) noexcept;
    BufferManager &operator=(BufferManager &&other) noexcept;
    explicit BufferManager(const std::set<CommonBufferType> &activeBuffers = {POSITION, NORMAL, UV, COLOR, INDEX});
    explicit BufferManager(bool materials, const std::set<CommonBufferType> &extras = {} );

    int bufferLength(CommonBufferType type) const;
    size_t bufferSize(CommonBufferType type) const { return bufferLength(type) * bufferElementSize(type); }
    void *firstElementAddress(CommonBufferType type) const;
    bool isActive(CommonBufferType type) const { return activeBuffers.contains(type); }
    bool hasMaterial() const { return isActive(MATERIAL1); }

    int addTriangleVertexIndices(glm::ivec3 ind, int shift = 0);

    int addStdAttributesFromVertex(vec3 pos, vec3 norm, vec2 uv, vec4 col);
    int addMaterialBufferData(const MaterialPhong &mat);
    int addFullVertexData(vec3 pos, vec3 norm, vec2 uv, vec4 col);
    int addFullVertexData(const Vertex &v);

    void reserveSpace(int targetSize);
    void reserveSpaceForIndex(int targetSize) { indices->reserve(targetSize); }
    void reserveAdditionalSpace(int extraStorage) { reserveSpace(bufferLength(POSITION) + extraStorage); }
    void reserveAdditionalSpaceForIndex(int extraStorage) { reserveSpaceForIndex(bufferLength(INDEX) + extraStorage); }
    void initialiseExtraBufferSlot(int slot);

    vec3 getPosition(int index) const { return stds->positions[index]; }
    vec3 getNormal(int index) const { return stds->normals[index]; }
    vec2 getUV(int index) const { return stds->uvs[index]; }
    vec4 getColor(int index) const { return stds->colors[index]; }
    vec4 getExtra(int index, int slot = 1) const;
    float getExtraSlot(int index, int slot = 1, int component = 3) const { return getExtra(index, slot)[component]; }
    glm::ivec3 getFaceIndices(int index) const { return (*indices)[index]; }
    Vertex getVertex(int index) const { return Vertex(getPosition(index), getUV(index), getNormal(index), getColor(index)); }

    void setPosition(int index, vec3 value) { stds->positions[index] = value; }
    void setNormal(int index, vec3 value) { stds->normals[index] = value; }
    void setUV(int index, vec2 value) { stds->uvs[index] = value; }
    void setColor(int index, vec4 value) { stds->colors[index] = value; }
    void setColor(int index, float value, int component) { stds->colors[index][component] = value; }
    void setMaterial(int index, mat4 value);
	void setFaceIndices(int index, const ivec3 &in) { indices->at(index) = in; }


    void setExtra(int index, vec4 value, int slot = 1);
    void setExtra(int index, vec3 value, int slot = 1);
    void setExtra(int index, float value, int slot = 1, int component = 3);
};





class BufferedVertex {
  BufferManager &bufferBoss;
  int index;
public:
  BufferedVertex(BufferManager &bufferBoss, int index) : bufferBoss(bufferBoss), index(index) {}
  BufferedVertex(const BufferedVertex &other) = default;
  BufferedVertex(BufferedVertex &&other) noexcept : bufferBoss(other.bufferBoss), index(other.index) {}
  BufferedVertex(BufferManager &bufferBoss, const Vertex &v);

  int getIndex() const { return index; }
  vec3 getPosition() const { return bufferBoss.getPosition(index); }
  vec3 getNormal() const { return bufferBoss.getNormal(index); }
  vec2 getUV() const { return bufferBoss.getUV(index); }
  vec4 getColor() const { return bufferBoss.getColor(index); }
  vec4 getExtra(int slot = 1) const { return bufferBoss.getExtra(index, slot); }
  Vertex getVertex() const { return Vertex(getPosition(), getUV(), getNormal(), getColor()); }
	vec2 getSurfaceParams() const {return vec2(getColor());}

  void setPosition(vec3 value) { bufferBoss.setPosition(index, value); }
  void setNormal(vec3 value) { bufferBoss.setNormal(index, value); }
  void setUV(vec2 value) { bufferBoss.setUV(index, value); }
  void setColor(vec4 value) { bufferBoss.setColor(index, value); }
  void setColor(float value, int i) { bufferBoss.setColor(index, value, i); }
  void setMaterial(const mat4 &value) { bufferBoss.setMaterial(index, value); }
  void setExtra(vec4 value, int slot = 1) { bufferBoss.setExtra(index, value, slot); }
  void setExtra(vec3 value, int slot = 1) { bufferBoss.setExtra(index, value, slot); }
  void setExtra(float value, int slot = 1, int component = 3) { bufferBoss.setExtra(index, value, slot, component); }
  void applyFunction(const SpaceEndomorphism &f);
  void setVertex(const Vertex &v) { setPosition(v.getPosition()); setUV(v.getUV()); setNormal(v.getNormal()); setColor(v.getColor()); }

};

class WeakSuperMesh;





class IndexedTriangleDeprecated {
  int metaIndex;
  glm::ivec3 indicesWithinPolygroup;

public:
  explicit IndexedTriangleDeprecated(glm::ivec3 indices);
  IndexedTriangleDeprecated(glm::ivec3 indices, const std::vector<BufferedVertex> &arrayWithinPoly, BufferManager &bufferBoss);


  glm::ivec3 bufferIndices (const std::vector<BufferedVertex> &arrayWithinPoly) const;
  void addToBuffer(const std::vector<BufferedVertex> &arrayWithinPoly, BufferManager &bufferBoss);
  glm::ivec3 getLocalIndices() const { return indicesWithinPolygroup; }
  mat2 barMatrix(const std::vector<BufferedVertex> &arrayWithinPoly) const;
  mat3 orthonormalFrame(const std::vector<BufferedVertex> &arrayWithinPoly) const;
  vec3 fromPlanar(vec2, const std::vector<BufferedVertex> &arrayWithinPoly) const;
  vec2 toPlanar(vec3, const std::vector<BufferedVertex> &arrayWithinPoly) const;
  vec3 fromBars(vec2, const std::vector<BufferedVertex> &arrayWithinPoly) const;
  vec2 toBars(vec3, const std::vector<BufferedVertex> &arrayWithinPoly) const;
  std::array<vec3, 3> borderTriangle(float width, const std::vector<BufferedVertex> &arrayWithinPoly) const;
  Vertex getVertex(int i, const std::vector<BufferedVertex> &arrayWithinPoly) const { return arrayWithinPoly[indicesWithinPolygroup[i]].getVertex(); }
  vec3 faceNormal(const std::vector<BufferedVertex> &arrayWithinPoly) const { return normalize(cross(getVertex(1, arrayWithinPoly).getPosition() - getVertex(0, arrayWithinPoly).getPosition(), getVertex(2, arrayWithinPoly).getPosition() - getVertex(0, arrayWithinPoly).getPosition())); }
};



class IndexedTriangle {
    int index;
    BufferManager &bufferBoss;

public:
//    IndexedTriangle(BufferManager &bufferBoss, int index) : index(index), bufferBoss(bufferBoss) {}
    IndexedTriangle(const IndexedTriangle &other) : index(other.index), bufferBoss(other.bufferBoss) {}
    IndexedTriangle(IndexedTriangle &&other) noexcept : index(other.index), bufferBoss(other.bufferBoss) {}
    IndexedTriangle(BufferManager &bufferBoss, glm::ivec3 index, int shift) : index(bufferBoss.addTriangleVertexIndices(index, shift)), bufferBoss(bufferBoss) {}

    ivec3 getVertexIndices() const { return bufferBoss.getFaceIndices(index); }
    Vertex getVertex(int i) const { return bufferBoss.getVertex(getVertexIndices()[i]); }
    mat2 barMatrix() const;
    mat3 orthonormalFrame() const;
    vec3 fromPlanar(vec2 v) const;
    vec2 toPlanar(vec3 v) const;
    vec3 fromBars(vec2 v) const;
    vec2 toBars(vec3 v) const;
    std::array<vec3, 3> borderTriangle(float width) const;
    vec3 faceNormal() const { return normalize(cross(getVertex(1).getPosition() - getVertex(0).getPosition(), getVertex(2).getPosition() - getVertex(0).getPosition())); }
    vec3 center() const { return (getVertex(0).getPosition() + getVertex(1).getPosition() + getVertex(2).getPosition()) / 3.f; }
    float area() const { return 0.5f * length(cross(getVertex(1).getPosition() - getVertex(0).getPosition(), getVertex(2).getPosition() - getVertex(0).getPosition())); }

	bool containsEdge(int i, int j) const;
	bool containsEdge(ivec2 edge) const { return containsEdge(edge.x, edge.y); }
	void setVertexIndices(const ivec3 &in) { bufferBoss.setFaceIndices(index, in); }
	void changeOrientation() { bufferBoss.setFaceIndices(index, ivec3(getVertexIndices().z, getVertexIndices().y, getVertexIndices().x)); }
};



class SmoothParametricSurface;

class WeakSuperMesh {
protected:
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
	WeakSuperMesh(const SmoothParametricSurface &surf, int tRes, int uRes) : WeakSuperMesh(surf, tRes, uRes, randomID()	) {}
  void addUniformSurface(const SmoothParametricSurface &surf, int tRes, int uRes, const PolyGroupID &id);
	void addUniformSurface(const SmoothParametricSurface &surf, int tRes, int uRes) {return addUniformSurface(surf, tRes, uRes, randomID());}
  void merge (const WeakSuperMesh &other);
	void mergeAndKeepID(const WeakSuperMesh &other);

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
  BufferedVertex& getAnyVertexFromPolyGroup(const PolyGroupID &id) { return vertices.at(id).front(); }

  void deformPerVertex(const PolyGroupID &id, const std::function<void(BufferedVertex&)> &deformation) { for (auto &v : vertices.at(id)) deformation(v);  }
  void deformPerVertex(const std::function<void(BufferedVertex&)> &deformation) { for (auto id: getPolyGroupIDs()) for (auto &v : vertices.at(id)) deformation(v);  }
  void deformPerVertex(const PolyGroupID &id, const std::function<void(int, BufferedVertex&)> &deformation);
  void deformPerId(const std::function<void(BufferedVertex&, PolyGroupID)> &deformation) { for (auto id: getPolyGroupIDs()) for (auto &v : vertices.at(id)) deformation(v, id);  }

  vec2 getSurfaceParameters(const BufferedVertex &v) const;
  void encodeSurfacePoint(BufferedVertex &v, const SmoothParametricSurface &surf, vec2 tu);
  void adjustToNewSurface(const SmoothParametricSurface &surf, const PolyGroupID &id);
  void adjustToNewSurface(const SmoothParametricSurface &surf);

  void moveAlongVectorField(const PolyGroupID &id, VectorFieldR3 X, float delta=1);
  void deformWithAmbientMap(const PolyGroupID &id, SpaceEndomorphism f);
  void deformWithAmbientMap(const SpaceEndomorphism &f) { for (auto id: getPolyGroupIDs()) deformWithAmbientMap(id, f); }
  void initGlobalTextures() {if (hasGlobalTextures()) material->initTextures();}

  void affineTransform(const mat3 &M, vec3 v, const PolyGroupID &id) {deformWithAmbientMap(id, SpaceEndomorphism::affine(M, v));}
  void affineTransform(const mat3 &M, vec3 v) {for (auto& name: getPolyGroupIDs()) affineTransform(M, v, name);}
  void shift(vec3 v, const PolyGroupID &id) { affineTransform(mat3(1), v, id); }
  void shift(vec3 v) { affineTransform(mat3(1), v); }
  void scale(float s, const PolyGroupID &id) { affineTransform(mat3(s), vec3(0), id); }
  void scale(float s) { affineTransform(mat3(s), vec3(0)); }
  void addGlobalMaterial(const MaterialPhong &mat) { material = std::make_shared<MaterialPhong>(mat); }

  void flipNormals(const PolyGroupID &id) { deformPerVertex(id, [](BufferedVertex &v) { v.setNormal(-v.getNormal()); }); }
  void flipNormals() { for (auto &name: getPolyGroupIDs()) flipNormals(name); }
  void pointNormalsInDirection(vec3 dir, const PolyGroupID &id);
  void pointNormalsInDirection(vec3 dir) { for (auto &name: getPolyGroupIDs()) pointNormalsInDirection(dir, name); }


  WeakSuperMesh subdivideBarycentric(const PolyGroupID &id) const;
  WeakSuperMesh subdivideEdgecentric(const PolyGroupID &id) const;
  WeakSuperMesh wireframe(PolyGroupID id, PolyGroupID targetId, float width, float heightCenter, float heightSide) const;

  std::vector<Vertex> getVertices(const PolyGroupID &id) const;
  std::vector<BufferedVertex> getBufferedVertices(const PolyGroupID &id) const { return vertices.at(id); }
  std::vector<glm::ivec3> getIndices(const PolyGroupID &id) const;
  std::vector<IndexedTriangle> getTriangles(const PolyGroupID &id) const { return triangles.at(id); }
  vec4 getIntencities() const { return material->compressIntencities(); }
  MaterialPhong getMaterial() const { return *material; }

	vector<int> findVertexNeighbours(int i, const PolyGroupID &id) const;
	vector<int> findVertexParentTriangles(int i, const PolyGroupID &id) const;
	void recalculateNormal(int i, const PolyGroupID &id);
	void recalculateNormalsNearby(int i, const PolyGroupID &id);
	void recalculateNormals(const PolyGroupID &id);
	void recalculateNormals();
	void orientFaces(const PolyGroupID &id);
	void orientFaces();


  vector<int> findNeighboursSorted(int i, const PolyGroupID &id) const;
	bool checkIfHasCompleteNeighbourhood(int i, const PolyGroupID &id) const;
	float meanCurvature(int i, const PolyGroupID &id) const;
	vec3 meanCurvatureVector(int i, const PolyGroupID &id) const;
	void meanCurvatureFlowDeform(float dt, const PolyGroupID &id);
	float GaussCurvature(int i, const PolyGroupID &id) const;
	void paintMeanCurvature(const PolyGroupID &id);
	void paintMeanCurvature();

  template<typename T>
    T integrateOverTriangles(const std::function<T(const IndexedTriangle &)> &f, PolyGroupID id) const;

    vec3 centerOfMass(PolyGroupID id) const;
	vec3 centerOfMass() const;

	mat3 inertiaTensorCMAppBd(PolyGroupID id) const;
	mat3 inertiaTensorAppBd(PolyGroupID id, vec3 p) const;
};


WeakSuperMesh _wireframe(const SmoothParametricSurface &s, float width, int n, int m, int curve_res_rad, int curve_res_hor);




std::function<void(float, float)> deformationOperator (const std::function<void(BufferedVertex&, float, float)> &deformation, WeakSuperMesh &mesh, const PolyGroupID &id);
std::function<void(float)> deformationOperator (const std::function<void(BufferedVertex&, float)> &deformation, WeakSuperMesh &mesh, const PolyGroupID &id);
std::function<void(float, float)> moveAlongCurve(const SmoothParametricCurve &curve, WeakSuperMesh &mesh, const PolyGroupID &id);

template<typename T>
T WeakSuperMesh::integrateOverTriangles(const std::function<T(const IndexedTriangle &)> &f, PolyGroupID id) const {
    T sum = T(0);
    for (const IndexedTriangle &t: triangles.at(id))
        sum += f(t)*t.area();
    return sum;
}


class Wireframe : public WeakSuperMesh {
	SmoothParametricSurface surf;
	float width;
	int n, m, curve_res_rad, curve_res_hor;
public:
	Wireframe(const SmoothParametricSurface &surf, float width, int n, int m, int curve_res_rad, int curve_res_hor);

	void changeBaseSurface(const SmoothParametricSurface &newsurf);
	vec2 getSurfaceParameters(const BufferedVertex &v) const;
};
