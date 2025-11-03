#pragma once

#include "renderingUtils.hpp"
#include "pde.hpp"
#include "randomUtils.hpp"
#include "SDFObjects.hpp"


struct Stds {
	vector<vec3> positions;
	vector<vec3> normals;
	vector<vec2> uvs;
	vector<vec4> colors;
};

namespace BufferNames {
	constexpr string COL_BUF = "color";
	constexpr string POS_BUF = "position";
	constexpr string NORM_BUF = "normal";
	constexpr string UV_BUF = "uv";
	constexpr string IND_BUF = "index";
	const vector RESERVED_BUF = {COL_BUF, POS_BUF, NORM_BUF, UV_BUF, IND_BUF};
}




class BufferManager {
	unique_ptr<Stds> stds;
	unique_ptr<vector<vec4>> extra0;
	unique_ptr<vector<vec4>> extra1;
	unique_ptr<vector<vec4>> extra2;
	unique_ptr<vector<vec4>> extra3;
	unique_ptr<vector<vec4>> extra4;
	unique_ptr<vector<ivec3>> indices;
	vector<string> extraBufferNames;


	unique_ptr<vector<vec4>>& extraBufferPtrAt(int slot);
	vec4 extraElementAt(int index, int slot) const;
	unique_ptr<vector<vec4>>& extraBufferPtrByName(const string& name);
	void setElementAt(int index, float value, int slot, int component);
	void setElementAt(int index, vec4 value, int slot);

public:
	explicit BufferManager(const vector<string>& extra_names = {});
	explicit BufferManager(const string& extra_name);
	explicit BufferManager(const string& extra_name0, const string& extra_name1);
	explicit BufferManager(const string& extra_name0, const string& extra_name1, const string& extra_name2);
	explicit BufferManager(const string& extra_name0, const string& extra_name1, const string& extra_name2, const string& extra_name3);
	explicit BufferManager(const string& extra_name0, const string& extra_name1, const string& extra_name2, const string& extra_name3, const string& extra_name4);

	BufferManager(const BufferManager& other);
	BufferManager& operator=(const BufferManager& other);
	BufferManager(BufferManager&& other) noexcept;
	BufferManager& operator=(BufferManager&& other) noexcept;

	vs_dim getAttributeDimension(const string& attributeName) const;
	array_len attributeDataLength() const;
	array_len indexDataLength() const;
	byte_size attributeDataSize(const string& attributeName) const;
	byte_size indexDataSize() const;
	raw_data_ptr attributeDataPtr(const string& attributeName) const;
	raw_data_ptr indexDataPtr() const;

	bool isActive(const string& attributeName) const;
	bool isActive(int i) const;
	string getExtraBufferName(int slot) const;
	vector<string> getExtraBufferNames() const;
	vs_dim numberOfExtraBuffers() const;

	int indexOfExtraBuffer(const string& name) const;

	int addTriangleVertexIndices(ivec3 ind, int shift = 0) const;
	int addFullVertexData(const Vertex& v) const;

	void reserveSpace(int targetSize) const;
	void reserveSpaceForIndex(int targetSize) const;
	void reserveAdditionalSpace(int extraStorage) const;
	void reserveAdditionalSpaceForIndex(int extraStorage);

	vec3 getPosition(int index) const;
	vec3 getNormal(int index) const;
	vec2 getUV(int index) const;
	vec4 getColor(int index) const;
	vec4 getExtra(int index, const string& name) const;

	ivec3 getFaceIndices(int index) const;
	Vertex getVertex(int index) const;

	void setPosition(int index, vec3 value);
	void setNormal(int index, vec3 value);
	void setUV(int index, vec2 value);
	void setUV(int index, float value, int coord);
	void setColor(int index, vec4 value);
	void setColor(int index, float value, int component);
	void setExtra(int index, vec4 value, const string& name);
	void setExtra(int index, vec3 value, const string& name);
	void setExtra(int index, float value, const string& name, int component);
	void setFaceIndices(int index, const ivec3& in);

};


class BufferedVertex {
	BufferManager& bufferBoss;
	int index;

public:
	BufferedVertex();
	BufferedVertex(BufferManager& bufferBoss, int index);

	BufferedVertex(const BufferedVertex& other);
	BufferedVertex(BufferedVertex&& other) noexcept;
	BufferedVertex(BufferManager& bufferBoss, const Vertex& v);
	BufferedVertex& operator=(const BufferedVertex& other);
	BufferedVertex& operator=(BufferedVertex&& other) noexcept;

	int getIndex() const;
	vec3 getPosition() const;
	vec3 getNormal() const;
	vec2 getUV() const;
	vec4 getColor() const;
	vec4 getExtra(const string& slot) const;
	Vertex getVertex() const;

	void setPosition(vec3 value) const;
	void setNormal(vec3 value) const;
	void setUV(vec2 value) const;
	void setUV(float value, int i) const;
	void setColor(vec4 value) const;
	void setColor(float value, int i) const;
	void setExtra(vec4 value, const string& name);
	void setExtra(vec3 value, const string& name);
	void setExtra(float value, const string& name, int component);
	void setVertex(const Vertex& v);

	void applyFunction(const SpaceEndomorphism& f);

};


class IndexedTriangle {
	int index;
	BufferManager& bufferBoss;

public:
	//    IndexedTriangle(BufferManager &bufferBoss, int index) : index(index), bufferBoss(bufferBoss) {}
	IndexedTriangle(const IndexedTriangle& other);
	IndexedTriangle(IndexedTriangle&& other) noexcept;
	IndexedTriangle(BufferManager& bufferBoss, ivec3 index, int shift);
	IndexedTriangle& operator=(const IndexedTriangle& other);
	IndexedTriangle& operator=(IndexedTriangle&& other) noexcept;

	ivec3 getVertexIndices() const;
	Vertex getVertex(int i) const;
	mat2 barMatrix() const;
	mat3 orthonormalFrame() const;
	vec3 fromPlanar(vec2 v) const;
	vec2 toPlanar(vec3 v) const;
	vec3 fromBars(vec2 v) const;
	vec2 toBars(vec3 v) const;
	array<vec3, 3> borderTriangle(float width) const;
	vec3 faceNormal() const;
	vec3 center() const;
	float area() const;

	bool containsEdge(int i, int j) const;
	bool containsEdge(ivec2 edge) const;
	void setVertexIndices(const ivec3& in) const;
	void changeOrientation() const;
};

class SmoothParametricSurface;


class IndexedMesh {
protected:
	unique_ptr<BufferManager> boss;
	unordered_map<PolyGroupID, int> polygroupIndexOrder = {};
	vector<vector<BufferedVertex>> vertices = {};
	vector<vector<IndexedTriangle>> triangles = {};

public:
	virtual ~IndexedMesh() = default;

	IndexedMesh();
	IndexedMesh(const vector<Vertex>& hardVertices, const vector<ivec3>& faceIndices, const PolyGroupID& id);
	IndexedMesh(const char* filename, const PolyGroupID& id);
	IndexedMesh(const SmoothParametricSurface& surf, int tRes, int uRes, const PolyGroupID& id = randomID());

	void addNewPolygroup(const vector<Vertex>& hardVertices, const vector<ivec3>& faceIndices, const PolyGroupID& id);
	void addNewPolygroup(const char* filename, const PolyGroupID& id);

	IndexedMesh& operator=(const IndexedMesh& other);
	IndexedMesh(IndexedMesh&& other) noexcept;
	IndexedMesh& operator=(IndexedMesh&& other) noexcept;
	IndexedMesh(const IndexedMesh& other);

	void addUniformSurface(const SmoothParametricSurface& surf, int tRes, int uRes, const PolyGroupID& id = randomID());

	void merge(const IndexedMesh& other);
	void mergeAndKeepID(const IndexedMesh& other);
	void copyPolygroup(const IndexedMesh& other, const PolyGroupID& id, const PolyGroupID& newId);
	void copyPolygroup(const PolyGroupID& id, const PolyGroupID& newId);

	vector<string> getActiveExtraBuffers() const;

	raw_data_ptr bufferIndexLocation() const;
	raw_data_ptr getBufferLocation(const string& name) const;

	byte_size bufferIndexSize() const;
	array_len bufferIndexLength() const;
	array_len getBufferLength() const;
	size_t getBufferSize(const string& name) const;

	vector<PolyGroupID> getPolyGroupIDs() const;
	BufferManager& getBufferBoss() const;

	BufferedVertex& getAnyVertexFromPolyGroup(const PolyGroupID& id);

	void deformPerVertex(const PolyGroupID& id, const HOM(BufferedVertex&, void)& deformation);
	void deformPerVertex(const HOM(BufferedVertex&, void)& deformation);
	void deformPerVertex(const PolyGroupID& id, const std::function<void(int, BufferedVertex&)>& deformation);
	void deformPerId(const std::function<void(BufferedVertex&, PolyGroupID)>& deformation);

	static vec2 getSurfaceParameters(const BufferedVertex& v);
	static void encodeSurfacePoint(BufferedVertex& v, const SmoothParametricSurface& surf, vec2 tu);
	void adjustToNewSurface(const SmoothParametricSurface& surf, const PolyGroupID& id);
	void adjustToNewSurface(const SmoothParametricSurface& surf);

	void moveAlongVectorField(const PolyGroupID& id, const VectorField& X, float delta = 1);
	void deformWithAmbientMap(const PolyGroupID& id, const SpaceEndomorphism& f);
	void deformWithAmbientMap(const SpaceEndomorphism& f);


	void affineTransform(const mat3& M, vec3 v, const PolyGroupID& id);
	void affineTransform(const mat3& M, vec3 v);
	void shift(vec3 v, const PolyGroupID& id);
	void shift(vec3 v);
	void scale(float s, const PolyGroupID& id);
	void scale(float s);

	void flipNormals(const PolyGroupID& id);
	void flipNormals();
	void pointNormalsInDirection(vec3 dir, const PolyGroupID& id);
	void pointNormalsInDirection(vec3 dir);

	IndexedMesh subdivideBarycentric(const PolyGroupID& id) const;
	IndexedMesh subdivideEdgecentric(const PolyGroupID& id) const;
	IndexedMesh wireframe(PolyGroupID id, PolyGroupID targetId, float width, float heightCenter, float heightSide) const;

	vector<Vertex> getVertices(const PolyGroupID& id) const;
	vector<BufferedVertex> getBufferedVertices(const PolyGroupID& id) const;
	vector<ivec3> getIndices(const PolyGroupID& id) const;
	vector<IndexedTriangle> getTriangles(const PolyGroupID& id) const;


	vector<int> findVertexNeighbours(int i, const PolyGroupID& id) const;
	vector<int> findVertexParentTriangles(int i, const PolyGroupID& id) const;
	void recalculateNormal(int i, const PolyGroupID& id);
	void recalculateNormalsNearby(int i, const PolyGroupID& id);
	void recalculateNormals(const PolyGroupID& id);
	void recalculateNormals();
	void orientFaces(const PolyGroupID& id);
	void orientFaces();

	vector<int> findNeighboursSorted(int i, const PolyGroupID& id) const;
	bool checkIfHasCompleteNeighbourhood(int i, const PolyGroupID& id) const;
	float meanCurvature(int i, const PolyGroupID& id) const;
	vec3 meanCurvatureVector(int i, const PolyGroupID& id) const;
	float GaussCurvature(int i, const PolyGroupID& id) const;

	template <typename T>
	T integrateOverTriangles(const HOM(const IndexedTriangle&, T)& f, PolyGroupID id) const;

	vec3 centerOfMass(PolyGroupID id) const;
	vec3 centerOfMass() const;

	mat3 inertiaTensorCMAppBd(PolyGroupID id) const;
	mat3 inertiaTensorAppBd(PolyGroupID id, vec3 p) const;
};


std::function<void(float, float)> deformationOperator(const std::function<void(BufferedVertex&, float, float)>& deformation, IndexedMesh& mesh, const PolyGroupID& id);
std::function<void(float)> deformationOperator(const std::function<void(BufferedVertex&, float)>& deformation, IndexedMesh& mesh, const PolyGroupID& id);
std::function<void(float, float)> moveAlongCurve(const SmoothParametricCurve& curve, IndexedMesh& mesh, const PolyGroupID& id);

template <typename T>
T IndexedMesh::integrateOverTriangles(const std::function<T(const IndexedTriangle&)>& f, PolyGroupID id) const {
	T sum = T(0);
	for (const IndexedTriangle& t : triangles.at(polygroupIndexOrder.at(id)))
		sum += f(t) * t.area();
	return sum;
}


class Wireframe : public IndexedMesh {
	SmoothParametricSurface surf;
	float width;
	int n, m, curve_res_rad, curve_res_hor;

public:
	Wireframe(const SmoothParametricSurface& surf, float width, int n, int m, int curve_res_rad, int curve_res_hor);
	void changeBaseSurface(const SmoothParametricSurface& newsurf);
	vec2 getSurfaceParameters(const BufferedVertex& v) const;
};


class PlanarFlowLines : public IndexedMesh {
	VectorFieldR2 X;
	float dt;
	int steps;
	std::function<float(float, float, float, vec2, vec2)> width; //w(t, t0, speed, x, x0)
	std::function<vec4(float, float, float, vec2, vec2)> color;;
	vector<vec2> startPoints;
	vector<float> startTimes;
	vector<PolyGroupID> ids;

	// uv = (t, w)
	// pos = (x, y, 0)
	// n = (0, 0, 1)
	// col = (color, speed)
	// extra1 = (t0, x0, y0, len)

public:
	PlanarFlowLines(const VectorFieldR2& X, float dt, int steps, const std::function<float(float, float, float, vec2, vec2)>& width,
					const std::function<vec4(float, float, float, vec2, vec2)>& color);
	void generateGrid(vec2 v_min, vec2 v_max, ivec2 res);
	void generateRandomUniform(vec2 v_min, vec2 v_max, int n);
	void generateStartTimesAll0();
	void generateStartTimesUniform(float t_max);
	void generateLine(int i);
	void generateLines();

	static float getTimeRelative(const BufferedVertex& v);
	static float getT0(const BufferedVertex& v);
	static float getTimeAbsolute(const BufferedVertex& v);
	static vec2 getPos(const BufferedVertex& v);
	static vec2 getStartPoint(const BufferedVertex& v);
	static float getSpeed(const BufferedVertex& v);
	static float getLength(const BufferedVertex& v);
	static vec4 getColor(const BufferedVertex& v);
	static float getWidth(const BufferedVertex& v);
};


class PlanarDiffusedInterval : public IndexedMesh {
	VectorFieldR2 X;
	float dt;
	int steps;
	vec4 color;
	HOM(float, float) time_dump;
	HOM(float, float) width_ratio_dump;
	vec2 a, b;
	float t0;
	PolyGroupID id = randomID();

	// uv = (t, t0)
	// pos = (x, y, 0)
	// n = (0, 0, 1)
	// col = (color)

public:
	PlanarDiffusedInterval(const VectorFieldR2& X, float dt, int steps, vec4 color, const HOM(float, float)& time_dump, const HOM(float, float)& width_ratio_dump, vec2 a, vec2 b,
						   float t0);
};


class PlanarDiffusedCurve : public IndexedMesh {
	VectorFieldR2 X;
	float dt;
	int steps;
	vec4 color;
	HOM(float, float) time_dump;
	HOM(float, float) width_ratio_dump;
	SmoothParametricPlaneCurve curve;
	vec2 domain;
	int resolution;
	float t0;

	// uv = (t, t0)
	// pos = (x, y, 0)
	// n = (0, 0, 1)
	// col = (color)

public:
	PlanarDiffusedCurve(const VectorFieldR2& X, float dt, int steps, vec4 color, const HOM(float, float)& time_dump, const HOM(float, float)& width_ratio_dump,
						const SmoothParametricPlaneCurve& curve, vec2 domain, int resolution, float t0);
};


class PlanarDiffusedPatterns : public IndexedMesh {
	// uv = (t, t0)
	// pos = (x, y, 0)
	// n = (0, 0, 1)
	// col = (color)

public:
	PlanarDiffusedPatterns(const VectorFieldR2& X, float dt, int steps, const vector<vec4>& colors, const vector<vec2>& shifts, const HOM(float, float)& time_dump,
						   const HOM(float, float)& width_ratio_dump, const SmoothParametricPlaneCurve& curve_pattern, vec2 domain, int resolution, float t0);
};


struct PIPE_SETTINGS {
	float radius;
	int horRes;
	int radialRes;
	dict(string, vec4) extra_defaults = {{"EXTRA0", vec4(0, 0, 0, 0)}};
	bool bounding = false;
	vec3 bound_min = vec3(-1, -1, -1);
	vec3 bound_max = vec3(1, 1, 1);
	vector<float> discontinuities = {};
};


/**
 *\class PipeCurveVertexShader
 *\brief curve with GPU pipe surface generating from flat mesh
 *
 *
 * \POS  curve point
 * \NORMAL: curve normal
 * \UV: (curve parameter \t, polar angle \theta)
 * \COLOR: (curve \binormal | pipe \radius)
 *
 *extra0 loaded by default, filled with 0s
 */

class PipeCurveVertexShader : public IndexedMesh {
	PIPE_SETTINGS settings;
	PolyGroupID id;

public:
	PipeCurveVertexShader(const SmoothParametricCurve& curve, float r, int horRes, int radialRes, const PolyGroupID& id = randomID(), dict(string, vec4) extra_defaults = {});
	PipeCurveVertexShader(const RealFunction& plot, vec2 dom, float r, int horRes, int radialRes, const PolyGroupID& id = randomID());
	PipeCurveVertexShader(const DiscreteRealFunction& plot, float r, int radialRes, const PolyGroupID& id = randomID());
	PipeCurveVertexShader(const DiscreteRealFunctionNonUniform& plot, float r, int radialRes, const PolyGroupID& id = randomID());

	PipeCurveVertexShader(const SmoothParametricCurve& curve, const PIPE_SETTINGS& s, const PolyGroupID& id = randomID());
	PipeCurveVertexShader(const RealFunction& plot, vec2 dom, const PIPE_SETTINGS& s, const PolyGroupID& id = randomID());
	PipeCurveVertexShader(const DiscreteRealFunction& plot, const PIPE_SETTINGS& s, const PolyGroupID& id = randomID());
	PipeCurveVertexShader(const DiscreteRealFunctionNonUniform& plot, const PIPE_SETTINGS& s, const PolyGroupID& id = randomID());

	void duplicateCurve(const PolyGroupID& copy_id);

	void updateCurve(const DiscreteRealFunction& plot);
	void updateCurve(const DiscreteRealFunction& plot, const DiscreteRealFunction& df_precomputed);
	void updateCurve(const DiscreteRealFunctionNonUniform& plot);
	void updateRadius(const HOM(float, float)& r);
	void updateRadius(float r);

	void updateCurve(const SmoothParametricCurve& curve);

	static vec3 getBinormal(const BufferedVertex& v);
	static float getAngle(const BufferedVertex& v);
	static float getRadius(const BufferedVertex& v);
	static float getParameter(const BufferedVertex& v);
	static vec4 getExtra(const BufferedVertex& v);
	static void setBinormal(BufferedVertex& v, vec3 b);
	static void setAngle(BufferedVertex& v, float theta);
	static void setRadius(BufferedVertex& v, float r);
	static void setParameter(BufferedVertex& v, float t);
	static void setExtra(BufferedVertex& v, float value, int extra_index);

	void transform(SpaceAutomorphism F);
};


class SurfacePlotDiscretisedMesh : public IndexedMesh {
public:
	explicit SurfacePlotDiscretisedMesh(const DiscreteRealFunctionR2& plot);
	void transform(SpaceAutomorphism F);
};


class SurfacePolarPlotDiscretisedMesh : public IndexedMesh {
public:
	explicit SurfacePolarPlotDiscretisedMesh(const DiscreteRealFunctionR2& plot, float r = 1, float rot_speed = 0);
};


class FoliatedParametricSurfaceMesh : public IndexedMesh {
	ParametricSurfaceFoliation foliation;
	int special_leaves, continuous_leaves, leaf_radial_res, leaf_hor_res;
	HOM(float, vec4) color_map, extra1_map;

	vector<vec4> special_leaf_colors, special_extra1s;

	HOM(float, float) leaf_radius_map;
	vector<float> special_leaf_radii;
	string id_prefix = randomStringLetters(10);

public:
	FoliatedParametricSurfaceMesh(ParametricSurfaceFoliation foliation, int special_leaves, int continuous_leaves, int leaf_radial_res, int leaf_hor_res,
								  HOM(float, vec4) color_map, const vector<vec4>& special_leaf_colors, HOM(float, float) leaf_radius_map, const vector<float>& special_leaf_radii,
								  HOM(float, vec4) extra1_map, const vector<vec4>& special_extra1s);

	FoliatedParametricSurfaceMesh(const FoliatedParametricSurfaceMesh& other);
	FoliatedParametricSurfaceMesh(FoliatedParametricSurfaceMesh&& other) noexcept;
	FoliatedParametricSurfaceMesh& operator=(const FoliatedParametricSurfaceMesh& other);
	FoliatedParametricSurfaceMesh& operator=(FoliatedParametricSurfaceMesh&& other) noexcept;
};


template <typename VertexStruct>
struct GeneralTriangle {
	VertexStruct v0, v1, v2;

	GeneralTriangle(VertexStruct v0, VertexStruct v1, VertexStruct v2)
	: v0(v0), v1(v1), v2(v2) {}

	VertexStruct& operator[](int i) {
		if (i == 0)
			return v0;
		if (i == 1)
			return v1;
		if (i == 2)
			return v2;
		throw std::out_of_range("Index out of range in GeneralTriangle");
	}

	VertexStruct operator[](int i) const {
		if (i == 0)
			return v0;
		if (i == 1)
			return v1;
		if (i == 2)
			return v2;
		throw std::out_of_range("Index out of range in GeneralTriangle");
	}
};

template <typename VertexStruct>
struct GeneralIndexedMeshData {
	vector<VertexStruct> vertices;
	vector<ivec3> indices;

	GeneralIndexedMeshData()
	: vertices() {}

	GeneralIndexedMeshData(GeneralIndexedMeshData& other) = delete;

	GeneralIndexedMeshData(const vector<VertexStruct>& vertices, const vector<ivec3>& indices)
	: vertices(vertices), indices(indices) {}
};


template <typename VertexStruct>
class GeneralIndexedMesh {
	unique_ptr<GeneralIndexedMeshData<VertexStruct>> data;
};
