#ifndef GEOMETRY_HPP
#define GEOMETRY_HPP

#include <vector>
#include <array>
#include <variant>
#include <map>

#include <glm/glm.hpp>
#include <string>

#include "func.hpp"
#include <gl/glew.h>

#define End1P std::function<SpaceEndomorphism(float)>
#define End2P std::function<SpaceEndomorphism(float, float)>

template<typename vec> 
std::vector<float> vecToVecHeHe(vec v);


enum GLSLType {
	FLOAT,
	INT,
	VEC2,
	VEC3,
	VEC4,
	MAT2,
	MAT3,
	MAT4,
	SAMPLER1D,
	SAMPLER2D,
	SAMPLER3D
};

class Texture {
public:
	int width, height;
	unsigned char* data;
	int size;
	GLuint textureID;
	GLenum textureSlot;
	const char* samplerName;
	GLuint frameBufferID;

	Texture(int width, int height, int slot = 0, const char* sampler = "tex");
	explicit Texture(const char* filename, int slot = 0, const char* sampler = "tex");
	~Texture();

	void addFilters(GLenum minFilter, GLenum magFilter, GLenum wrapS, GLenum wrapT);
	void bind();
	void bindToFrameBuffer();
	void calculateMipmap();
};


class MaterialPhong {
public:
	Texture* texture;
	glm::vec4 ambientColor;
	glm::vec4 diffuseColor;
	glm::vec4 specularColor;
	float ambientIntensity;
	float diffuseIntensity;
	float specularIntensity;
	float shininess;

	MaterialPhong();
	MaterialPhong(glm::vec4 ambient, glm::vec4 diffuse, glm::vec4 specular,
		float ambientIntensity, float diffuseIntensity, float specularIntensity,
		float shininess, Texture* texture=nullptr);

	bool textured();
	glm::mat4 compressToMatrix();
};


class TriangleR3 {
public:
	std::vector<glm::vec3> vertices;
	std::vector<glm::vec2> uvs;
	std::vector<glm::vec3> normals;
	std::vector<glm::vec4> vertexColors;
	std::map<std::string, std::array<glm::vec4, 3>> extraData = {};
	std::optional<MaterialPhong> material = std::nullopt;


	TriangleR3(std::vector<glm::vec3> vertices, std::vector<glm::vec3> normals, std::vector<glm::vec4> colors, std::vector<glm::vec2> uvs, MaterialPhong material);
	TriangleR3(std::vector<glm::vec3> vertices, std::vector<glm::vec3> normals, std::vector<glm::vec4> colors, std::vector<glm::vec2> uvs);
	TriangleR3(std::vector<glm::vec3> vertices, std::vector<glm::vec3> normals, std::vector<glm::vec4> colors);
	TriangleR3(std::vector<glm::vec3> vertices, std::vector<glm::vec4> colors);
	TriangleR3(std::vector<glm::vec3> vertices, std::vector<glm::vec3> normal);
	TriangleR3(std::vector<glm::vec3> vertices, std::vector<glm::vec3> normal, MaterialPhong material);
	explicit TriangleR3(std::vector<glm::vec3> vertices);
	TriangleR3(std::vector<glm::vec3> vertices, MaterialPhong material);
	TriangleR3(const std::vector<glm::vec3> &vertices, const std::vector<glm::vec2> &uvs);
	TriangleR3(std::vector<glm::vec3> vertices, glm::vec3 normal, std::vector<glm::vec4> colors, std::vector<glm::vec2> uvs);
	TriangleR3(std::vector<glm::vec3> vertices, glm::vec3 normal, std::vector<glm::vec4> colors);
	TriangleR3(std::vector<glm::vec3> vertices, glm::vec3 normal);
	TriangleR3(std::vector<glm::vec3> vertices, std::vector<glm::vec3> normals, glm::vec4 color, std::vector<glm::vec2> uvs);
	TriangleR3(std::vector<glm::vec3> vertices, std::vector<glm::vec3> normals, glm::vec4 color);
	TriangleR3(const std::vector<glm::vec3> &vertices, glm::vec4 color);
	TriangleR3(std::vector<glm::vec3> vertices, glm::vec3 normal, glm::vec4 color, std::vector<glm::vec2> uvs);
	TriangleR3(std::vector<glm::vec3> vertices, glm::vec3 normal, glm::vec4 color);
	TriangleR3();

	glm::vec3 operator[](int i) const;
	glm::vec3 barycenter() const;
	float area() const;
	glm::vec3 barycentricToWorld(glm::vec3 coords) const;
	glm::vec3 worldToBarycentric(glm::vec3 point) const;
	TriangleR3 operator+(glm::vec3 v) const;
	TriangleR3 operator*(const glm::mat4 &M) const;
	bool operator<(const TriangleR3& t) const;
	bool hasMaterial() const;

	void recalculateNormal();
	void flipNormals();

	void addExtraData(std::string name, std::array<glm::vec4, 3> data);
	std::array<glm::vec4, 3> getExtraData(std::string name);
	void addMaterial(MaterialPhong material);

};


class TriangleR2 {
public:
	std::vector<glm::vec2> vertices;
	std::vector<glm::vec2> uvs;
	std::vector<glm::vec4> vertexColors;

	TriangleR2(std::vector<glm::vec2> vertices, std::vector<glm::vec4> colors, std::vector<glm::vec2> uvs);
	TriangleR2(std::vector<glm::vec2> vertices, std::vector<glm::vec4> colors);
	explicit TriangleR2(std::vector<glm::vec2> vertices);
	TriangleR2(std::vector<glm::vec2> vertices, std::vector<glm::vec2> uvs);
	TriangleR2(glm::vec2 v1, glm::vec2 v2, glm::vec2 v3);
	TriangleR2();

	glm::vec2 operator[](int i);
	TriangleR2 operator+(glm::vec2 v);
	TriangleR3 embeddInR3(float z = 0);
	float area() const;
};

class TriangleComplex {
public:
	std::array<Complex, 3> vertices;
	std::array<glm::vec2, 3> uvs;
	std::array<glm::vec4, 3> vertexColors;

	TriangleComplex(std::array<Complex, 3> vertices, std::array<glm::vec4, 3> colors, std::array<glm::vec2, 3> uvs);
	TriangleComplex(std::array<Complex, 3> vertices, std::array<glm::vec2, 3> uvs);
	TriangleComplex(std::array<Complex, 3> vertices);

	operator TriangleR2() const;
	TriangleComplex operator+(Complex v) const;
	TriangleComplex operator*(Complex M) const;
	TriangleComplex operator*(Matrix<Complex, 2> M) const;
	TriangleComplex operator*(Meromorphism f) const;
	TriangleR3 embeddInR3(float z = 0) const;

	std::array<TriangleComplex, 3> subdivide() const;
	static std::vector<TriangleComplex> subdivideTriangulation(const std::vector<TriangleComplex> &triangles);
};


class Vertex {
public:
	glm::vec3 position;
	glm::vec3 normal;
	glm::vec2 uv;
	glm::vec4 color;
	Vertex(glm::vec3 position, glm::vec3 normal, glm::vec2 uv, glm::vec4 color);
	Vertex(glm::vec3 position, glm::vec3 normal, glm::vec2 uv);
	Vertex(glm::vec3 position, glm::vec3 normal, glm::vec4 color);
	Vertex(glm::vec3 position, glm::vec3 normal);
	Vertex(glm::vec3 position, glm::vec2 uv, glm::vec4 color);
	Vertex(glm::vec3 position, glm::vec2 uv);
	Vertex(glm::vec3 position, glm::vec4 color);
	explicit Vertex(glm::vec3 position);
	Vertex();

	std::string hash() const;

	bool operator<(const Vertex& v) const;
	bool operator==(const Vertex& v) const;
	void translate(glm::vec3 v);
	void transform(const glm::mat4 &M);
	Vertex operator+(glm::vec3 v) const;
	Vertex operator*(glm::mat4 M) const;
	Vertex operator*(glm::mat3 M) const;
	void operator+=(glm::vec3 v) { this->translate(v); }
	void operator*=(glm::mat4 M) { this->transform(M); }
	void operator*=(glm::mat3 M) { this->transform(glm::mat4(M)); }
};



MaterialPhong lerp(MaterialPhong m0, MaterialPhong m1, float t);



class MaterialFamily1P {
	MaterialPhong m0, m1;
public:
	MaterialFamily1P(MaterialPhong m0, MaterialPhong m1);
	MaterialFamily1P(glm::vec4 c1, glm::vec4 c2, float ambientIntensity, float diffuseIntensity, float specularIntensity, float shininess);
	MaterialPhong operator()(float t) const;
};


class PointLight {
public:
	glm::vec3 position;
	glm::vec4 color;
	float intensity;
	PointLight(glm::vec3 position, glm::vec4 color, float intensity);
	PointLight();
	glm::mat4 compressToMatrix();
};

enum MeshFormat {
	OBJ = 0
};


class TriangularMesh {
protected:
	std::vector <TriangleR3> triangles;

	std::vector <glm::vec3> posBuff;
	std::vector <glm::vec3> normBuff;
	std::vector <glm::vec4> colorBuff;
	std::vector <glm::vec2> uvBuff;

	std::vector <glm::vec3> calculatePositionBuffer() const;
	std::vector <glm::vec3> calculateNormalBuffer() const;
	std::vector <glm::vec4> calculateColorBuffer() const;
	std::vector <glm::vec2> calculateUVBuffer() const;

	std::map<std::string, std::vector<glm::vec4>> extraBuff = {};

public:
    std::array<void*, 4> bufferLocations;
    std::array<int, 4> bufferSizes;

	TriangularMesh();
	explicit TriangularMesh(const std::vector <TriangleR3> &triangles);
	TriangularMesh(std::vector <Vertex> vertices, std::vector <glm::ivec3> faceIndices);
	explicit TriangularMesh(const char* filename, MeshFormat format = OBJ);
	TriangularMesh(ComplexCurve* curve, int nSegments, float h_middle, float w_middle, float w_side);
	std::vector<TriangleR3> getTriangles() const;

	// --- transforms ---
	void translate(glm::vec3 v);
	void transform(const glm::mat4 &M);
	TriangularMesh operator+(glm::vec3 v) const;
	TriangularMesh operator*(const glm::mat4 &M) const;
	TriangularMesh operator*(const glm::mat3 &M) const;
	void operator+=(glm::vec3 v) { this->translate(v); }
	void operator*=(glm::mat4 M) { this->transform(M); }
	void operator*=(glm::mat3 M) { this->transform(glm::mat4(M)); }
    void applyMap(const std::function<glm::vec3(glm::vec3)> &f);
    void applyMap(const std::function<glm::vec3(glm::vec3)> &f,
             const std::function<glm::vec3(glm::vec3, glm::vec3)> &f_normal);
	void applyMap(const std::function<glm::vec3(glm::vec3)> &f, std::string customDomain);
	void applyMap(std::function<glm::vec3(glm::vec3)> f, std::function<glm::vec3(glm::vec3, glm::vec3)> f_normal, std::string customDomain);

	// --- buffer generators ---

	void precomputeBuffers();
	void precomputeExtraBuffer(std::string name);
	void recomputeBuffers();
	void recomputeBuffers(bool pos, bool norm, bool color, bool uv);

	// --- modifiers ---
	void randomizeFaceColors(glm::vec3 min = glm::vec3(0, 0, 0), glm::vec3 max = glm::vec3(1, 1, 1));
	void cleanUpEmptyTriangles();
	void flipNormals();
	void recalculateNormals();
};


class IndexedMesh {
public:
	std::vector <Vertex> vertices;
	std::vector<std::array<int, 3>> faceIndices;
	std::vector <glm::vec3> posBuff;
	std::vector <glm::vec3> normBuff;
	std::vector <glm::vec4> colorBuff;
	std::vector <glm::vec2> uvBuff;

	IndexedMesh(const std::vector <Vertex> &vertices, const std::vector<std::array<int, 3>> &faceIndices);
	explicit IndexedMesh(const char* filename, MeshFormat format = OBJ);
	explicit operator TriangularMesh();

	std::vector <glm::vec3> calculatePositionBuffer();
	std::vector <glm::vec3> calculateNormalBuffer();
	std::vector <glm::vec4> calculateColorBuffer();
	std::vector <glm::vec2> calculateUVBuffer();
	void precomputeBuffers();

	void translate(glm::vec3 v);
	void transform(glm::mat4 M);
	IndexedMesh operator+(glm::vec3 v) const;
	IndexedMesh operator*(glm::mat4 M) const;
	IndexedMesh operator*(glm::mat3 M) const;
	void operator+=(glm::vec3 v) { this->translate(v); }
	void operator*=(glm::mat4 M) { this->transform(M); }
	void operator*=(glm::mat3 M) { this->transform(glm::mat4(M)); }

	void randomiseVertexColors();
	void flipNormals();
};


class Model3D {
public:
	std::shared_ptr<TriangularMesh> mesh;
	std::shared_ptr<MaterialPhong> material;
	glm::mat4 transform;
	Model3D();
	Model3D(TriangularMesh &mesh, MaterialPhong &material, glm::mat4 transform);
	void addTransform(glm::mat4 transform);
};

enum BoundaryEmbeddingType {
	KERB,
	PIPE,
	FLAT,
	CURVE
};
struct BoundaryEmbeddingStyle {
	BoundaryEmbeddingType type;
	float width = 0;
	float height = 0;
	float outerMargin = 0;
	float skewness = 0;
	float width_middle = 0;
	float height_middle = 0;
	float width_side = 0;
	int nSegments = 60;
};

class PlanarMeshWithBoundary {
public:
	std::vector<TriangleR2> triangles;
	std::vector < std::vector<glm::vec2>> boundaries;
	std::vector<bool> boundaryCyclic;
	bool encodesVectorField = false;
	

	PlanarMeshWithBoundary();
	PlanarMeshWithBoundary(std::vector<TriangleR2> triangles, std::vector < std::vector<glm::vec2>> boundaries, std::vector<bool> boundaryCyclic);
	PlanarMeshWithBoundary(std::vector<TriangleR2> triangles, std::vector < std::vector<glm::vec2>> boundaries);
	PlanarMeshWithBoundary(std::vector<TriangleR2> triangles, std::vector<glm::vec2> boundary, bool cyclic = true);
	explicit PlanarMeshWithBoundary(std::vector<TriangleR2> triangles);
	
	void addVectorField(VectorFieldR2 vectorField);

	TriangularMesh embeddInR3(float z = 0) const;
	TriangularMesh embeddInR3LowerByArea(float z=0, float factor=2) const;

	PlanarMeshWithBoundary offsetBoundaryMesh(float w);
	std::vector<TriangularMesh> kerbBoundaryEmbedding(float w, float h, float outerMargin, float skewness);
	std::vector<TriangularMesh> stylisedBoundaryEmbedding(BoundaryEmbeddingStyle style); // only kerb for now implemented
};

class MeshFamily1P{
protected:
	float _t=0;
	std::shared_ptr<TriangularMesh> mesh;
	std::function<std::function<glm::vec3(glm::vec3)>(float, float)> time_operator;
	std::function<std::function<glm::vec3(glm::vec3, glm::vec3)>(float, float)> time_operator_normal;
	bool planar=false;
public:
	MeshFamily1P(std::shared_ptr<TriangularMesh> mesh, std::function<std::function<glm::vec3(glm::vec3)>(float, float)> time_operator, std::function<std::function<glm::vec3(glm::vec3, glm::vec3)>(float, float)> time_operator_normal, float t=0);
	MeshFamily1P(std::shared_ptr<TriangularMesh> embedded_mesh, std::function<PlaneSmoothEndomorphism(float, float)> time_operator, glm::vec3 embedding_shift=glm::vec3(0, 0, 0), float t=0);
	MeshFamily1P(std::shared_ptr<TriangularMesh> embedded_mesh, std::function<PlaneAutomorphism(float)> time_operator_from0, glm::vec3 embedding_shift=glm::vec3(0, 0, 0));

	float time() const;
	virtual void transformMesh(float new_t);
	virtual void meshDeformation(float dt);
};


class MeshFamily1PExtraDomain : public MeshFamily1P {
protected:
	std::string domain;
public:
	MeshFamily1PExtraDomain(std::shared_ptr<TriangularMesh> mesh, std::string domain, std::function<std::function<glm::vec3(glm::vec3)>(float, float)> time_operator, std::function<std::function<glm::vec3(glm::vec3, glm::vec3)>(float, float)> time_operator_normal, float t=0);
	MeshFamily1PExtraDomain(std::shared_ptr<TriangularMesh> embedded_mesh, std::string domain, std::function<PlaneSmoothEndomorphism(float, float)> time_operator, glm::vec3 embedding_shift=glm::vec3(0, 0, 0), float t=0);
	MeshFamily1PExtraDomain(std::shared_ptr<TriangularMesh> embedded_mesh, std::string domain, std::function<PlaneAutomorphism(float)> time_operator_from0, glm::vec3 embedding_shift=glm::vec3(0, 0, 0));

	void transformMesh(float new_t) override;
	void meshDeformation(float dt) override;
};



#define PolyGroupID std::variant<int, std::string>

inline std::string polyGroupIDtoString(PolyGroupID id) {
    if (std::holds_alternative<int>(id)) {
        return std::to_string(std::get<int>(id));
    } else {
        return std::get<std::string>(id);
    }
}

inline PolyGroupID prefix(PolyGroupID id, std::string prefix) {
    return PolyGroupID(prefix + polyGroupIDtoString(id));
}

inline std::string randomString() {
    int a = rand();
	return std::to_string(a);
}

inline PolyGroupID make_unique_id(PolyGroupID id) {
    return prefix(id, randomString());
}




inline PolyGroupID bdGroup(int n) { return PolyGroupID("bd" + std::to_string(n));}

inline PolyGroupID curveGroup(int n) { return PolyGroupID("curva" + std::to_string(n));}

struct MaterialBuffers {
    std::vector<glm::vec4> ambientColors;
	std::vector<glm::vec4> diffuseColors;
	std::vector<glm::vec4> specularColors;
	std::vector<glm::vec4> intencitiesAndShininess;
};

struct StdAttributeBuffers {
	std::vector<glm::vec3>  positions;
	std::vector<glm::vec3> normals;
	std::vector<glm::vec4> colors;
	std::vector<glm::vec2>  uvs;
};



const BoundaryEmbeddingStyle STD_KERB = BoundaryEmbeddingStyle{KERB, .01f, .02f, .01f, .2};
const BoundaryEmbeddingStyle STD_CURVE = BoundaryEmbeddingStyle{CURVE, .0f, .0f, .0f, .0f, .005f, .005f, .005f};


class SuperMesh {
	std::map<PolyGroupID, std::vector<TriangleR3>> triangleGroups={};
	std::map<PolyGroupID, std::vector<TriangleR3>> boundaryGroups={}; // "bdPoint": extra attribute for unextruded pt
	std::map<PolyGroupID, std::vector<TriangleR3>> embedded_curves={}; // "curvePoint": extra attribute for unextruded pt
	std::map<PolyGroupID, std::vector<TriangleR3>> embedded_points={}; // "loc": extra const attribute for unextruded pt
	std::map<PolyGroupID, MaterialPhong> materials={};

	MaterialBuffers materialBuffers = {};
	StdAttributeBuffers stdAttributeBuffers = {};
	std::map<std::string, std::vector<glm::vec4>> extraBuff = {};
	std::map<std::string, int> extraBufferIndices = {};

public:
	std::vector<void*> bufferLocations;
	std::vector<int> bufferSizes;
	SuperMesh();
	explicit SuperMesh(TriangularMesh &mesh);
	SuperMesh(const char* filename, MeshFormat format = OBJ);
	SuperMesh(TriangularMesh &mesh, const MaterialPhong &material);
	SuperMesh(const char* filename, MaterialPhong &material, MeshFormat format = OBJ);
	SuperMesh(PlanarMeshWithBoundary &mesh, MaterialPhong &material, MaterialPhong &material_bd, BoundaryEmbeddingStyle embeddingStyle = STD_KERB);
	SuperMesh(const std::vector <TriangleR3> &triangles, const MaterialPhong &material);
	SuperMesh(const std::map<PolyGroupID, std::vector<TriangleR3>> &triangleGroups, const std::map<PolyGroupID, MaterialPhong> &materials);

	void addPolyGroup(PolyGroupID id, const std::vector<TriangleR3> &triangles, const MaterialPhong &material);
	void addPolyGroup(const std::vector<TriangleR3> &triangles, const MaterialPhong &material);
	void addPolyGroup(const TriangularMesh &mesh, const MaterialPhong &material);
	void addBdGroup(PolyGroupID id, const std::vector<TriangleR3> &triangles, const MaterialPhong &material);
	void addBdGroup(const std::vector<TriangleR3> &triangles, const MaterialPhong &material);
	void addEmbeddedCurve(PolyGroupID id, const std::vector<TriangleR3> &triangles, const MaterialPhong &material);
	void addEmbeddedCurve(const std::vector<TriangleR3> &triangles, const MaterialPhong &material);
	void addEmbeddedPoint(PolyGroupID id, const std::vector<TriangleR3> &triangles, const MaterialPhong &material);
	void addEmbeddedPoint(const std::vector<TriangleR3> &triangles, const MaterialPhong &material);
	void embedCurve(ComplexCurve* curve, int nSegments, float h_middle, float w_middle, float w_side, const MaterialPhong &material);
	void precomputeBuffers(bool materials=true, bool extra=true);
	void merge(const SuperMesh &other);
	void precomputeExtraBuffer(std::string name);

	void actOnPositions(std::function<glm::vec3(glm::vec3)> f);
	void actAtEmbeddedPlane(Meromorphism f);
	void actOnPositionsWithCustomShift(std::function<glm::vec3(glm::vec3)> f, std::map<PolyGroupID, std::string> useShiftOfCustomDomain = {});
	void translate(glm::vec3 v);
	void doPerTriangle(std::function<void(TriangleR3&)> f);
	void doPerTriangle(PolyGroupID grp, std::function<void(TriangleR3&)> f);

	void actOnEmbeddedCurve(SpaceEndomorphism f);


	void randomizeMaterials(MaterialPhong &min, MaterialPhong &max);
	void randomiseMaterialsDynamically(float stepMax);
};





class SuperPencilPlanar : public SuperMesh{
	float _t=0;
	std::function<Biholomorphism(float)> _time_operator = [](float t){return Id;};

public:
	using SuperMesh::SuperMesh;
	void makePencil(std::function<Biholomorphism(float)> _time_operator, float t=0);
	void transformMesh(float new_t);
	void deformMesh(float dt);

};

class CurveSample {
	glm::vec3 position;
	glm::vec3 normal;
	glm::vec3 tangent;
	MaterialPhong material;
	float width;
	glm::vec4 extraInfo = glm::vec4(0, 0, 0, 0);

public:
	CurveSample(glm::vec3 position, glm::vec3 normal, glm::vec3 tangent, MaterialPhong material, float width);
	glm::vec3 getPosition() const {return position;}
	glm::vec3 getNormal() const {return normal;}
	glm::vec3 getTangent() const {return tangent;}
	MaterialPhong getMaterial() const {return material;}
	float getWidth() const {return width;}
	glm::vec3 getBinormal() const {return glm::cross(tangent, normal);}
	void addExtra(glm::vec4 extra) {extraInfo = extra;}
	void addExtra(glm::vec3 extra) {extraInfo = glm::vec4(extra, 0);}
	void addExtra(float extra, int position=4) {extraInfo[position] = extra;}
	glm::vec4 readExtra() const {return extraInfo;}
	glm::vec3 readExtraXYZ() const {return glm::vec3(extraInfo);}
	float readExtra(int position) const {return extraInfo[position];}
	float readExtraLast() const	{return extraInfo.w;}
	void updatePosition(glm::vec3 new_position) {position = new_position;}
	void updateNormal(glm::vec3 new_normal) {normal = new_normal;}
	void updateTangent(glm::vec3 new_tangent) {tangent = new_tangent;}
	void updateMaterial(MaterialPhong new_material) {material = new_material;}
	void updateWidth(float new_width) {width = new_width;}
};

std::vector<CurveSample> sampleCurve(SmoothParametricCurve curve, std::function<float(float)> width,
	std::function<MaterialPhong(float)> material, float t0=0, float t1=1, int n=100, bool periodic=false);



std::vector<CurveSample> sampleCurve(SmoothParametricCurve curve, float width,
                                     MaterialPhong material, float t0=0, float t1=1, int n=100, bool periodic=false);


class SuperCurve{
	float t0, t1;
	std::vector<CurveSample> samples;
	std::shared_ptr<SuperMesh> _mesh = nullptr;
public:
	SuperCurve(const SmoothParametricCurve &curve, const std::function<float(float)> &width, const std::function<MaterialPhong(float)> &material, float t0=0, float t1=1, int segments=100, bool periodic=true);
	SuperCurve(SmoothParametricCurve curve, float width, MaterialPhong material, float t0=0, float t1=1, int nSegments=100, bool periodic=true);
	SuperMesh mesh(int radialSegments) const;
	std::shared_ptr<SuperMesh> associateMesh(int radialSegments);

	void transformMeshByAmbientMap(
            const SpaceEndomorphism &f); // just composition
	void internalDeformation(End1P);
	void updateCurveMeshOnly(const SmoothParametricCurve &new_curve);
	void updateCurve(const SmoothParametricCurve &new_curve);


};

class SuperPencilCurve : public SuperCurve {
	float _t = 0;
	std::unique_ptr<End2P> _ambient_operator = nullptr;
	std::unique_ptr<std::function<SmoothParametricCurve(float)>> _parametric_operator = nullptr;
public:
	using SuperCurve::SuperCurve;
	void addAmbientDeformation(End2P _ambient_operator, float t=0);
	void addLocalDeformation(End1P _local_operator, float t=0);
	void addDeformationAlongVectorField(VectorFieldR3 vectorField, float t=0);
	void addPencil(std::function<SmoothParametricCurve(float)> family,
                       float t = 0);
        void updateCurveMeshOnly(SmoothParametricCurve &new_curve);
	float time() const;
	void transformMesh(float new_t);
};



#endif

