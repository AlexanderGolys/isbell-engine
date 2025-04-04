#pragma once
#include "../geometry/smoothImplicit.hpp"

#include <array>
#include <map>
#include <variant>
#include <vector>

#include <string>
#include <glm/glm.hpp>

#include <memory>
#include <optional>


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
    bool alpha=false;

	Texture(int width, int height, int slot = 0, const char* sampler = "tex");
    explicit Texture(vec3 color, int slot = 0, const char* sampler = "tex");
    explicit Texture(vec4 color, int slot = 0, const char* sampler = "tex");
	explicit Texture(const char* filename, int slot, const char* sampler, bool alpha=false);
	explicit Texture(const std::string& filename, int slot, const char* sampler, bool alpha=false);

	void deleteTexture();
	~Texture() { deleteTexture(); }

	static void addFilters(GLenum minFilter, GLenum magFilter, GLenum wrapS, GLenum wrapT);
	void bind() const;
	void bindToFrameBuffer() const;
	void calculateMipmap() const;
    void load();
};


class MaterialPhong {
	float ambientIntensity;
	float diffuseIntensity;
	float specularIntensity;
	float shininess;


public:
	std::shared_ptr<Texture> texture_ambient;
	std::shared_ptr<Texture> texture_diffuse;
	std::shared_ptr<Texture> texture_specular;

	virtual ~MaterialPhong() = default;

	MaterialPhong() : MaterialPhong(vec4(0), vec4(0), vec4(0), 0, 0, 0, 0) {}

	static std::shared_ptr<Texture> constAmbientTexture(vec4 color);
	static std::shared_ptr<Texture> constDiffuseTexture(vec4 color);
	static std::shared_ptr<Texture> constSpecularTexture(vec4 color);
	static std::shared_ptr<Texture> constAmbientTexture(vec3 color) { return constAmbientTexture(vec4(color, 1)); }
	static std::shared_ptr<Texture> constDiffuseTexture(vec3 color) { return constDiffuseTexture(vec4(color, 1)); }
	static std::shared_ptr<Texture> constSpecularTexture(vec3 color) { return constSpecularTexture(vec4(color, 1)); }
	static std::shared_ptr<Texture> makeAmbientTexture(const std::string &filename, bool alpha=false) { return std::make_shared<Texture>(filename, 0, "texture_ambient", alpha); }
	static std::shared_ptr<Texture> makeDiffuseTexture(const std::string &filename, bool alpha=false) { return std::make_shared<Texture>(filename, 1, "texture_diffuse", alpha); }
	static std::shared_ptr<Texture> makeSpecularTexture(const std::string &filename, bool alpha=false) { return std::make_shared<Texture>(filename, 2, "texture_specular", alpha); }

  MaterialPhong(const std::shared_ptr<Texture> &texture_ambient, const std::shared_ptr<Texture> &texture_diffuse, const std::shared_ptr<Texture> &texture_specular,
                float ambientIntensity, float diffuseIntensity, float specularIntensity, float shininess);



	MaterialPhong(const std::string& textureFilenameAmbient, const std::string& textureFilenameDiffuse, const std::string& textureFilenameSpecular,
		float ambientIntensity, float diffuseIntensity, float specularIntensity, float shininess,
		bool alphaAmbient=false, bool alphaDiffuse=false, bool alphaSpecular=false)
						:  MaterialPhong(makeAmbientTexture(textureFilenameAmbient, alphaAmbient),
								makeDiffuseTexture(textureFilenameDiffuse, alphaDiffuse),
								makeSpecularTexture(textureFilenameSpecular, alphaSpecular),
										ambientIntensity, diffuseIntensity, specularIntensity, shininess) {}



	MaterialPhong(vec4 ambient, vec4 diffuse, vec4 specular, float ambientIntensity, float diffuseIntensity, float specularIntensity, float shininess)
	: MaterialPhong(constAmbientTexture(ambient), constDiffuseTexture(diffuse), constSpecularTexture(specular), ambientIntensity, diffuseIntensity, specularIntensity, shininess) {}

	MaterialPhong(vec4 ambient, float ambientIntensity, float diffuseIntensity, float specularIntensity, float shininess)
		: MaterialPhong(ambient, ambient, ambient, ambientIntensity, diffuseIntensity, specularIntensity, shininess) {}

	vec4 compressIntencities() const;
    void initTextures() { texture_ambient->load(); texture_diffuse->load(); texture_specular->load(); }

	virtual bool textured() const {return true;}

	void setAmbientTexture(const std::shared_ptr<Texture> &texture) {
	    texture_ambient->deleteTexture();
    	texture_ambient = texture;
    	texture_ambient->load();
    }
	void setDiffuseTexture(const std::shared_ptr<Texture> &texture) {
	    texture_diffuse->deleteTexture();
    	texture_diffuse = texture;
    	texture_diffuse->load();
    }
	void setSpecularTexture(const std::shared_ptr<Texture> &texture) {
	    texture_specular->deleteTexture();
    	texture_specular = texture;
    	texture_specular->load();
    }

	void addAmbientTexture (const std::string &filename, bool alpha=false) { setAmbientTexture(makeAmbientTexture(filename, alpha));}
	void addDiffuseTexture (const std::string &filename, bool alpha=false) { setDiffuseTexture(makeDiffuseTexture(filename, alpha));}
	void addSpecularTexture(const std::string &filename, bool alpha=false) { setSpecularTexture(makeSpecularTexture(filename, alpha));}

	void addAmbientTexture (vec4 color) { setAmbientTexture(constAmbientTexture(color)); }
	void addDiffuseTexture (vec4 color) { setDiffuseTexture(constDiffuseTexture(color)); }
	void addSpecularTexture(vec4 color) { setSpecularTexture(constSpecularTexture(color)); }
};

class MaterialPhongConstColor : public MaterialPhong {
public:
	vec4 ambientColor, diffuseColor, specularColor;
	MaterialPhongConstColor(vec4 ambient, vec4 diffuse, vec4 specular,
		float ambientIntensity, float diffuseIntensity, float specularIntensity, float shininess)
		: MaterialPhong(ambient, diffuse, specular, ambientIntensity, diffuseIntensity, specularIntensity, shininess),
		ambientColor(ambient), diffuseColor(diffuse), specularColor(specular) {}

	explicit MaterialPhongConstColor(mat4 mat) : MaterialPhongConstColor(mat[0], mat[1], mat[2], mat[3][0], mat[3][1], mat[3][2], mat[3][3]) {}
	MaterialPhongConstColor() : MaterialPhongConstColor(vec4(0), vec4(0), vec4(0), 0, 0, 0, 0) {}

	mat4 compressToMatrix() const {return mat4(ambientColor, diffuseColor, specularColor, compressIntencities());}
	bool textured() const override {return false;}
};


class TriangleR3 {
public:
	std::vector<vec3> vertices;
	std::vector<vec2> uvs;
	std::vector<vec3> normals;
	std::vector<vec4> vertexColors;
	std::map<std::string, std::array<vec4, 3>> extraData = {};
	std::optional<MaterialPhong> material = std::nullopt;


	TriangleR3(std::vector<vec3> vertices, std::vector<vec3> normals, std::vector<vec4> colors, std::vector<vec2> uvs, MaterialPhong material);
	TriangleR3(std::vector<vec3> vertices, std::vector<vec3> normals, std::vector<vec4> colors, std::vector<vec2> uvs);
	TriangleR3(std::vector<vec3> vertices, std::vector<vec3> normals, std::vector<vec4> colors);
	TriangleR3(std::vector<vec3> vertices, std::vector<vec4> colors);
	TriangleR3(std::vector<vec3> vertices, std::vector<vec3> normal);
	TriangleR3(std::vector<vec3> vertices, std::vector<vec3> normal, MaterialPhong material);
	explicit TriangleR3(std::vector<vec3> vertices);
	TriangleR3(std::vector<vec3> vertices, MaterialPhong material);
	TriangleR3(const std::vector<vec3> &vertices, const std::vector<vec2> &uvs);
	TriangleR3(std::vector<vec3> vertices, vec3 normal, std::vector<vec4> colors, std::vector<vec2> uvs);
	TriangleR3(std::vector<vec3> vertices, vec3 normal, std::vector<vec4> colors);
	TriangleR3(std::vector<vec3> vertices, vec3 normal);
	TriangleR3(std::vector<vec3> vertices, std::vector<vec3> normals, vec4 color, std::vector<vec2> uvs);
	TriangleR3(std::vector<vec3> vertices, std::vector<vec3> normals, vec4 color);
	TriangleR3(const std::vector<vec3> &vertices, vec4 color);
	TriangleR3(std::vector<vec3> vertices, vec3 normal, vec4 color, std::vector<vec2> uvs);
	TriangleR3(std::vector<vec3> vertices, vec3 normal, vec4 color);
	TriangleR3();

	vec3 operator[](int i) const;
	vec3 barycenter() const;
	float area() const;
	vec3 barycentricToWorld(vec3 coords) const;
	vec3 worldToBarycentric(vec3 point) const;
	TriangleR3 operator+(vec3 v) const;
	TriangleR3 operator*(const mat4 &M) const;
	bool operator<(const TriangleR3& t) const;
	bool hasMaterial() const;

	void recalculateNormal();
	void flipNormals();

	void addExtraData(std::string name, std::array<vec4, 3> data);
    void setExtraData(std::string name, int i, vec4 data) { extraData[name][i] = data; }
	std::array<vec4, 3> getExtraData(std::string name);
	void addMaterial(MaterialPhong material);

};


class TriangleR2 {
public:
	std::vector<vec2> vertices;
	std::vector<vec2> uvs;
	std::vector<vec4> vertexColors;

	TriangleR2(std::vector<vec2> vertices, std::vector<vec4> colors, std::vector<vec2> uvs);
	TriangleR2(std::vector<vec2> vertices, std::vector<vec4> colors);
	explicit TriangleR2(std::vector<vec2> vertices);
	TriangleR2(std::vector<vec2> vertices, std::vector<vec2> uvs);
	TriangleR2(vec2 v1, vec2 v2, vec2 v3);
	TriangleR2();

	vec2 operator[](int i);
	TriangleR2 operator+(vec2 v);
	TriangleR3 embeddInR3(float z = 0);
	float area() const;
};

class Meromorphism;

class TriangleComplex {
public:
	std::array<Complex, 3> vertices;
	std::array<vec2, 3> uvs;
	std::array<vec4, 3> vertexColors;

	TriangleComplex(std::array<Complex, 3> vertices, std::array<vec4, 3> colors, std::array<vec2, 3> uvs);
	TriangleComplex(std::array<Complex, 3> vertices, std::array<vec2, 3> uvs);
	TriangleComplex(std::array<Complex, 3> vertices);

	operator TriangleR2() const;
	TriangleComplex operator+(Complex v) const;
	TriangleComplex operator*(Complex M) const;
	TriangleComplex operator*(const Matrix<Complex> &M) const;
	TriangleComplex operator*(Meromorphism f) const;
	TriangleR3 embeddInR3(float z = 0) const;

	std::array<TriangleComplex, 3> subdivide() const;
	static std::vector<TriangleComplex> subdivideTriangulation(const std::vector<TriangleComplex> &triangles);
};

struct MaterialBuffers {
    std::vector<vec4> ambientColors;
	std::vector<vec4> diffuseColors;
	std::vector<vec4> specularColors;
	std::vector<vec4> intencitiesAndShininess;
};

struct StdAttributeBuffers {
	std::vector<vec3>  positions;
	std::vector<vec3> normals;
	std::vector<vec4> colors;
	std::vector<vec2>  uvs;
};

class IndexedTriangle;
class CurveSample;

class Vertex {
	vec3 position;
	vec3 normal;
	vec2 uv;
	vec4 color;
    maybeMaterial material;
    std::map<std::string, vec4> extraData = {};
    int index = -1;
    std::vector<std::pair<std::weak_ptr<IndexedTriangle>, int>> triangles = {};
public:
	Vertex(vec3 position,
	       vec2 uv,
	       vec3 normal=e3,
	       vec4 color=BLACK,
	       maybeMaterial material=std::nullopt,
	       std::map<std::string, vec4> extraData = {})
	{ this->position = position; this->normal = normal; this->uv = uv;
		this->color = color; this->material = material; this->extraData = extraData; }

	Vertex(vec3 position, vec3 normal=vec3(0,0,1),  vec4 color=vec4(0,0,0,1), // set projected position as uv
	       maybeMaterial material=std::nullopt) : Vertex(position, vec2(position.x, position.y), normal, color, material) {}

	std::string hash() const;
	bool operator<(const Vertex& v) const;
	bool operator==(const Vertex& v) const;
    void setIndex(int i) { index = i; }
    int getIndex() const { return index; }
    bool hasIndex() const { return index != -1; }
    void gotAddedAsVertex (std::weak_ptr<IndexedTriangle> triangle, int i) { triangles.push_back({triangle, i}); }
    std::vector<std::pair<std::weak_ptr<IndexedTriangle>, int>> getTriangles() const { return triangles; }
    std::vector<std::weak_ptr<Vertex>> getNeighbours() const;
    void recomputeNormals(bool weightByArea=true);

	vec3 getPosition() const { return position; }
	vec3 getNormal() const { return normal; }
	vec2 getUV() const { return uv; }
	vec4 getColor() const { return color; }
    MaterialPhong getMaterial() const { return material.value(); }
    mat4 getMaterialMat() const { return material.value().compressToMatrix(); }
    vec4 getExtraData(const std::string &name) { return extraData[name]; }
    vec3 getExtraData_xyz(const std::string &name) { return vec3(extraData[name]); }
    float getExtraData(const std::string &name, int i) { return extraData[name][i]; }
    float getExtraLast(const std::string &name) { return extraData[name][3]; }
    bool hasExtraData() { return !extraData.empty(); }
    std::vector<std::string> getExtraDataNames();
    bool hasMaterial() { return material.has_value(); }

    void addExtraData(const std::string &name, vec4 data) { extraData[name] = data; }
    void addExtraData(const std::string &name, vec3 data);
    void addExtraData(const std::string &name, float data, int i=3);
	void translate(vec3 v) { this->position += v; }
	void transform(const SpaceEndomorphism &M);
	Vertex operator+(vec3 v) const;
	void operator+=(vec3 v) { this->translate(v); }
	void setMaterial(MaterialPhongConstColor material) { this->material = material; }
	void setPosition(vec3 position) { this->position = position; };
    void setNormal(vec3 normal) { this->normal = normal; }
	void setUV(vec2 uv) { this->uv = uv; }
	void setColor(vec4 color) { this->color = color; }
    void appendToBuffers(StdAttributeBuffers &buffers, MaterialBuffers &materialBuffers);
    void appendToBuffers(StdAttributeBuffers &buffers);
    void appendExtraDataToBuffer(const std::string &name, std::vector<vec4> &buffer) const;

    void appendToList(std::vector<Vertex> &list);
    void addTriangle(std::shared_ptr<IndexedTriangle> triangle, int i);

    void setCurveParameter(float t) { addExtraData("curvePoint", t); }
    void setCurvePosition(vec3 pos) { addExtraData("curvePoint", pos); }
    void setCurveTangent(vec3 v) { addExtraData("curveTangent", v); }
    void setCurveNormal(vec3 v) { addExtraData("curveNormal", v); }
    void setCurveNormalAngle(float a) { addExtraData("curveTangent", a); }
    void setAllParametricCurveExtras(float t, CurveSample &sample);

    float getCurveParameter() { return getExtraLast("curvePoint"); }
    vec3 getCurvePosition() { return getExtraData_xyz("curvePoint"); }
    vec3 getCurveTangent() { return getExtraData_xyz("curveTangent"); }
    vec3 getCurveNormal() { return getExtraData_xyz("curveNormal"); }
    float getCurveNormalAngle() { return getExtraLast("curveTangent"); }
    float getCurveWidth() { return norm(getPosition() - getCurvePosition()); }
};

Vertex barycenter(Vertex v1, Vertex v2, Vertex v3);
Vertex center(Vertex v1, Vertex v2);




MaterialPhong lerp(MaterialPhong m0, MaterialPhong m1, float t);



class MaterialFamily1P {
	MaterialPhong* m0 = nullptr;
    MaterialPhong *m1 = nullptr;
public:
	MaterialFamily1P(MaterialPhong *m0, MaterialPhong *m1);
	MaterialFamily1P(vec4 c1, vec4 c2, float ambientIntensity, float diffuseIntensity, float specularIntensity, float shininess);
	MaterialPhong operator()(float t) const;
};



class Light {
protected:
	vec3 position;
	vec4 color;
	vec4 intensities;
	float mode;
	float softShadowRadius;
public:
	Light(vec3 position, vec4 color, vec4 intensities, float mode, float r) : position(position), color(color), intensities(intensities), mode(mode), softShadowRadius(r) {}
	virtual ~Light() = default;
	virtual mat4 compressToMatrix() const { return mat4(vec4(position, 0), color, intensities, vec4(mode, softShadowRadius, 0, 0)); }
	vec3 getPosition() const { return position; }
	vec4 getColor() const { return color; }
};

class PointLightQuadric : public Light {
	float intensity;
public:
	PointLightQuadric(vec3 position, vec4 color, float intensity, float r=.1) : Light(position, color, vec4(intensity, 0, 0, 0), 0, r), intensity(intensity) {}
	PointLightQuadric() : PointLightQuadric(vec3(0, 0, 0), vec4(1, 1, 1, 1), 1) {}
	using Light::compressToMatrix;
};

class PointLight: public Light {

public:
	PointLight(vec3 position, float intensity_constant, float intensity_linear, float intensity_quadratic, float softShadowRadius=.1f, vec4 color=WHITE);

	PointLight(vec3 position, vec3 intensities, float r=.1f, vec4 color=WHITE) : PointLight(position, intensities[0], intensities[1], intensities[2], r, color) {}
	PointLight(vec3 position, float intensity_linear, float intensity_quadratic, vec4 color=WHITE) : PointLight(position, 1, intensity_linear, intensity_quadratic, .1f, color) {}
};





enum MeshFormat {
	OBJ = 0
};

std::map<std::string, int> countEstimatedBufferSizesInOBJFile(const char *filename);




class TriangularMesh {
public:
	std::vector <TriangleR3> triangles;

	std::vector <vec3> posBuff;
	std::vector <vec3> normBuff;
	std::vector <vec4> colorBuff;
	std::vector <vec2> uvBuff;

	std::vector <vec3> calculatePositionBuffer() const;
	std::vector <vec3> calculateNormalBuffer() const;
	std::vector <vec4> calculateColorBuffer() const;
	std::vector <vec2> calculateUVBuffer() const;

	std::map<std::string, std::vector<vec4>> extraBuff = {};

    std::array<void*, 4> bufferLocations;
    std::array<int, 4> bufferSizes;

	TriangularMesh();
	explicit TriangularMesh(const std::vector <TriangleR3> &triangles);
	TriangularMesh(std::vector <Vertex> vertices, std::vector <ivec3> faceIndices);
	explicit TriangularMesh(const char* filename, MeshFormat format = OBJ);
	TriangularMesh(ComplexCurve* curve, int nSegments, float h_middle, float w_middle, float w_side);
	std::vector<TriangleR3> getTriangles() const;

	// --- transforms ---
	void translate(vec3 v);
	void transform(const mat4 &M);
	TriangularMesh operator+(vec3 v) const;
	TriangularMesh operator*(const mat4 &M) const;
	TriangularMesh operator*(const mat3 &M) const;
	void operator+=(vec3 v) { this->translate(v); }
	void operator*=(mat4 M) { this->transform(M); }
	void operator*=(mat3 M) { this->transform(mat4(M)); }
    void applyMap(const std::function<vec3(vec3)> &f);
    void applyMap(const std::function<vec3(vec3)> &f, const std::function<vec3(vec3, vec3)> &f_normal);
	void applyMap(const std::function<vec3(vec3)> &f, std::string customDomain);
	void applyMap(std::function<vec3(vec3)> f, std::function<vec3(vec3, vec3)> f_normal, std::string customDomain);

	// --- buffer generators ---

	void precomputeBuffers();
	void precomputeExtraBuffer(std::string name);
	void recomputeBuffers();
	void recomputeBuffers(bool pos, bool norm, bool color, bool uv);

	// --- modifiers ---
	void randomizeFaceColors(vec3 min = vec3(0, 0, 0), vec3 max = vec3(1, 1, 1));
	void cleanUpEmptyTriangles();
	void flipNormals();
	void recalculateNormals();
};




class Model3D {
public:
	std::shared_ptr<TriangularMesh> mesh;
	std::shared_ptr<MaterialPhong> material;
	mat4 transform;
	Model3D();
	Model3D(TriangularMesh &mesh, MaterialPhong &material, const mat4 &transform);
	void addTransform(const mat4 &transform);
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
	std::vector < std::vector<vec2>> boundaries;
	std::vector<bool> boundaryCyclic;
	bool encodesVectorField = false;


	PlanarMeshWithBoundary();
	PlanarMeshWithBoundary(std::vector<TriangleR2> triangles, std::vector < std::vector<vec2>> boundaries, std::vector<bool> boundaryCyclic);
	PlanarMeshWithBoundary(std::vector<TriangleR2> triangles, std::vector < std::vector<vec2>> boundaries);
	PlanarMeshWithBoundary(std::vector<TriangleR2> triangles, std::vector<vec2> boundary, bool cyclic = true);
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
	std::function<std::function<vec3(vec3)>(float, float)> time_operator;
	std::function<std::function<vec3(vec3, vec3)>(float, float)> time_operator_normal;
	bool planar=false;
public:
	MeshFamily1P(std::shared_ptr<TriangularMesh> mesh, std::function<std::function<vec3(vec3)>(float, float)> time_operator, std::function<std::function<vec3(vec3, vec3)>(float, float)> time_operator_normal, float t=0);
	MeshFamily1P(std::shared_ptr<TriangularMesh> embedded_mesh, std::function<PlaneSmoothEndomorphism(float, float)> time_operator, vec3 embedding_shift=vec3(0, 0, 0), float t=0);
	MeshFamily1P(std::shared_ptr<TriangularMesh> embedded_mesh, std::function<PlaneAutomorphism(float)> time_operator_from0, vec3 embedding_shift=vec3(0, 0, 0));

	float time() const;
	virtual void transformMesh(float new_t);
	virtual void meshDeformation(float dt);
};


class MeshFamily1PExtraDomain : public MeshFamily1P {
protected:
	std::string domain;
public:
	MeshFamily1PExtraDomain(std::shared_ptr<TriangularMesh> mesh, std::string domain, std::function<std::function<vec3(vec3)>(float, float)> time_operator, std::function<std::function<vec3(vec3, vec3)>(float, float)> time_operator_normal, float t=0);
	MeshFamily1PExtraDomain(std::shared_ptr<TriangularMesh> embedded_mesh, std::string domain, std::function<PlaneSmoothEndomorphism(float, float)> time_operator, vec3 embedding_shift=vec3(0, 0, 0), float t=0);
	MeshFamily1PExtraDomain(std::shared_ptr<TriangularMesh> embedded_mesh, std::string domain, std::function<PlaneAutomorphism(float)> time_operator_from0, vec3 embedding_shift=vec3(0, 0, 0));

	void transformMesh(float new_t) override;
	void meshDeformation(float dt) override;
};








const BoundaryEmbeddingStyle STD_KERB = BoundaryEmbeddingStyle{KERB, .01f, .02f, .01f, .2};
const BoundaryEmbeddingStyle STD_CURVE = BoundaryEmbeddingStyle{CURVE, .0f, .0f, .0f, .0f, .005f, .005f, .005f};


class SuperMesh {
	std::map<PolyGroupID, std::vector<TriangleR3>> triangleGroups={};
	std::map<PolyGroupID, std::vector<TriangleR3>> boundaryGroups={}; // "bdPoint": extra attribute for unextruded pt
	std::map<PolyGroupID, std::vector<TriangleR3>> embedded_curves={}; // "curvePoint": extra attribute for unextruded pt
	std::map<PolyGroupID, std::vector<TriangleR3>> embedded_points={}; // "loc": extra const attribute for unextruded pt
	std::map<PolyGroupID, MaterialPhongConstColor> materials = std::map<PolyGroupID, MaterialPhongConstColor>();

	MaterialBuffers materialBuffers = {{}, {}, {}, {}};
	StdAttributeBuffers stdAttributeBuffers = {{}, {}, {}, {}};
	std::map<std::string, std::vector<vec4>> extraBuff = {};
	std::map<std::string, int> extraBufferIndices = {};

public:
	std::vector<void*> bufferLocations;
	std::vector<int> bufferSizes;
	SuperMesh();
	explicit SuperMesh(TriangularMesh &mesh);
	SuperMesh(const char* filename, MeshFormat format = OBJ);
	SuperMesh(TriangularMesh &mesh, const MaterialPhongConstColor &material);
	SuperMesh(const char* filename, const MaterialPhongConstColor &material, MeshFormat format = OBJ);
	SuperMesh(PlanarMeshWithBoundary &mesh, const MaterialPhongConstColor &material, MaterialPhongConstColor &material_bd, BoundaryEmbeddingStyle embeddingStyle = STD_KERB);
	SuperMesh(const std::vector <TriangleR3> &triangles, const MaterialPhongConstColor &material);
	SuperMesh(const std::map<PolyGroupID, std::vector<TriangleR3>> &triangleGroups, const std::map<PolyGroupID, MaterialPhongConstColor> &materials);

	void addPolyGroup(PolyGroupID id, const std::vector<TriangleR3> &triangles, const MaterialPhongConstColor &material);
	void addPolyGroup(const std::vector<TriangleR3> &triangles, const MaterialPhongConstColor &material);
	void addPolyGroup(const TriangularMesh &mesh, const MaterialPhongConstColor &material);
	void addBdGroup(PolyGroupID id, const std::vector<TriangleR3> &triangles, const MaterialPhongConstColor &material);
	void addBdGroup(const std::vector<TriangleR3> &triangles, const MaterialPhongConstColor &material);
	void addEmbeddedCurve(PolyGroupID id, const std::vector<TriangleR3> &triangles, const MaterialPhongConstColor &material);
	void addEmbeddedCurve(const std::vector<TriangleR3> &triangles, const MaterialPhongConstColor &material);
	void addEmbeddedPoint(PolyGroupID id, const std::vector<TriangleR3> &triangles, const MaterialPhongConstColor &material);
	void addEmbeddedPoint(const std::vector<TriangleR3> &triangles, const MaterialPhongConstColor &material);
	void embedCurve(ComplexCurve* curve, int nSegments, float h_middle, float w_middle, float w_side, const MaterialPhongConstColor &material);
	void precomputeBuffers(bool materials=true, bool extra=true);
	void merge(const SuperMesh &other);
	void precomputeExtraBuffer(std::string name);

	void actOnPositions(std::function<vec3(vec3)> f);
	void actAtEmbeddedPlane(const Meromorphism &f);
	void actOnPositionsWithCustomShift(std::function<vec3(vec3)> f, std::map<PolyGroupID, std::string> useShiftOfCustomDomain = {});
	void translate(vec3 v);
	void doPerTriangle(std::function<void( TriangleR3&)> f);
	void doPerTriangle(PolyGroupID grp, std::function<void(TriangleR3&)> f);

	void actOnEmbeddedCurve(SpaceEndomorphism f, bool buffOnly=false);


	void randomizeMaterials(MaterialPhongConstColor &min, MaterialPhongConstColor &max);
	void randomiseMaterialsDynamically(float stepMax);
};



class SuperPencilPlanar : public SuperMesh{
	float _t=0;
	std::function<Biholomorphism(float)> _time_operator = [](float t){return Biholomorphism::linear(1, 0);};

public:
	using SuperMesh::SuperMesh;
	void makePencil(std::function<Biholomorphism(float)> _time_operator, float t=0);
	void transformMesh(float new_t);
	void deformMesh(float dt);

};

class CurveSample {
	vec3 position;
	vec3 normal;
	vec3 tangent;
	mat4 material;
	float width;
	vec4 extraInfo = vec4(0, 0, 0, 0);

public:
	CurveSample(vec3 position, vec3 normal, vec3 tangent, MaterialPhongConstColor material, float width) : position(position), normal(normal), tangent(tangent), material(material.compressToMatrix()), width(width) {}
    CurveSample(const CurveSample &other);
    CurveSample(CurveSample &&other) noexcept;
    CurveSample &operator=(const CurveSample &other);
    CurveSample &operator=(CurveSample &&other) noexcept;

    float getWidth() const {return width;}
	vec3 getBinormal() const {return cross(tangent, normal);}
    vec3 getPosition() const {return position;}
    vec3 getNormal() const {return normal;}
    vec3 getTangent() const {return tangent;}
    MaterialPhongConstColor getMaterial() const { return MaterialPhongConstColor(material); }
    mat4 getMaterialMatrix() const { return material; }
	float readExtra(int position) const {return extraInfo[position];}
	float readExtraLast() const	{return extraInfo.w;}
    vec4 readExtra() const {return extraInfo;}
    vec3 readExtraXYZ() const {return vec3(extraInfo);}

	void updatePosition(vec3 new_position) {position = new_position;}
	void updateNormal(vec3 new_normal) {normal = new_normal;}
	void updateTangent(vec3 new_tangent) {tangent = new_tangent;}
	void updateMaterial(const MaterialPhongConstColor &new_material) {material = new_material.compressToMatrix();}
	void updateWidth(float new_width) {width = new_width;}
    void updateExtra(vec4 extra) {extraInfo = extra;}
    void updateExtra(vec3 extra) {extraInfo = vec4(extra, extraInfo.w);}
    void updateExtra(float extra, int position=3) {extraInfo[position] = extra;}
};

std::vector<CurveSample> sampleCurve(SmoothParametricCurve curve, std::function<float(float)> width,
	std::function<MaterialPhongConstColor(float)> material, float t0=0, float t1=1, int n=100, bool periodic=false);



std::vector<CurveSample> sampleCurve(SmoothParametricCurve curve, float width,
                                     MaterialPhongConstColor material, float t0=0, float t1=1, int n=100, bool periodic=false);



enum CurveEmbeddingTypeID {
  PLANAR,
  TUBE,
  NOT_EMBEDDED
};

const std::map<CurveEmbeddingTypeID, std::string> curveEmbeddingTypes = {
        {PLANAR, "PLANAR"}, {TUBE, "TUBE"}, {NOT_EMBEDDED, "NOT EMBEDDED"}};

inline std::string embeddingTypeName(CurveEmbeddingTypeID type) {
  if (curveEmbeddingTypes.contains(type))
    return curveEmbeddingTypes.at(type);
  throw UnknownVariantError("Embedding type with id " + std::to_string(type) + " not recognised as valid variant.");
};


class SuperCurve{
  float t0, t1;
  std::vector<CurveSample> samples;
  std::shared_ptr<SuperMesh> _mesh = nullptr;
  PolyGroupID id;
  CurveEmbeddingTypeID embeddingType = NOT_EMBEDDED;

  SuperMesh generateMeshTube(int radialSegments) const;
  void initMesh() { _mesh = std::make_shared<SuperMesh>(); }

public:
  SuperCurve(const SmoothParametricCurve &curve, const std::function<float(float)> &width, const std::function<MaterialPhongConstColor(float)> &material, int segments, float t0=0, float t1=1,  bool periodic=true);
  SuperCurve(const SmoothParametricCurve &curve, float width, MaterialPhongConstColor material, int nSegments, float t0=0, float t1=1, bool periodic=true);
  // SuperCurve(const SuperCurve &other) : t0(other.t0), t1(other.t1), samples(other.samples), _mesh(other._mesh), _weak_mesh(other._weak_mesh), id(other.id), embeddingType(other.embeddingType) {}
  // SuperCurve(SuperCurve &&other) noexcept : t0(other.domain().x), t1(other.domain().y), samples(std::move(other.samples)), _mesh(std::move(other._mesh)) {}

  void generateMesh(int radialSegments, CurveEmbeddingTypeID type);

  void transformMeshByAmbientMap( const SpaceEndomorphism &f); // just composition
  vec2 domain() const {return vec2(t0, t1);}

  void precomputeBuffers();

  void updateCurve(const SmoothParametricCurve &new_curve);
  std::function<void(float)> pencilDeformerWeak(std::function<SmoothParametricCurve(float)> pencil);
};



class SuperPencilCurve : public SuperCurve {
  float _t = 0;
  std::unique_ptr<End2P> _ambient_operator = nullptr;
  std::unique_ptr<std::function<SmoothParametricCurve(float)>> _parametric_operator = nullptr;
public:
  using SuperCurve::SuperCurve;
  SuperPencilCurve(const SuperCurve &other) : SuperCurve(other) {}
  SuperPencilCurve(SuperCurve &&other) noexcept : SuperCurve(std::move(other)) {}

  void addAmbientDeformation(End2P _ambient_operator, float t=0);
  void addLocalDeformation(End1P _local_operator, float t=0);
  void addDeformationAlongVectorField(VectorField vectorField, float t=0);
  void addPencil(std::function<SmoothParametricCurve(float)> family, float t = 0);
  float time() const;
  void transformMesh(float new_t);
};
