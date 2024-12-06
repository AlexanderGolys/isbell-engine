#pragma once

#include "src/common/indexedRendering.hpp"


class TriangulatedManifold : public WeakSuperMesh {
	vector<ivec2> edgesVert;
	vector<ivec3> facesEdg;
	PolyGroupID id;

public:
	vec3 vertexPosition(int i) const { return getVertices(id)[i].getPosition(); }
	mat2x3 edgeVertices (int i) const { return {
	mat2x3(vertexPosition(edgesVert[i].x), vertexPosition(edgesVert[i].y))}; }
	mat3 faceVertices(int i);
	vec3 edgeVector(int i) const { return vertexPosition(edgesVert[i].y) - vertexPosition(edgesVert[i].x); }
	mat3 edgesVectorsOfFace(int i) const { return mat3(edgeVector(facesEdg[i].x), edgeVector(facesEdg[i].y), edgeVector(facesEdg[i].z)); }
	mat3 tangentNormalFrameOfFace(int i) const {return GramSchmidtProcess(edgesVectorsOfFace(i));}

	ivec2 edgeVerticesIndices(int i) const { return edgesVert[i]; }
	ivec3 faceEdgesIndices(int i) const { return facesEdg[i]; }

	mat2x3 orthoFaceTangents(int i) const;
	mat2x3 stdFaceTangents(int i) const;
	vec2 stdTangentToOrtho(vec2 v, int i) const;
	vec2 orthoTangentToStd(vec2 v, int i) const;

	int numEdges() const { return edgesVert.size(); }
	int numFaces() const { return facesEdg.size(); }
	int numVertices() const { return getVertices(id).size(); }
};


class Discrete1Form;


class DiscreteRealFunction {
	BigVector values;
	std::shared_ptr<TriangulatedManifold> domain;
public:
	explicit DiscreteRealFunction(const vector<float> &values, const std::shared_ptr<TriangulatedManifold> &domain) : values(values), domain(domain) {}
	explicit DiscreteRealFunction(const BigVector &values, const std::shared_ptr<TriangulatedManifold> &domain) : values(values), domain(domain) {}

	float operator()(int i) { return values[i]; }
	DiscreteRealFunction operator+(const DiscreteRealFunction &f) const { return DiscreteRealFunction(values + f.values, domain); }
	DiscreteRealFunction operator-(const DiscreteRealFunction &f) const { return DiscreteRealFunction(values - f.values, domain); }
	DiscreteRealFunction operator*(float a) const { return DiscreteRealFunction(values*a, domain); }
	DiscreteRealFunction operator/(float a) const { return DiscreteRealFunction(values/a, domain); }
	DiscreteRealFunction operator-() const { return DiscreteRealFunction(-values, domain); }
	BigVector toBumpBasis() const;
	static DiscreteRealFunction bump(int i);

	float edgeLength(int i) const { return length(domain->edgeVector(i)); }
	float faceArea(int i) const { return abs(det(domain->faceVertices(i))); }
	Discrete1Form exteriorDerivative() const;
};

class Discrete2Form;

class Discrete1Form {
	BigVector edgeValues;
	std::shared_ptr<TriangulatedManifold> domain;

public:
	explicit Discrete1Form(const vector<float> &values, const std::shared_ptr<TriangulatedManifold> &domain) : edgeValues(values), domain(domain) {}
	explicit Discrete1Form(const BigVector &values, const std::shared_ptr<TriangulatedManifold> &domain) :  edgeValues(values), domain(domain) {}

	float operator()(int i) { return  edgeValues[i]; }
	Discrete1Form operator+(const Discrete1Form &f) const { return Discrete1Form( edgeValues + f. edgeValues, domain); }
	Discrete1Form operator-(const Discrete1Form &f) const { return Discrete1Form( edgeValues - f. edgeValues, domain); }
	Discrete1Form operator*(float a) const { return  Discrete1Form( edgeValues*a, domain); }
	Discrete1Form operator/(float a) const { return Discrete1Form( edgeValues/a, domain); }
	Discrete1Form operator-() const { return Discrete1Form(- edgeValues, domain); }
	BigVector toBumpOrthoBasis() const;
	static Discrete1Form bumpOrtho(int i, int j);
	static Discrete1Form bumpStd(int i, int j);
	float integrate(const vector<int> &edgePath) const;
	Discrete2Form exteriorDerivative() const;
//	Discrete2Form d() const { return exteriorDerivative(); }
};

class Discrete2Form {
	BigVector faceValues;
	std::shared_ptr<TriangulatedManifold> domain;
public:
	explicit Discrete2Form(const vector<float> &values, const std::shared_ptr<TriangulatedManifold> &domain) : faceValues(values), domain(domain) {}
	explicit Discrete2Form(const BigVector &values, const std::shared_ptr<TriangulatedManifold> &domain) : faceValues(values), domain(domain) {}
	Discrete2Form operator+(const Discrete2Form &f) const { return Discrete2Form( faceValues + f. faceValues, domain); }
	Discrete2Form operator-(const Discrete2Form &f) const { return Discrete2Form( faceValues - f. faceValues, domain); }
	Discrete2Form operator*(float a) const { return  Discrete2Form( faceValues*a, domain); }
	Discrete2Form operator/(float a) const { return Discrete2Form(  faceValues/a, domain); }
	Discrete2Form operator-() const { return Discrete2Form(- faceValues, domain); }
	float integrate();

	BigVector toBumpOrthoBasis() const;
	static Discrete2Form bumpOrtho(int i, int j);
	static Discrete2Form bumpStd(int i, int j);
};
