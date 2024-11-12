#pragma once

#include "indexedGeometry.hpp"


const Biholomorphism EXP = Biholomorphism::_EXP();
const Biholomorphism LOG = Biholomorphism::_LOG();
const Biholomorphism Id = Biholomorphism::linear(ONE, ZERO);
const Biholomorphism ADD1 = Biholomorphism::linear(ONE, ONE);
const Biholomorphism SQUARE = Biholomorphism::power(2);
const Biholomorphism SQRT = Biholomorphism::power(.5f);
const Biholomorphism CAYLEY = Biholomorphism::mobius(Matrix<Complex, 2>(ONE, -I, ONE, I));
const VectorFieldR3 dabbaX = VectorFieldR3::constant(glm::vec3(1, 0, 0));
const VectorFieldR3 dabbaY = VectorFieldR3::constant(glm::vec3(0, 1, 0));
const VectorFieldR3 dabbaZ = VectorFieldR3::constant(glm::vec3(0, 0, 1));

#define Mob Matrix<Complex, 2>

const Mob CayleyTransform = Mob(ONE, -I, ONE, I);
const Mob CayleyTransformInv = ~Mob(ONE, -I, ONE, I);
const Mob Imob = Mob(1, ZERO, ZERO, 1);


class BoxMesh : public TriangularMesh {
public:
	BoxMesh(glm::vec3 minCorner, glm::vec3 maxCorner);
	BoxMesh(glm::vec3 minCorner, glm::vec3 maxCorner, glm::vec4 color);
	BoxMesh(glm::vec3 minCorner, glm::vec3 maxCorner, std::array<glm::vec4, 6> faceColors);


	void assignFaceColors(std::array<glm::vec4, 6> faceColors);
	void randomizeFaceColors();
};

class BallMesh : public TriangularMesh {
public:
	glm::vec3 center;
	float radius;
	std::vector<glm::vec3> bdVertices;
	BallMesh(glm::vec3 center, float radius, int radialSegments, int verticalSegments, glm::vec3 normal);
	BallMesh(glm::vec3 center, float radius, int radialSegments, int verticalSegments=3);
	BallMesh(int radialSegments, int verticalSegments=3);
	BallMesh();
	TriangularMesh extrudeAlongNormal(float h);
};

class CircularRing : public TriangularMesh {
public:
	glm::vec3 center;
	float radiusBig;
	float radiusSmall;
	std::vector<glm::vec3> bdVerticesInner;
	std::vector<glm::vec3> bdVerticesOuter;


	CircularRing(glm::vec3 center, float radiusBig, float radiusSmall, int radialSegments, int verticalSegments, glm::vec3 normal);
	CircularRing(glm::vec3 center, float radiusBig, float radiusSmall, int radialSegments, int verticalSegments = 3);
	CircularRing(int radialSegments, int verticalSegments = 3);
	CircularRing();
	TriangularMesh extrudeAlongNormal(float h);
};


class PlanarUnitDisk : public PlanarMeshWithBoundary {
public:
	PlanarUnitDisk(int radial_res, int vertical_res);
};

class PlanarConvexPolygon: public PlanarMeshWithBoundary{
public:
	std::vector<glm::vec2> vertices;
	PlanarConvexPolygon(std::vector<glm::vec2> verts);
};

class PlanarRing : public PlanarMeshWithBoundary {
public:
	glm::vec2 center;
	float radiusBig;
	float radiusSmall;
	PlanarRing(int radial_res, int vertical_res, glm::vec2 center, float radiusBig, float radiusSmall);
	void addRotationalField(float power);
};

class COLOR_PALETTE {
public:
	glm::vec4 mainColor;
	glm::vec4 second;
	glm::vec4 third;
	glm::vec4 accent;
	glm::vec4 accent2;

	COLOR_PALETTE(glm::vec4 mainColor, glm::vec4 second, glm::vec4 third, glm::vec4 accent, glm::vec4 accent2);
	COLOR_PALETTE(glm::vec3 mainColor, glm::vec3 second, glm::vec3 third, glm::vec3 accent, glm::vec3 accent2);
	COLOR_PALETTE(glm::ivec3 mainColor, glm::ivec3 second, glm::ivec3 third, glm::ivec3 accent, glm::ivec3 accent2);
	std::vector<glm::vec4> colors();
	glm::vec4 operator[] (int i);

};

class COLOR_PALETTE10 {
public:
	std::array<glm::vec4, 10> cls;

	COLOR_PALETTE10(std::array<glm::vec4, 10> colors) : cls(colors) {}
	COLOR_PALETTE10(COLOR_PALETTE p1, COLOR_PALETTE p2) : cls({p1[0], p1[1], p1[2], p1[3], p1[4], p2[0], p2[1], p2[2], p2[3], p2[4]}) {}
	COLOR_PALETTE10(glm::ivec3 c1, glm::ivec3 c2, glm::ivec3 c3, glm::ivec3 c4, glm::ivec3 c5, glm::ivec3 c6, glm::ivec3 c7, glm::ivec3 c8, glm::ivec3 c9, glm::ivec3 c10);
	std::vector<glm::vec4> colors() const { return std::vector<glm::vec4>(cls.begin(), cls.end()); }
	glm::vec4 operator[] (int i) const { return cls[i]; }
};

SmoothParametricPlaneCurve circle(float r, glm::vec2 center=PLANE_ORIGIN, float eps=.01);
SmoothParametricPlaneCurve ellipse(float a, float b, glm::vec2 center=PLANE_ORIGIN, float eps=.01);
SmoothParametricPlaneCurve epitrochoid(float r, float R, float d, float eps = .01);
inline SmoothParametricPlaneCurve epicycloid(float r, float R, float eps = .01) { return epitrochoid(r, R, r, eps); }
inline SmoothParametricPlaneCurve cardioid(float r, float eps=.01) { return epicycloid(r, r, eps); }
inline SmoothParametricPlaneCurve nephroid(float r, float eps=.01) { return epicycloid(r, 2*r, eps); }
inline SmoothParametricPlaneCurve trefoiloid(float r, float eps=.01) { return epicycloid(r, 3*r, eps); }
inline SmoothParametricPlaneCurve quatrefoloid(float r, float eps=.01) { return epicycloid(r, 4*r, eps); }
SmoothParametricPlaneCurve hypotrochoid(float r, float R, float d, float eps = .01);
inline SmoothParametricPlaneCurve hypocycloid(float r, float R, float eps = .01) { return hypotrochoid(r, R, r, eps); }
inline SmoothParametricPlaneCurve astroid(float r, float eps=.01) { return hypocycloid(r, 4*r, eps); }
inline SmoothParametricPlaneCurve deltoid(float r, float eps=.01) { return hypocycloid(r, 3*r, eps); }
inline SmoothParametricPlaneCurve pentoid(float r, float eps=.01) { return hypocycloid(r, 5*r, eps); }
inline SmoothParametricPlaneCurve exoid(float r, float eps=.01) { return hypocycloid(r, 6*r, eps); }
SmoothParametricPlaneCurve LissajousCurve(float a, float b, float delta=PI/2, float r1=1, float r2=1, float eps = .01);
SmoothParametricCurve VivaniCurve(float r, float eps = .01);
SmoothParametricCurve sphericalSpiral(float a, float t_max, PolyGroupID id, float eps = .001);
SmoothParametricCurve sphericalSpiral(float a, float r, float t_max, PolyGroupID id, float eps = .001);


SmoothParametricCurve circle(float r, glm::vec3 center, glm::vec3 v1 = e1, glm::vec3 v2 = e2, float eps = .01);
SuperCurve circle(float r, std::function<float(float)> w, const std::function<MaterialPhong(float)> &mat, int n, glm::vec3 center=ORIGIN, glm::vec3 v1=e1, glm::vec3 v2=e2, float eps=.01);


WeakSuperMesh singleTrig(glm::vec3 v0, glm::vec3 v1, glm::vec3 v2, MaterialPhong &material, PolyGroupID id);
WeakSuperMesh singleTrig(glm::vec3 v0, glm::vec3 v1, glm::vec3 v2, MaterialPhong &material1, MaterialPhong &material2, MaterialPhong &material3, PolyGroupID id);

WeakSuperMesh singleQuadShadeSmooth(glm::vec3 outer1, glm::vec3 inner1, glm::vec3 inner2, glm::vec3 outer2, MaterialPhong &material, PolyGroupID id);

WeakSuperMesh singleQuadShadeSmooth(glm::vec3 outer1, glm::vec3 inner1, glm::vec3 inner2, glm::vec3 outer2,
                                    MaterialPhong &material1in, MaterialPhong &material2in, MaterialPhong &material1out, MaterialPhong &material2out, PolyGroupID id);

WeakSuperMesh singleQuadShadeFlat(glm::vec3 outer1, glm::vec3 inner1, glm::vec3 inner2, glm::vec3 outer2, MaterialPhong &material, PolyGroupID id);
WeakSuperMesh singleQuadShadeFlat(glm::vec3 outer1, glm::vec3 inner1, glm::vec3 inner2, glm::vec3 outer2, MaterialPhong &material1, MaterialPhong &material2, PolyGroupID id);
WeakSuperMesh singleQuad(glm::vec3 outer1, glm::vec3 inner1, glm::vec3 inner2, glm::vec3 outer2, MaterialPhong &materiaInner1, MaterialPhong &materialInner2, MaterialPhong &materialOuter1, MaterialPhong &materialOuter2, bool shadeSmooth, PolyGroupID id);
WeakSuperMesh icosahedron(float r, glm::vec3 center, MaterialPhong &material, PolyGroupID id);
WeakSuperMesh icosphere(float r, int n, glm::vec3 center, MaterialPhong &material, PolyGroupID id);

SmoothParametricSurface sphere(float r, glm::vec3 center=ORIGIN, float cutdown=0, float eps=.01);
