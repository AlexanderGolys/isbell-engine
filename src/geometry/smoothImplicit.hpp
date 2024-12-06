#pragma once
#include <array>

#include <iosfwd>
#include <utility>
#include <vector>

//#include "src/common/indexedRendering.hpp"

#include "smoothParametric.hpp"

class AffineLine;
class TriangulatedImplicitSurface;
class WeakSuperMesh;
// R3 -> R
class SmoothImplicitSurface {
    RealFunctionR3 _F;
public:
    explicit SmoothImplicitSurface(const RealFunctionR3 &F);
    float operator()(vec3 p) const;
};

class AffinePlane : public SmoothImplicitSurface {
    vec3 n;
    float d; // (n, x) - d = 0
    vec3 pivot, v1, v2;


public:


    AffinePlane(vec3 n, float d);
    AffinePlane(vec3 pivot, vec3 v1, vec3 v2);
    static AffinePlane spanOfPts(vec3 p0, vec3 p1, vec3 p2);
    AffineLine intersect(const AffinePlane &p) const;
    float distance (vec3 p) const { return dot(p, n) - d; }
    vec3 orthogonalProjection(vec3 p) const { return p - distance(p) * n; }
    bool contains(vec3 p, float eps=1e-6) const { return abs(distance(p)) < eps; }

    std::pair<vec3, float> equationCoefs() const { return {n, d}; }
    vec3 intersect(const AffineLine &l) const;
    RealFunctionR3 distanceField() const {return RealFunctionR3([this](vec3 p) {return distance(p);} ); }
    mat3 pivotAndBasis() const;
    vec2 localCoordinates(vec3 p) const;

    vec3 normal() const { return n; }
    float getD() const { return d; }
};



class ImplicitSurfacePoint {
	vec3 p;
	mat3 orthoFrame;
	int birth_time;
	float angle = 0;
	bool requiresRecalculation = true;


public:
	bool border = false;

	ImplicitSurfacePoint(vec3 p, const mat3 &frame, int birth_time, bool border) : p(p), orthoFrame(frame), birth_time(birth_time), border(border) {}


	vec3 projectOnTangentPlane(vec3 q) const;
	vec3 rotateAroundNormal(vec3 q, float angle) const;

	vec3 coords_in_frame(vec3 v) const { return inverse(orthoFrame) * v; }
	vec3 getPosition() const { return p; }
	vec3 getNormal() const { return orthoFrame[2]; }
	std::pair<vec3, vec3> getTangents() const { return {orthoFrame[0], orthoFrame[1]}; }
	vec3 getTangent1() const { return  orthoFrame[0]; }
	vec3 getTangent2() const { return  orthoFrame[1]; }

	void markForRecalculation() { requiresRecalculation = true; }
	void markAsCalculated() { requiresRecalculation = false; }
	bool needsRecalculation() const { return requiresRecalculation; }
	float getAngle() const { return angle; }
	void setAngle(float a) { angle = a; requiresRecalculation = false; }
};


class TriangulatedImplicitSurface {
	RealFunctionR3 F;
	HOM(vec3, bool) borderCheck;
	int max_iter;
	int NewtonMaxSteps;
	float NewtonEps;
	int age = 0;
	std::vector<ImplicitSurfacePoint> points = {};
	vector<ivec3> triangles = {};
	vector<vector<int>> polygons = {};

	float angle(int i);

public:
	TriangulatedImplicitSurface(const RealFunctionR3 &F, int max_iter, HOM(vec3, bool) borderCheck, int NewtonMaxSteps, float NewtonEps);

	ImplicitSurfacePoint projectOnSurface(vec3 p);
	ImplicitSurfacePoint constructPoint(vec3 p);

	void initialiseHexagon(vec3 p0, float len);
	void expandFrontPolygon();
	void calculateMissingAngles();
	void calculateAngle(int i);
	int minAngleIndex();

	void checkFrontPolygonSelfIntersections();
	void checkFrontPolygonIntersectionWithPolygon(int i);
	void expandFrontPolygon(int i);
	void step();

	void generate();

};
