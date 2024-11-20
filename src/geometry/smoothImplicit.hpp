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
    float operator()(glm::vec3 p) const;
};

class AffinePlane : public SmoothImplicitSurface {
    glm::vec3 n;
    float d; // (n, x) - d = 0
    glm::vec3 pivot, v1, v2;
public:
    AffinePlane(glm::vec3 n, float d);
    AffinePlane(glm::vec3 pivot, glm::vec3 v1, glm::vec3 v2);
    static AffinePlane spanOfPts(glm::vec3 p0, glm::vec3 p1, glm::vec3 p2);
    AffineLine intersect(const AffinePlane &p) const;
    float distance (glm::vec3 p) const { return dot(p, n) - d; }
    glm::vec3 orthogonalProjection(glm::vec3 p) const { return p - distance(p) * n; }
    bool contains(glm::vec3 p, float eps=1e-6) const { return abs(distance(p)) < eps; }

    std::pair<glm::vec3, float> equationCoefs() const { return {n, d}; }
    glm::vec3 intersect(const AffineLine &l) const;
    RealFunctionR3 distanceField() const {return RealFunctionR3([this](glm::vec3 p) {return distance(p);} ); }
    glm::mat3 pivotAndBasis() const;
    glm::vec2 localCoordinates(glm::vec3 p) const;

    glm::vec3 normal() const { return n; }
    float getD() const { return d; }
};


class ImplicitSurfacePoint {
	vec3 p;
	mat3 orthoFrame;
	bool border;
	float front_angle = 2137;
	bool angle_changed = false;
	int birth_time;
public:
	ImplicitSurfacePoint(vec3 p, const mat3 &orthoFrame, bool border, int birth_time) : p(p), orthoFrame(orthoFrame), border(border), birth_time(birth_time){}
	ImplicitSurfacePoint(vec3 p, vec3 normal, bool border);
	vec3 projectOnTangentPlane(vec3 q) const;
	vec3 rotateAroundNormal(vec3 q, float angle) const;

	vec3 coords_in_frame(vec3 v) const { return inverse(orthoFrame) * v; }
	vec3 getPosition() const { return p; }
	vec3 getNormal() const { return orthoFrame[2]; }
	float getAngle() const { return front_angle; }
	std::pair<vec3, vec3> getTangents() const { return {orthoFrame[0], orthoFrame[1]}; }
	vec3 getTangent1() const { return  orthoFrame[0]; }
	vec3 getTangent2() const { return  orthoFrame[1]; }
	void setAngle(float angle) { front_angle = angle; }
	void angleChanged() { angle_changed = true; }
	bool isBorder() const { return border; }
	void angleRecomputed() { angle_changed = false; }
};


class FrontPolygon {
	std::vector<ImplicitSurfacePoint> points;
	float side;
	std::vector<int> excluded_for_distance_check = {};
public:
	FrontPolygon(std::vector<ImplicitSurfacePoint> points, float side) : points(std::move(points)), side(side) {}
	void recalculate_angle(int index);
	void recalculate_angles();
	int argminAngle() const;
	float minAngle() const;
	void addPoint(ImplicitSurfacePoint p, int index);

	bool removePoint(int index);

	void expandVertex(int index, TriangulatedImplicitSurface &surface);
	void step(TriangulatedImplicitSurface &surface);
	bool checkForSelfIntersections();
	bool merge(FrontPolygon &other, int ind_self, int ind_other);
	bool checkForCrossIntersections(const FrontPolygon &other);
	int size() const { return points.size(); }

	void noExcluded() { excluded_for_distance_check.clear(); }
};


class TriangulatedImplicitSurface {
	RealFunctionR3 F;
	int max_iter = 1;
	std::vector<ImplicitSurfacePoint> points = {};
	std::vector<std::array<vec3, 3>> triangles = {};
	std::vector<FrontPolygon> polygons = {};
	int current_step = 0;
	float eps = 0.0001;
	int front_polygon = 0;
public:
	explicit TriangulatedImplicitSurface(RealFunctionR3 F) : F(std::move(F)){}
	void removeDegeneratePolygons();
	ImplicitSurfacePoint findNearbyPoint(glm::vec3 p) const;
	void initHexagon(vec3 p, float side);

	bool step();
	void addTriangle(const std::array<vec3, 3> &t) { triangles.push_back(t); }
	WeakSuperMesh compute(int max_iter, vec3 p0, float side);
};
