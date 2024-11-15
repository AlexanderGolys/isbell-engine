#pragma once
#include "smoothParametric.hpp"

class AffineLine;
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
