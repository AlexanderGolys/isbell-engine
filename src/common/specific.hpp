#pragma once

// #include "buffer_utils.hpp"
#include <glm/detail/_vectorize.hpp>

#include "glsl_utils.hpp"
// #include "renderingUtils.hpp"


const Biholomorphism EXP = Biholomorphism::_EXP();
const Biholomorphism LOG = Biholomorphism::_LOG();
const Biholomorphism IdC = Biholomorphism::linear(ONE, ZERO);
const Biholomorphism ADD1 = Biholomorphism::linear(ONE, ONE);
const Biholomorphism SQUARE = Biholomorphism::power(2);
const Biholomorphism SQRT = Biholomorphism::power(.5f);
const Biholomorphism CAYLEY = Biholomorphism::mobius(Matrix<Complex, 2>(ONE, -I, ONE, I));
const VectorFieldR3 dabbaX = VectorFieldR3::constant(glm::vec3(1, 0, 0));
const VectorFieldR3 dabbaY = VectorFieldR3::constant(glm::vec3(0, 1, 0));
const VectorFieldR3 dabbaZ = VectorFieldR3::constant(glm::vec3(0, 0, 1));






PlanarMeshWithBoundary PlanarUnitDisk(int radial_res, int vertical_res);
PlanarMeshWithBoundary PlanarConvexPolygon(const std::vector<glm::vec2> &verts);
PlanarMeshWithBoundary PlanarRing(int radial_res, int vertical_res, glm::vec2 center, float radiusBig, float radiusSmall);


SmoothParametricPlaneCurve circle(float r, glm::vec2 center=PLANE_ORIGIN, float eps=.01);
SmoothParametricPlaneCurve ellipse(float a, float b, glm::vec2 center=PLANE_ORIGIN, float eps=.01);
SmoothParametricPlaneCurve epitrochoid(float r, float R, float d, float eps = .01);
SmoothParametricPlaneCurve epicycloid(float r, float R, float eps = .01);
SmoothParametricPlaneCurve cardioid(float r, float eps=.01);
SmoothParametricPlaneCurve nephroid(float r, float eps=.01);
SmoothParametricPlaneCurve trefoiloid(float r, float eps=.01);
SmoothParametricPlaneCurve quatrefoloid(float r, float eps=.01);
SmoothParametricPlaneCurve hypotrochoid(float r, float R, float d, float eps = .01);
SmoothParametricPlaneCurve hypocycloid(float r, float R, float eps = .01);
SmoothParametricPlaneCurve astroid(float r, float eps=.01);
SmoothParametricPlaneCurve deltoid(float r, float eps=.01);
SmoothParametricPlaneCurve pentoid(float r, float eps=.01);
SmoothParametricPlaneCurve exoid(float r, float eps=.01);
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
WeakSuperMesh icosahedron(float r, glm::vec3 center, PolyGroupID id);
WeakSuperMesh icosphere(float r, int n, glm::vec3 center, PolyGroupID id);
WeakSuperMesh disk3d(float r, glm::vec3 center, glm::vec3 v1, glm::vec3 v2,int radial_res, int vertical_res, const PolyGroupID &id);

class Disk3D : public WeakSuperMesh {
    glm::vec3 center;
    glm::vec3 forward;
    glm::vec3 down;
    glm::vec3 normal;
    float radius;
    PolyGroupID id;
public:
    Disk3D(const std::vector<Vertex> &nodes, const std::vector<glm::ivec3> &faceInds, glm::vec3 center, glm::vec3 forward, glm::vec3 down, PolyGroupID id);
    Disk3D(const char* filename, glm::vec3 center, glm::vec3 forward, glm::vec3 down, PolyGroupID id);
    Disk3D(float r, glm::vec3 center, glm::vec3 forward, glm::vec3 down, int radial_res, int vertical_res, const PolyGroupID &id);
    void move(glm::vec3 center, glm::vec3 forward, glm::vec3 down, bool scaleWidth);



    float moveRotate(glm::vec3 center, glm::vec3 forward, glm::vec3 down);
    void rotate(float angle);
    static float angle(const BufferedVertex &v) { return v.getColor().x; }
    static float rParam(const BufferedVertex &v) { return v.getColor().y; }
    static float width(const BufferedVertex &v) { return v.getColor().z; }
    static float widthNormalised(const BufferedVertex &v) { return v.getColor().w; }
    float getR() const {return radius;}

    float rReal(const BufferedVertex &v);
    void scaleR(float r, bool scaleWidth);
    void setR(float r);;
    static void setAbsoluteWidth(BufferedVertex &v, float w) {v.setColor(w, 2);}
    static void setRelativeWidth(BufferedVertex &v, float w) {v.setColor(w, 3);}
    void setEmpiricalRadius();
    void setColorInfo();
};

inline float angleBetween(glm::vec3 vec3, glm::vec3 down) {return  acos(glm::dot(vec3, down)/glm::length(vec3));}

SmoothParametricSurface sphere(float r, glm::vec3 center=ORIGIN, float cutdown=0, float eps=.01);
