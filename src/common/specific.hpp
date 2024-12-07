#pragma once

// #include "buffer_utils.hpp"
#include <chrono>
#include <iosfwd>
#include <vector>


#include "glsl_utils.hpp"
// #include "renderingUtils.hpp"


const Biholomorphism EXP = Biholomorphism::_EXP();
const Biholomorphism LOG = Biholomorphism::_LOG();
const Biholomorphism IdC = Biholomorphism::linear(ONE, ZERO);
const Biholomorphism ADD1 = Biholomorphism::linear(ONE, ONE);
const Biholomorphism SQUARE = Biholomorphism::power(2);
const Biholomorphism SQRT = Biholomorphism::power(.5f);
const Biholomorphism CAYLEY = Biholomorphism::mobius(Matrix<Complex, 2>(ONE, -I, ONE, I));
const VectorFieldR3 dabbaX = VectorFieldR3::constant(vec3(1, 0, 0));
const VectorFieldR3 dabbaY = VectorFieldR3::constant(vec3(0, 1, 0));
const VectorFieldR3 dabbaZ = VectorFieldR3::constant(vec3(0, 0, 1));






PlanarMeshWithBoundary PlanarUnitDisk(int radial_res, int vertical_res);
PlanarMeshWithBoundary PlanarConvexPolygon(const std::vector<vec2> &verts);
PlanarMeshWithBoundary PlanarRing(int radial_res, int vertical_res, vec2 center, float radiusBig, float radiusSmall);


SmoothParametricPlaneCurve circle(float r, vec2 center=PLANE_ORIGIN, float eps=.01);
SmoothParametricPlaneCurve ellipse(float a, float b, vec2 center=PLANE_ORIGIN, float eps=.01);
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
SmoothParametricCurve trefoil(float r, float R, float eps = .01);
SmoothParametricCurve torusKnot23(float scale, float R, float eps = .01);



SmoothParametricCurve circle(float r, vec3 center, vec3 v1 = e1, vec3 v2 = e2, float eps = .01);
SuperCurve circle(float r, std::function<float(float)> w, const std::function<MaterialPhong(float)> &mat, int n, vec3 center=ORIGIN, vec3 v1=e1, vec3 v2=e2, float eps=.01);

WeakSuperMesh singleTrig(vec3 v0, vec3 v1, vec3 v2, MaterialPhong &material, PolyGroupID id);
WeakSuperMesh singleTrig(vec3 v0, vec3 v1, vec3 v2, MaterialPhong &material1, MaterialPhong &material2, MaterialPhong &material3, PolyGroupID id);
WeakSuperMesh singleTrig(vec3 v0, vec3 v1, vec3 v2, PolyGroupID id);

WeakSuperMesh singleQuadShadeSmooth(vec3 outer1, vec3 inner1, vec3 inner2, vec3 outer2, MaterialPhong &material, PolyGroupID id);
WeakSuperMesh singleQuadShadeSmooth(vec3 outer1, vec3 inner1, vec3 inner2, vec3 outer2,
                                    MaterialPhong &material1in, MaterialPhong &material2in, MaterialPhong &material1out, MaterialPhong &material2out, PolyGroupID id);

WeakSuperMesh singleQuadShadeFlat(vec3 outer1, vec3 inner1, vec3 inner2, vec3 outer2, MaterialPhong &material, PolyGroupID id);
WeakSuperMesh singleQuadShadeFlat(vec3 outer1, vec3 inner1, vec3 inner2, vec3 outer2, MaterialPhong &material1, MaterialPhong &material2, PolyGroupID id);
WeakSuperMesh singleQuad(vec3 outer1, vec3 inner1, vec3 inner2, vec3 outer2, MaterialPhong &materiaInner1, MaterialPhong &materialInner2, MaterialPhong &materialOuter1, MaterialPhong &materialOuter2, bool shadeSmooth, PolyGroupID id);
WeakSuperMesh icosahedron(float r, vec3 center, PolyGroupID id);
WeakSuperMesh icosphere(float r, int n, vec3 center, PolyGroupID id, vec4 color=BLACK);
WeakSuperMesh disk3d(float r, vec3 center, vec3 v1, vec3 v2, int radial_res, int vertical_res, const PolyGroupID &id);

WeakSuperMesh singleQuadShadeFlat(vec3 outer1, vec3 inner1, vec3 inner2, vec3 outer2, PolyGroupID id);




class Disk3D : public WeakSuperMesh {
    vec3 center;
    vec3 forward;
    vec3 down;
    vec3 normal;
    float radius;
    PolyGroupID id;
public:
    Disk3D(const std::vector<Vertex> &nodes, const std::vector<glm::ivec3> &faceInds, vec3 center, vec3 forward, vec3 down, PolyGroupID id);
    Disk3D(const char* filename, vec3 center, vec3 forward, vec3 down, PolyGroupID id);
    Disk3D(float r, vec3 center, vec3 forward, vec3 down, int radial_res, int vertical_res, const PolyGroupID &id);
    void move(vec3 center, vec3 forward, vec3 down, bool scaleWidth);



    float moveRotate(vec3 center, vec3 forward, vec3 down);
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

inline float angleBetween(vec3 vec3, glm::vec3 down) {return  acos(glm::dot(vec3, down)/length(vec3));}

SmoothParametricSurface sphere(float r, vec3 center=ORIGIN, float cutdown=0, float eps=.01);
SmoothParametricSurface DupinCyclide(float a, float b, float d, float eps=.01);
SmoothParametricSurface disk(float r, vec3 center, vec3 v1, vec3 v2, float eps);
SmoothParametricSurface cylinder(float r, vec3 c1, vec3 c2, vec3 v1, vec3 v2, float eps);
SmoothParametricSurface hyperbolic_helicoid(float a, float eps=.01); // todo check out

SmoothParametricSurface LawsonTwist(float alpha, Quaternion q, vec2 range_u, float eps=.01);
SmoothParametricSurface coolLawson(float eps=0.01);

SmoothParametricSurface sudaneseMobius(float eps=0.01);

SmoothParametricSurface LawsonKleinBottle(float eps=0.01);
SmoothParametricSurface catenoid(float a, float eps=.01);
SmoothParametricSurface helicoid(float a, float eps=.01);
SmoothParametricSurface enneper(float k, float l, float eps=.01);
SmoothParametricSurface twistedTorus(float a, float m, float n, int dommul1, int dommul2, float eps=.01);

SmoothParametricSurface polyTorus(float r, float n, float alpha, float eps=.01);
SmoothParametricSurface ellipticTorus(float c, float eps=.01);
SmoothParametricSurface limpetTorus(float eps=.01);
SmoothParametricSurface fig8(float c, float eps=.01);
SmoothParametricSurface doubleTorus(float eps=.01);
SmoothParametricSurface saddleTorus(float eps=.01);
SmoothParametricSurface kinkyTorus(float eps=.01);
SmoothParametricSurface GraysKlein(float a, float n, float m, float eps=.01);
SmoothParametricSurface bowTie(float eps=.01);
SmoothParametricSurface bohemianDome(float a, float b, float c, float eps=.01);
SmoothParametricSurface horn(float eps=.01);
SmoothParametricSurface crescent(float eps=.01);
SmoothParametricSurface seaShell(int n, float a, float b, float c, float eps=.01);

SmoothParametricSurface cone(const SmoothParametricCurve &base, vec3 apex, float eps);
SmoothParametricSurface coneSide(float r, float h, vec3 center, vec3 v1, vec3 v2, float eps);

WeakSuperMesh arrow(vec3 start, vec3 head, float radius, float head_len, float head_radius, int radial, int straight, float eps, std::variant<int, std::string> id);
WeakSuperMesh drawArrows(const vector<vec3> &points, const vector<vec3> &directions, float radius, float head_len, float head_radius, int radial, int straight, float eps, const std::variant<int, std::string> &id);
WeakSuperMesh drawVectorFieldArrows(const VectorFieldR3 &field, const vector<vec3> &points, const HOM(float, float)& len, const HOM(float, float) &radius, const HOM(float, float) &head_len, const HOM(float, float)& head_radius, int radial, int straight, float eps, const std::variant<int, std::string> &id);
vec3 getArrayHead(const BufferedVertex &v);
vec3 getArrayStart(const BufferedVertex &v);
vec3 getArrayDirection(const BufferedVertex &v);
vec3 getArrayHead(PolyGroupID& id, WeakSuperMesh &mesh);
vec3 getArrayStart(PolyGroupID &id, WeakSuperMesh &mesh);
vec3 getArrayDirection(PolyGroupID &id, WeakSuperMesh &mesh);

inline void setArrayStart(BufferedVertex &v, vec3 head) {v.setColor(head.x, 0); v.setColor(head.y, 1); v.setColor(head.z, 2);}
inline void setArrayHead(BufferedVertex &v, vec3 start) {v.setUV(vec2(start.x, start.y)); v.setColor(start.z, 3);}


SmoothParametricCurve PLCurve(std::vector<vec3> points);
SmoothParametricCurve segment(vec3 p0, vec3 p1, float t0=0, float t1=1);

inline SmoothParametricCurve segment(vec3 p0, vec3 p1, float t0, float t1) {
	return SmoothParametricCurve([p0, p1, t0, t1](float t) { return (p1*(t-t0) + p0*(t1-t))/(t1-t0); }, "seg", t0, t1, false, .01);
}

SmoothParametricPlaneCurve GernsterWave(float a, float b, float k, float c);
VectorFieldR3 PousevillePlanarFlow(float h, float nabla_p, float mu, float v0, float eps=.01);
VectorFieldR3 PousevillePipeFlow(float nabla_p, float mu, float c1, float c2, float eps=.01);

vec3 freeSurfaceComponent(float amplitude, float phase, vec2 waveNumber, float angFrequency, float h, vec2 ab, float t);
SurfaceParametricPencil freeSurface(vector<float> a_m, vector<float> phi_m, vector<vec2> k_m, vector<float> omega_m, float h);



WeakSuperMesh particles(int n, vec3 bound1, vec3 bound2, float radius);
WeakSuperMesh box(vec3 size, vec3 center,PolyGroupID id);
