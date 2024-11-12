#pragma once

#include "func.hpp"



// R2 -> R3
class SmoothParametricSurface {
  Foo113 _f;
  Foo113 _df_t;
  Foo113 _df_u;
  float t0, t1, u0, u1;
  bool t_periodic, u_periodic;
    float epsilon=.01;

public:
  SmoothParametricSurface(Foo113 f, Foo113 df_t, Foo113 df_u, glm::vec2 t_range, glm::vec2 u_range, bool t_periodic=false, bool u_periodic=false);
  SmoothParametricSurface(Foo113 f, glm::vec2 t_range, glm::vec2 u_range, bool t_periodic=false, bool u_periodic=false, float epsilon=.01);
  SmoothParametricSurface(std::function<SmoothParametricCurve(float)> pencil, glm::vec2 t_range, glm::vec2 u_range, bool t_periodic=false, bool u_periodic=false, float eps=.01);
  SmoothParametricSurface(SmoothParametricCurve canal, Fooo width, float eps=.01);
  SmoothParametricSurface(SmoothParametricCurve pipe, float width, float eps=.01) : SmoothParametricSurface(pipe, [width](float t) {return width; }, eps) {}


  glm::vec3 operator()(float t, float s) const;
  glm::vec3 operator()(glm::vec2 tu) const;
  glm::vec3 parametersNormalised(glm::vec2 tu) const { return operator()(t0 + tu.x*(t1-t0), u0 + tu.y*(u1-u0)); }
  glm::vec3 parametersNormalised(float t, float u) const { return operator()(t0 + t*(t1-t0), u0 + u*(u1-u0)); }
  SmoothParametricCurve precompose(SmoothParametricPlaneCurve c, PolyGroupID id=420) const;
  SmoothParametricCurve restrictToInterval(glm::vec2 p0, glm::vec2 p1, PolyGroupID id=420) const;

  glm::vec2 boundsT() const { return glm::vec2(t0, t1); }
  glm::vec2 boundsU() const { return glm::vec2(u0, u1); }
  float tMin() const { return t0; }
  float tMax() const { return t1; }
  float uMin() const { return u0; }
  float uMax() const { return u1; }
  bool isPeriodicT() const { return t_periodic; }
  bool isPeriodicU() const { return u_periodic; }
  float periodT() const { return t_periodic ? t1 - t0 : 0; }
  float periodU() const { return u_periodic ? u1 - u0 : 0; }

  glm::vec3 normal(float t, float s) const;
  glm::vec3 normal(glm::vec2 tu) const { return normal(tu.x, tu.y); }


};


// R3 -> R
class SmoothImplicitSurface {
  SmoothRealFunctionR3 _F;
public:
  SmoothImplicitSurface(SmoothRealFunctionR3 F);
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
  AffineLine intersection(AffinePlane &p) const;
  float distance (glm::vec3 p) const { return dot(p, n) - d; }
  glm::vec3 orthogonalProjection(glm::vec3 p) const { return p - distance(p) * n; }
  bool contains(glm::vec3 p, float eps=1e-6) const { return abs(distance(p)) < eps; }

  std::pair<glm::vec3, float> equationCoefs() const { return {n, d}; }
  glm::vec3 intersection(AffineLine &l) const;
  SmoothRealFunctionR3 distanceField() const {return SmoothRealFunctionR3([this](glm::vec3 p) {return distance(p);} ); }
  glm::mat3 pivotAndBasis() const;
  glm::vec2 localCoordinates(glm::vec3 p) const;

  glm::vec3 normal() const { return n; }
  float getD() const { return d; }
};

