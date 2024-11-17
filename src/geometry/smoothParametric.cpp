#include "smoothParametric.hpp"
#include <map>

using namespace glm;
using std::vector, std::string, std::map, std::shared_ptr, std::unique_ptr, std::pair, std::make_unique, std::make_shared, std::function;

SmoothParametricSurface::SmoothParametricSurface(std::function<vec3(float, float)> f, std::function<vec3(float, float)> df_t,
                                                 std::function<vec3(float, float)> df_u, vec2 t_range, vec2 u_range,
                                                 bool t_periodic, bool u_periodic) {
    _f = f;
    _df_t = df_t;
    _df_u = df_u;
    t0 = t_range.x;
    t1 = t_range.y;
    u0 = u_range.x;
    u1 = u_range.y;
    this->t_periodic = t_periodic;
    this->u_periodic = u_periodic;
    epsilon = .01;

}
SmoothParametricSurface::SmoothParametricSurface(std::function<vec3(float, float)> f, vec2 t_range, vec2 u_range,
                                                 bool t_periodic, bool u_periodic, float epsilon)
        : SmoothParametricSurface(f,
           [g=f, epsilon](float t, float u) {return (g(t+epsilon, u) - g(t-epsilon, u))/(epsilon*2); },
           [g=f, epsilon](float t, float u) {return (g(t, u+epsilon) - g(t, u-epsilon))/(epsilon*2); },
           t_range, u_range, t_periodic, u_periodic) {
    this->epsilon = epsilon;
}

SmoothParametricSurface::SmoothParametricSurface(std::function<SmoothParametricCurve(float)> pencil, vec2 t_range, vec2 u_range,
                                                 bool t_periodic, bool u_periodic, float eps) :
    SmoothParametricSurface([pencil](float t, float u) { return pencil(t)(u); }, t_range, u_range, t_periodic, u_periodic, eps) {}

SmoothParametricSurface::SmoothParametricSurface(SmoothParametricCurve canal, std::function<float(float)> width, float eps)
    : SmoothParametricSurface([canal, width](float t, float theta) {
    return canal(t) + (canal.normal(t)*cos(theta) + canal.binormal(t)*sin(theta))*width(t); },
         canal.bounds(), vec2(0, TAU), canal.isPeriodic(), true, eps) {}


SmoothParametricCurve SmoothParametricSurface::precompose(SmoothParametricPlaneCurve c, PolyGroupID id) const {
    return SmoothParametricCurve([f=*this, g=c](float t) {return f(g(t)); }, id, c.bounds().x, c.bounds().y, c.isPeriodic(), epsilon);
}

SmoothParametricCurve SmoothParametricSurface::restrictToInterval(vec2 p0, vec2 p1, PolyGroupID id) const {
    return SmoothParametricCurve([f=*this, p0, p1](float t) {return f(p1*t + p0*(1-t)); }, id, 0, 1, false, epsilon);
}


vec3 SmoothParametricSurface::operator()(float t, float s) const {
	return _f(t, s);
}

vec3 SmoothParametricSurface::operator()(vec2 tu) const { return _f(tu.x, tu.y); }

vec3 SmoothParametricSurface::normal(float t, float u) const {
    vec3 vn = cross(_df_t(t, u), _df_u(t, u));
    if (norm(vn) < .02) {
        float e_t = (t1 - t0) / 40;
        float e_u = (u1 - u0) / 40;
        vec3 new_n = (cross(_df_t(t + e_t, u), _df_u(t + e_t, u)) + cross(_df_t(t, u + e_u), _df_u(t, u + e_u)) +
                      cross(_df_t(t - e_t, u), _df_u(t - e_t, u)) + cross(_df_t(t, u - e_u), _df_u(t, u - e_u))) /
                     4.f;

        if (norm(new_n) < .03) {
            e_t = (t1 - t0) / 20;
            e_u = (u1 - u0) / 20;
            new_n = (cross(_df_t(t + e_t, u), _df_u(t + e_t, u)) + cross(_df_t(t, u + e_u), _df_u(t, u + e_u)) +
                     cross(_df_t(t - e_t, u), _df_u(t - e_t, u)) + cross(_df_t(t, u - e_u), _df_u(t, u - e_u))) /
                    4.f;
            return norm(new_n) < .001 ? vec3(0, 0, 1) : normalise(new_n);
        }

        return normalise(new_n);
    }
    return normalize(vn);
}
void SmoothParametricSurface::normaliseDomainToI2() {
    _f = [f=_f, t0=t0, t1=t1, u0=u0, u1=u1](float t, float u) {return f(lerp(t0, t1, t), lerp(u0, u1, u)); };
    _df_t = [df=_df_t, t0=t0, t1=t1, u0=u0, u1=u1](float t, float u) {return (t1-t0)*df(lerp(t0, t1, t), lerp(u0, u1, u)); };
    _df_u = [df=_df_u, t0=t0, t1=t1, u0=u0, u1=u1](float t, float u) {return (u1-u0)*df(lerp(t0, t1, t), lerp(u0, u1, u)); };
    t0 = 0;
    t1 = 1;
    u0 = 0;
    u1 = 1;
    epsilon /= max(t1-t0, u1-u0);
}

SmoothParametricCurve SmoothParametricPlaneCurve::embedding(vec3 v1, vec3 v2, vec3 pivot) const {
    auto aff = mat3(v1, v2, cross(v1, v2));
    return SmoothParametricCurve([pivot, aff, f=this->_f](float t) {return aff*vec3(f(t), 0) + pivot; },
                                           [aff, d=this->_df](float t) {return aff*vec3(d(t), 0); },
                                           [aff, dd=this->_ddf](float t) {return aff*vec3(dd(t), 0); });
}
