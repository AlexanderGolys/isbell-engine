#include "smoothParametric.hpp"
#include <map>
#include <glm/gtc/constants.hpp>

using namespace glm;
using std::vector, std::string, std::map, std::shared_ptr, std::unique_ptr, std::pair, std::make_unique, std::make_shared, std::function;





SmoothParametricCurve SmoothParametricSurface::precompose(SmoothParametricPlaneCurve c, PolyGroupID id) const {
    return SmoothParametricCurve([f=*this, g=c](float t) {return f(g(t)); }, id, c.bounds().x, c.bounds().y, c.isPeriodic(), epsilon);
}

SmoothParametricCurve SmoothParametricSurface::restrictToInterval(vec2 p0, vec2 p1, PolyGroupID id) const {
    return SmoothParametricCurve([f=*this, p0, p1](float t) {return f(p1*t + p0*(1-t)); }, id, 0, 1, false, epsilon);
}




SmoothParametricSurface AffineLine::tube(float radius, float t0, float t1) const {
	auto w = orthogonalComplementBasis(direction());
	SmoothParametricCurve circle = SmoothParametricCurve([radius, w, p=(*this)(t0)](float t) {return  p + radius*cos(t)*w.first + radius*sin(t)*w.second; });
	return circle.cylinder(direction(), t1 - t0);
}

SmoothParametricSurface::SmoothParametricSurface(function<vec3(float, float)> f, glm::vec2 t_range, glm::vec2 u_range, bool t_periodic, bool u_periodic, float epsilon) {
	_f = f;
	t0 = t_range.x;
	t1 = t_range.y;
	u0 = u_range.x;
	u1 = u_range.y;
	this->t_periodic = t_periodic;
	this->u_periodic = u_periodic;
	this->epsilon = epsilon;
	_df_t = [f, epsilon](float t, float u) {return (f(t + epsilon, u) - f(t - epsilon, u)) / (2 * epsilon); };
	_df_u = [f, epsilon](float t, float u) {return (f(t, u + epsilon) - f(t, u - epsilon)) / (2 * epsilon); };
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

SmoothParametricSurface SmoothParametricCurve::surfaceOfRevolution(const AffineLine& axis) const
{
	Foo113 param = [f = _f, axis](float t, float s) {
		vec3 p0 = axis.orthogonalProjection(f(t));
		vec3 v1 = f(t) - p0;
		vec3 v2 = cross(axis.direction(), v1);
		return p0 + v1*cos(s) + v2*sin(s);
	};
	return SmoothParametricSurface(param, vec2(t0, t1), vec2(0, TAU), periodic, true, eps);
}

SmoothParametricSurface SmoothParametricCurve::screwMotion(float speed, int iterations) const {
	Foo113 param = [f = _f, speed](float t, float s) {
		return vec3(f(t).x*cos(s) - f(t).y*sin(s), f(t).x*sin(s) + f(t).y*cos(s), f(t).z + speed*s);
	};
	return SmoothParametricSurface(param, vec2(t0, t1), vec2(0, TAU*iterations), periodic, false, eps);
}

SmoothParametricSurface SmoothParametricCurve::cylinder(vec3 direction, float h) const {
	Foo113 param = [f = _f, direction](float t, float s) {
		return f(t) + normalise(direction)*s;
	};
	return SmoothParametricSurface(param, vec2(t0, t1), vec2(0, h), periodic, false, eps);
}

SmoothParametricSurface SmoothParametricCurve::pipe(float r) const {
	Foo113 param = [f=_f, df=_df, eps=this->eps, r](float t, float s) {
		vec3 tangent = normalise(df(t));
		if (norm(tangent) < 1e-6)
			tangent = normalise(df(t+eps));
		auto v = orthogonalComplementBasis(tangent);
		return f(t) + r*cos(s)*v.first + r*sin(s)*v.second;
	};
	return SmoothParametricSurface(param, vec2(t0, t1), vec2(0, TAU), periodic, true, eps);
}
SmoothParametricSurface SmoothParametricCurve::canal(function<float(float)> radius) const {
	Foo113 param = [c=*this, r=radius, dr=derivativeOperator(radius, eps)](float t, float s) {

		return c(t) - r(t)*dr(t)/norm2(c.df(t))*c.df(t) + (c.normal(t)*cos(s) + c.binormal(t)*sin(s))*r(t)*sqrt(1 - dr(t)*dr(t)/norm2(c.df(t)));
	};
	return SmoothParametricSurface(param, vec2(t0, t1), vec2(0, TAU), periodic, true, eps);
}
