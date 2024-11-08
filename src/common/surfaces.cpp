#include "surfaces.hpp"


using namespace glm;


SmoothParametricSurface::SmoothParametricSurface(std::function<glm::vec3(float, float)> f, std::function<glm::vec3(float, float)> df_t,
                                                 std::function<glm::vec3(float, float)> df_u, glm::vec2 t_range, glm::vec2 u_range,
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
SmoothParametricSurface::SmoothParametricSurface(std::function<glm::vec3(float, float)> f, glm::vec2 t_range, glm::vec2 u_range,
                                                 bool t_periodic, bool u_periodic, float epsilon)
        : SmoothParametricSurface(f,
           [g=f, epsilon](float t, float u) {return (g(t+epsilon, u) - g(t-epsilon, u))/(epsilon*2); },
           [g=f, epsilon](float t, float u) {return (g(t, u+epsilon) - g(t, u-epsilon))/(epsilon*2); },
           t_range, u_range, t_periodic, u_periodic) {
    this->epsilon = epsilon;
}

SmoothParametricSurface::SmoothParametricSurface(std::function<SmoothParametricCurve(float)> pencil, glm::vec2 t_range, glm::vec2 u_range,
                                                 bool t_periodic, bool u_periodic, float eps) :
    SmoothParametricSurface([pencil](float t, float u) { return pencil(t)(u); }, t_range, u_range, t_periodic, u_periodic, eps) {}

SmoothParametricSurface::SmoothParametricSurface(SmoothParametricCurve canal, std::function<float(float)> width, float eps)
    : SmoothParametricSurface([canal, width](float t, float theta) {
    return canal(t) + (canal.normal(t)*cos(theta) + canal.binormal(t)*sin(theta))*width(t); },
         canal.bounds(), vec2(0, TAU), canal.isPeriodic(), true, eps) {}


SmoothParametricCurve SmoothParametricSurface::precompose(SmoothParametricPlaneCurve c) const {
    return SmoothParametricCurve([f=*this, g=c](float t) {return f(g(t)); }, c.bounds().x, c.bounds().y, c.isPeriodic(), epsilon);
}

SmoothParametricCurve SmoothParametricSurface::restrictToInterval(glm::vec2 p0, glm::vec2 p1) const {
    return SmoothParametricCurve([f=*this, p0, p1](float t) {return f(p1*t + p0*(1-t)); }, 0, 1, false, epsilon);
}


vec3 SmoothParametricSurface::operator()(float t, float s) const {
	return _f(t, s);
}

vec3 SmoothParametricSurface::operator()(vec2 tu) const { return _f(tu.x, tu.y); }

glm::vec3 SmoothParametricSurface::normal(float t, float u) const {
vec3 vn = cross(_df_t(t, u), _df_u(t, u));
    if (norm(vn) < .02) {
        float e_t = (t1 - t0)/40;
        float e_u = (u1 - u0)/40;
        vec3 new_n = (cross(_df_t(t+e_t, u), _df_u(t+e_t, u)) +
                      cross(_df_t(t, u+e_u), _df_u(t, u+e_u)) +
                      cross(_df_t(t-e_t, u), _df_u(t-e_t, u)) +
                      cross(_df_t(t, u-e_u), _df_u(t, u-e_u)))/4.f;

        if(norm(new_n) < .03) {
            e_t = (t1 - t0)/20;
            e_u = (u1 - u0)/20;
            new_n = (cross(_df_t(t+e_t, u), _df_u(t+e_t, u)) +
                     cross(_df_t(t, u+e_u), _df_u(t, u+e_u)) +
                     cross(_df_t(t-e_t, u), _df_u(t-e_t, u)) +
                     cross(_df_t(t, u-e_u), _df_u(t, u-e_u)))/4.f;
            return norm(new_n) < .001 ? vec3(0, 0, 1) : normalise(new_n);  }

        return normalise(new_n); }
return normalize(vn); }



SmoothImplicitSurface::SmoothImplicitSurface(SmoothRealFunctionR3 F) {
	_F = F;
}

float SmoothImplicitSurface::operator()(vec3 p) const {
	return _F(p);
}

AffinePlane::AffinePlane(vec3 n, float d) :
SmoothImplicitSurface(SmoothRealFunctionR3(
	[n, d](vec3 p) {return dot(p, n) - d; }, [n](vec3 p) {return n; })) {
	this->n = normalise(n);
	this->d = d;
	auto tangents = orthogonalComplementBasis(n);
	this->v1 = tangents.first;
	this->v2 = tangents.second;
	this->pivot = this->n * d;
}

AffinePlane::AffinePlane(vec3 pivot, vec3 v1, vec3 v2):
SmoothImplicitSurface(SmoothRealFunctionR3(
	[v1, v2, pivot](vec3 p) {return dot(p, normalise(cross(v1, v2))) - dot(pivot, normalise(cross(v1, v2))); },
	[v1, v2](vec3 p) {return normalise(cross(v1, v2)); })) {
	this->pivot = pivot;
	this->v1 = v1;
	this->v2 = v2;
	this->n = normalise(cross(v1, v2));
	this->d = dot(pivot, n);
}


AffinePlane AffinePlane::spanOfPts(vec3 p0, vec3 p1, vec3 p2) {
	return AffinePlane(cross(p1-p0, p2-p0), dot(p0, cross(p1-p0, p2-p0)));
}

AffineLine AffinePlane::intersection(AffinePlane &p) const {
	vec3 v = cross(n, p.normal());
	float denominator = norm2(n)*norm2(p.normal()) - dot(n, p.normal())*dot(n, p.normal());
	vec3 pivot = n*(d*norm2(p.normal()) - p.getD()*dot(n, p.normal()))/denominator +
				p.normal()*(p.getD()*norm2(n) - d*dot(n, p.normal()))/denominator;
	return AffineLine(pivot, v);
}


vec3 AffinePlane::intersection(AffineLine &l) const {
	return l.pivot() + dot(pivot - l.pivot(), n) / dot(l.direction(), n) * l.direction();
}


AffinePlane::operator SmoothImplicitSurface() const {
	return SmoothImplicitSurface(SmoothRealFunctionR3(
		[this](vec3 p) {return distance(p); },
		[this](vec3 p) {return n; }));
}

mat3 AffinePlane::pivotAndBasis() const {
	return mat3(pivot, v1, v2);
}

vec2 AffinePlane::localCoordinates(vec3 p) const {
	return vec2(dot(p - pivot, v1), dot(p - pivot, v2));
}

AffinePlane SmoothParametricCurve::osculatingPlane(float t) const {
    return AffinePlane(binormal(t), dot(binormal(t), _f(t)));
}