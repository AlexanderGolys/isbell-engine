#include "func.hpp"
#include <functional>
#include <glm/detail/_vectorize.hpp>
#include <glm/gtx/transform.hpp>
#include <iosfwd>
#include <iostream>
#include <memory>
#include <vector>

using namespace glm;
using std::vector, std::shared_ptr, std::make_shared;




VectorFieldR2::VectorFieldR2() {
	this->field = [](vec2 v) {return vec2(NAN, NAN); };
}

VectorFieldR2::VectorFieldR2(std::function<vec2(vec2)> field)
{
	this->field = field;

}

vec2 VectorFieldR2::operator()(vec2 v) const {
	return this->field(v);
}


std::function<float(float)> derivativeOperator(std::function<float(float)> f, float epsilon) {
	return [f, epsilon](float x) {return (f(x + epsilon) - f(x - epsilon)) / (2 * epsilon); };
}

std::function<vec2(float)> derivativeOperator(std::function<vec2(float)> f, float epsilon) {
	auto fx = [f](float x) {return f(x).x; };
	auto fy = [f](float x) {return f(x).y; };
	return [fx, fy, epsilon](float x) {return vec2(derivativeOperator(fx, epsilon)(x), derivativeOperator(fy, epsilon)(x)); };
}
std::function<vec3(float)> derivativeOperator(std::function<vec3(float)> f, float epsilon) {
	auto fx = [f](float x) {return f(x).x; };
	auto fy = [f](float x) {return f(x).y; };
	auto fz = [f](float x) {return f(x).z; };
	return [fx, fy, fz, epsilon](float x) {return vec3(derivativeOperator(fx, epsilon)(x), derivativeOperator(fy, epsilon)(x), derivativeOperator(fz, epsilon)(x)); };
}




SmoothRealFunctionR3::SmoothRealFunctionR3() {
	this->_f = [](vec3 v) {return 0; };
	this->_df = [](vec3 v) {return vec3(0, 0, 0); };
}


SmoothRealFunctionR3::SmoothRealFunctionR3(std::function<float(vec3)> f, std::function<vec3(vec3)> df) {
	this->_f = f;
	this->_df = df;
}

SmoothRealFunctionR3::SmoothRealFunctionR3(std::function<float(vec3)> f, float epsilon) {
	this->_f = f;
	this->_df = [f, epsilon](vec3 v) {
		return  (f(v + vec3(epsilon, 0, 0)) - f(v - vec3(epsilon, 0, 0))) / (2 * epsilon) * vec3(1, 0, 0) +
			(f(v + vec3(0, epsilon, 0)) - f(v - vec3(0, epsilon, 0))) / (2 * epsilon) * vec3(0, 1, 0) +
			(f(v + vec3(0, 0, epsilon)) - f(v - vec3(0, 0, epsilon))) / (2 * epsilon) * vec3(0, 0, 1);
	};
}

float SmoothRealFunctionR3::operator()(vec3 v) const {
	return _f(v);
}

vec3 SmoothRealFunctionR3::df(vec3 v) const {
	return _df(v);
}

SmoothRealFunctionR3 SmoothRealFunctionR3::operator+(SmoothRealFunctionR3 g) const {
	return SmoothRealFunctionR3([this, g](vec3 v) {return this->_f(v) + g(v); },
		[this, g](vec3 v) {return this->_df(v) + g.df(v); });
}

SmoothRealFunctionR3 SmoothRealFunctionR3::operator-(SmoothRealFunctionR3 g) const {
	return SmoothRealFunctionR3([this, g](vec3 v) {return this->_f(v) - g(v); },
		[this, g](vec3 v) {return this->_df(v) - g.df(v); });
}

SmoothRealFunctionR3 SmoothRealFunctionR3::operator*(SmoothRealFunctionR3 g) const {
	return SmoothRealFunctionR3([this, g](vec3 v) {return this->_f(v) * g(v); },
		[this, g](vec3 v) {return this->_f(v) * g.df(v) + this->_df(v) * g(v); });
}

SmoothRealFunctionR3 SmoothRealFunctionR3::operator/(SmoothRealFunctionR3 g) const {
	return SmoothRealFunctionR3([this, g](vec3 v) {return this->_f(v) / g(v); },
		[this, g](vec3 v) {return (this->_df(v) * g(v) - this->_f(v) * g.df(v)) / (g(v) * g(v)); });
}

SmoothRealFunctionR3 SmoothRealFunctionR3::operator*(float a) const {
	return SmoothRealFunctionR3([this, a](vec3 v) {return this->_f(v) * a; },
		[this, a](vec3 v) {return this->_df(v) * a; });
}

SmoothRealFunctionR3 SmoothRealFunctionR3::operator+(float a) const {
	return SmoothRealFunctionR3([this, a](vec3 v) {return this->_f(v) + a; },
		[this](vec3 v) {return this->_df(v); });
}

SmoothRealFunctionR3 SmoothRealFunctionR3::operator-(float a) const {
	return SmoothRealFunctionR3([this, a](vec3 v) {return this->_f(v) - a; },
		[this](vec3 v) {return this->_df(v); });
}

SmoothRealFunctionR3 SmoothRealFunctionR3::operator/(float a) const {
	return SmoothRealFunctionR3([this, a](vec3 v) {return this->_f(v) / a; },
		[this, a](vec3 v) {return this->_df(v) / a; });
}

VectorFieldR3 SmoothRealFunctionR3::gradient() const {
	return VectorFieldR3([this](vec3 v) {return this->_df(v); });
}

SmoothRealFunctionR3 SmoothRealFunctionR3::linear(vec3 v) {
	return SmoothRealFunctionR3([v](vec3 x) {return dot(x, v); },
		[v](vec3 x) {return v; });
}

SmoothRealFunctionR3 SmoothRealFunctionR3::projection(int i) {
    return SmoothRealFunctionR3([i](vec3 x) { return x[i]; }, [i](vec3 x) { return vec3(0, 0, 0); });
}
SpaceEndomorphism &SpaceEndomorphism::operator=(const SpaceEndomorphism &other) {
    if (this == &other)
        return *this;
    _f = other._f;
    _df = other._df;
    return *this;
}
SpaceEndomorphism &SpaceEndomorphism::operator=(SpaceEndomorphism &&other) noexcept {
    if (this == &other)
        return *this;
    _f = std::move(other._f);
    _df = std::move(other._df);
    return *this;
}


SmoothRealFunctionR3 SmoothRealFunctionR3::constant(float a) {
	return SmoothRealFunctionR3([a](vec3 x) {return a; },
		[a](vec3 x) {return vec3(0, 0, 0); });
}


vec3 SmoothParametricCurve::operator()(float t) const { return this->_f(t); }


SmoothParametricCurve::SmoothParametricCurve(std::function<vec3(float)> f, std::function<vec3(float)> df,
                                             std::function<vec3(float)> ddf, std::variant<int, std::string> id,
                                             std::optional<float> t0, std::optional<float> t1, bool periodic, float epsilon) : SmoothParametricCurve(f, df, ddf, t0, t1, periodic, epsilon) {
    this->id = id;
}

SmoothParametricCurve::SmoothParametricCurve(std::function<vec3(float)> f, std::function<vec3(float)> df,
                                             std::function<vec3(float)> ddf, RP1 t0, RP1 t1, bool periodic, float epsilon) {
    this->_f = f;
    this->_df = df;
    this->_ddf = ddf;
    eps = epsilon;
    this->t0 = t0;
    this->t1 = t1;
    this->periodic = periodic;
    id = randomCurvaID();
}
SmoothParametricCurve::SmoothParametricCurve(std::function<vec3(float)> f, std::function<vec3(float)> df, std::optional<float> t0,
                                             std::optional<float> t1, bool periodic, float epsilon) : SmoothParametricCurve(f, df, derivativeOperator(df, epsilon), t0, t1, periodic, epsilon) {}

SmoothParametricCurve::SmoothParametricCurve(std::function<vec3(float)> f, std::function<vec3(float)> df,
                                             std::variant<int, std::string> id, std::optional<float> t0, std::optional<float> t1,
                                             bool periodic, float epsilon) : SmoothParametricCurve(f, df, derivativeOperator(df, epsilon), id, t0, t1, periodic, epsilon) {}


SmoothParametricCurve::SmoothParametricCurve(const SmoothParametricCurve &other) :
    _f(other._f), _df(other._df), _ddf(other._ddf), _der_higher(other._der_higher), eps(other.eps), t0(other.t0), t1(other.t1),
    periodic(other.periodic), id(other.id){}
SmoothParametricCurve::SmoothParametricCurve(SmoothParametricCurve &&other) noexcept :
    _f(std::move(other._f)), _df(std::move(other._df)), _ddf(std::move(other._ddf)), _der_higher(std::move(other._der_higher)),
    eps(other.eps), t0(std::move(other.t0)), t1(std::move(other.t1)), periodic(other.periodic), id(other.id) {}
SmoothParametricCurve &SmoothParametricCurve::operator=(const SmoothParametricCurve &other) {
    if (this == &other)
        return *this;
    _f = other._f;
    _df = other._df;
    _ddf = other._ddf;
    _der_higher = other._der_higher;
    eps = other.eps;
    t0 = other.t0;
    t1 = other.t1;
    periodic = other.periodic;
    id = other.id;
    return *this;
}
SmoothParametricCurve &SmoothParametricCurve::operator=(SmoothParametricCurve &&other) noexcept {
    if (this == &other)
        return *this;
    _f = std::move(other._f);
    _df = std::move(other._df);
    _ddf = std::move(other._ddf);
    _der_higher = std::move(other._der_higher);
    eps = other.eps;
    t0 = std::move(other.t0);
    t1 = std::move(other.t1);
    periodic = other.periodic;
    id = other.id;
    return *this;
}
SmoothParametricCurve::SmoothParametricCurve(std::function<vec3(float)> f, std::function<std::function<vec3(float)>(int)> derivativeOperator,
                                             RP1 t0, RP1 t1, bool periodic, float epsilon) {
	this->_f = f;
	this->_df = derivativeOperator(1);
	this->_ddf = derivativeOperator(2);
	this->eps = epsilon;
	this->_der_higher = derivativeOperator;
    this->t0 = t0;
    this->t1 = t1;
    this->periodic = periodic;
    id = randomCurvaID();
}

SmoothParametricCurve::SmoothParametricCurve(std::function<vec3(float)> f, std::vector<std::function<vec3(float)>> derivatives,
                                             RP1 t0, RP1 t1, bool periodic, float epsilon) {
	this->_f = f;
	if (derivatives.size() > 0)
		this->_df = derivatives[0];
	else
		this->_df = derivativeOperator(f, epsilon);

	if (derivatives.size() > 1)
		this->_ddf = derivatives[1];
	else
		this->_ddf = derivativeOperator(_df, epsilon);

	this->_der_higher = [derivatives, f](int n) { return n==0 ? f : n <= derivatives.size() ? derivatives[n-1] : derivativeOperator(derivatives[n-1], 0.01f);} ;
	this->eps = epsilon;
    this->t0 = t0;
    this->t1 = t1;
    this->periodic = periodic;
    id = randomCurvaID();
}

AffinePlane SmoothParametricCurve::osculatingPlane(float t) const {
	return AffinePlane(binormal(t), dot(binormal(t), _f(t)));
}

SmoothParametricCurve SmoothParametricCurve::constCurve(vec3 v) {
	return SmoothParametricCurve([v](float t) {return v; },inf, inf, false);
}

float SmoothParametricCurve::length(float t0, float t1, int n) const {
    float res = 0;
    float dt = (t1 - t0) / n;
    for (int i = 0; i < n; i++) {
        res += norm(_f(lerp(t0, t1, dt * i)) - _f(lerp(t0, t1, dt * (i + 1))));
    }
    return res;
}
SmoothParametricCurve SmoothParametricCurve::precompose(SpaceEndomorphism g_) const {
    return SmoothParametricCurve([f = this->_f, g=g_](float t) {return g(f(t)); },
                                  [f = this->_f, d = this->_df, g=g_](float t) {return g.df(f(t)) * d(t); },
                                  id, this->t0, this->t1, this->periodic, this->eps);
}
void SmoothParametricCurve::precomposeInPlace(SpaceEndomorphism g) {
    this->_f = [f = this->_f, g](float t) {return g(f(t)); };
    this->_df = [f = this->_f, d = this->_df, g](float t) {return g(f(t)) * d(t); };
    this->_ddf = [f = this->_f, d = this->_df, dd = this->_ddf, g](float t) {return g(f(t)) * dd(t) + g.df(f(t)) * d(t); };
}


mat3 SmoothParametricCurve::FrenetFrame(float t) const {
	return mat3(tangent(t), normal(t), binormal(t));

}

float SmoothParametricCurve::curvature(float t) const {
	return norm(cross(_df(t), _ddf(t))) / pow(norm(_df(t)), 3);
}

float SmoothParametricCurve::torsion(float t) const {
	return dot(cross(df(t), ddf(t)), higher_derivative(t, 3)) / norm2(cross(df(t), ddf(t)));
}

vec3 SmoothParametricCurve::curvature_vector(float t) const {
	return normal(t)*curvature(t);
}

AffineLine::AffineLine(vec3 p0, vec3 v) :
SmoothParametricCurve([p0, v](float t) {return p0 + t * v; },
					  [v](float t) {return v; },
					  [](float t) {return vec3(0, 0, 0); }) {
	this->p0 = p0;
	this->v = v;
}

AffineLine AffineLine::spanOfPts(vec3 p0, vec3 p1) {
	return AffineLine(p0, p1 - p0);
}


SmoothRealFunctionR3 AffineLine::distanceField() const {
	return SmoothRealFunctionR3([this](vec3 p) {return distance(p); });
}

vec3 AffineLine::orthogonalProjection(vec3 p) const {
	return p0 + dot(p - p0, v) * v / dot(v, v);
}

vec3 AffineLine::pivot() const {
	return p0;
}

vec3 AffineLine::direction() const {
	return v;
}

float AffineLine::distance(AffineLine &l) const {
	return dot((p0 - l.pivot()), cross(v, l.direction())) / norm(cross(v, l.direction()));
}

AffineLine AffineLine::operator+(vec3 v) const {
	return AffineLine(p0 + v, this->v);
}

SmoothParametricSurface::SmoothParametricSurface(std::function<vec3(float, float)> f,
	std::function<vec3(float, float)> df_t, std::function<vec3(float, float)> df_u, float t0, float t1,
	float u0, float u1, bool t_periodic, bool u_periodic) {
	_f = f;
	_df_t = df_t;
	_df_u = df_u;
	this->t0 = t0;
	this->t1 = t1;
	this->u0 = u0;
	this->u1 = u1;
	this->t_periodic = t_periodic;
	this->u_periodic = u_periodic;
}

SmoothParametricSurface::SmoothParametricSurface(std::function<vec3(float, float)> f, float t0, float t1, float u0,
	float u1, bool t_periodic, bool u_periodic, float epsilon) {
	_f = f;
	_df_t = [f, epsilon](float t, float u) {
		return (f(t + epsilon, u) - f(t - epsilon, u)) / (2 * epsilon);
	};
	_df_u = [f, epsilon](float t, float u) {
		return (f(t, u + epsilon) - f(t, u - epsilon)) / (2 * epsilon);
	};
	this->t0 = t0;
	this->t1 = t1;
	this->u0 = u0;
	this->u1 = u1;
	this->t_periodic = t_periodic;
	this->u_periodic = u_periodic;
}

SmoothParametricSurface::SmoothParametricSurface(std::function<SmoothParametricCurve(float)> pencil, float t0, float t1,
	float curve_t0, float curve_t1, bool periodic, bool curve_periodic, float eps) {
	_f = [pencil](float t, float u) {return pencil(u)(t); };
	_df_t = [pencil](float t, float u) {return pencil(u).derivative(t); };
	_df_u = [pencil, eps](float t, float u) {return (pencil(u+eps)(t) - pencil(u - eps)(t))/(2*eps); };
	this->u0 = t0;
	this->u1 = t1;
	this->t0 = curve_t0;
	this->t1 = curve_t1;
	u_periodic = periodic;
	t_periodic = curve_periodic;
}

vec3 SmoothParametricSurface::operator()(float t, float s) const {
	return _f(t, s);
}

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

float AffinePlane::distance(vec3 p) const {
	return (*this)(p);
}

vec3 AffinePlane::orthogonalProjection(vec3 p) const {
	return p - (*this)(p) * n;
}

bool AffinePlane::contains(vec3 p, float eps) const {
	return abs(distance(p)) < eps;
}

vec3 AffinePlane::normal() const {
	return n;
}

float AffinePlane::getD() const {
	return d;
}

std::pair<vec3, float> AffinePlane::equationCoefs() const {
	return std::pair<vec3, float>(n, d);
}

vec3 AffinePlane::intersection(AffineLine &l) const {
	return l.pivot() + dot(pivot - l.pivot(), n) / dot(l.direction(), n) * l.direction();
}

SmoothRealFunctionR3 AffinePlane::distanceField() const {
	return SmoothRealFunctionR3([this](vec3 p) {return distance(p); });
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

float AffineLine::distance(vec3 p) const {
	return norm(cross(p - p0, v)) / norm(v);
}

bool AffineLine::contains(vec3 p, float eps) const {
	return distance(p) < eps;
}

SmoothParametricPlaneCurve::SmoothParametricPlaneCurve(std::function<vec2(float)> f, std::function<vec2(float)> df,
                                                       std::function<vec2(float)> ddf, float t0, float t1, bool period, float epsilon) {
    this->_f = f;
    this->t0 = t0;
    this->t1 = t1;
    this->periodic = period;
    this->eps = epsilon;
    this->_df = df;
    this->_ddf =ddf;
}
SmoothParametricPlaneCurve::SmoothParametricPlaneCurve(const SmoothParametricPlaneCurve &other) :
    _f(other._f), _df(other._df), _ddf(other._ddf), _der_higher(other._der_higher), eps(other.eps), t0(other.t0), t1(other.t1),
    periodic(other.periodic) {}
SmoothParametricPlaneCurve::SmoothParametricPlaneCurve(SmoothParametricPlaneCurve &&other) noexcept :
    _f(std::move(other._f)), _df(std::move(other._df)), _ddf(std::move(other._ddf)), _der_higher(std::move(other._der_higher)),
    eps(other.eps), t0(std::move(other.t0)), t1(std::move(other.t1)), periodic(other.periodic) {}
SmoothParametricPlaneCurve &SmoothParametricPlaneCurve::operator=(const SmoothParametricPlaneCurve &other) {
    if (this == &other)
        return *this;
    _f = other._f;
    _df = other._df;
    _ddf = other._ddf;
    _der_higher = other._der_higher;
    eps = other.eps;
    t0 = other.t0;
    t1 = other.t1;
    periodic = other.periodic;
    return *this;
}
SmoothParametricPlaneCurve &SmoothParametricPlaneCurve::operator=(SmoothParametricPlaneCurve &&other) noexcept {
    if (this == &other)
        return *this;
    _f = std::move(other._f);
    _df = std::move(other._df);
    _ddf = std::move(other._ddf);
    _der_higher = std::move(other._der_higher);
    eps = other.eps;
    t0 = std::move(other.t0);
    t1 = std::move(other.t1);
    periodic = other.periodic;
    return *this;
}

SmoothParametricPlaneCurve::SmoothParametricPlaneCurve( std::function<vec2(float)> f,
                                                        const std::function<vec2(float)>& df, float t0, float t1,
                                                        bool period, float epsilon) :
                            SmoothParametricPlaneCurve(f, df, derivativeOperator(df, epsilon), t0, t1, period, epsilon) {}



SmoothParametricPlaneCurve::SmoothParametricPlaneCurve( const std::function<vec2(float)>& curve, float t0, float t1,
                                                        bool period, float epsilon) :
                            SmoothParametricPlaneCurve(curve, derivativeOperator(curve, epsilon), t0, t1, period, epsilon) {}


vector<vec2> SmoothParametricPlaneCurve::sample(float t0, float t1, int n) const {
	vector<vec2> result = std::vector<vec2>();
	result.reserve(n);
	float dt = (t1 - t0) / (n-1);
	for (int i = 0; i < n; i++)
	{
		result.push_back(_f(t0 + dt*i));
	}
	return result;
}

vector<vec3> SmoothParametricPlaneCurve::adjacency_lines_buffer(float t0,
                                                                float t1, int n,
                                                                float z) const {
  vector<vec3> result = vector<vec3>();
  result.reserve(4 * n - 4);
  float dt = (t1 - t0) / n;
  for (int i = 1; i < n; i++) {
    vec3 p1 = vec3(_f(t0 + (i - 1) * dt), z);
    vec3 p2 = vec3(_f(t0 + i * dt), z);
    vec3 p0 = p1 + vec3(normal(t0 + (i - 1) * dt), 0);
    vec3 p3 = p2 + vec3(normal(t0 + i * dt), 0);
    result.push_back(p0);
    result.push_back(p1);
    result.push_back(p2);
    result.push_back(p3);
  }
  return result;
}
SmoothParametricCurve SmoothParametricPlaneCurve::embedding(vec3 v1, vec3 v2, vec3 pivot) const {
  auto aff = mat3(v1, v2, cross(v1, v2));
  return SmoothParametricCurve([pivot, aff, f=this->_f](float t) {return aff*vec3(f(t), 0) + pivot; },
                                         [aff, d=this->_df](float t) {return aff*vec3(d(t), 0); },
                                         [aff, dd=this->_ddf](float t) {return aff*vec3(dd(t), 0); });
}

ComplexCurve::ComplexCurve(std::function<Complex(float)> curve, float t0, float t1, float period, float epsilon)
{
	this->f = std::make_unique<std::function<Complex(float)>>(curve);
	this->cyclic = period > 0;
	this->epsilon = epsilon;
	this->t0 = t0;
	this->t1 = t1;

	this->df = std::make_unique<std::function<Complex(float)>>([this](float t) {
		return ((*this)(t + this->epsilon) - (*this)(t)) / this->epsilon;
	});

	this->ddf = std::make_unique<std::function<Complex(float)>>([this](float t) {
		return ((*this)(t + this->epsilon) - (*this)(t)*2.f + (*this)(t - this->epsilon)) / (this->epsilon * this->epsilon);
	});

	this->N = std::make_unique<std::function<Complex(float)>>([this](float t) {
		return Complex(normalise<vec2>((*ddf)(t))); });

	this->period = period;
}


Complex ComplexCurve::operator()(float t) const {
    return (*f)(t);
}

std::vector<Complex> ComplexCurve::sample(float t0, float t1, int n)
{
	if (n == 1)
	{
		return { (*f)(t0) };
	}
	std::vector<Complex> result = std::vector<Complex>();
	result.reserve(n);
	float dt = (t1 - t0) / (n-1);
	for (int i = 0; i < n; i++)
	{
		result.push_back((*f)(t0 + i * dt));
	}
	printVector(result, "sampled complex curve from " + std::to_string(t0) + " to " + std::to_string(t1));
	return result;
}

std::vector<Complex> ComplexCurve::sample(int n)
{
	auto result = sample(t0, t1, n);
	return result;
}

ComplexCurve ComplexCurve::disjointUnion(ComplexCurve &other)
{
    return ComplexCurve([this, &other](float t) {
		if (t < this->t1)
		{
			return (*this->f)(t);
		}
		else
		{
			return (*other.f)(t - this->t1+other.t0);
		}
	}, this->t0, this->t1 + other.t1-other.t0, 0, this->epsilon);
}

ComplexCurve ComplexCurve::line(Complex z0, Complex z1)
{
	return ComplexCurve([&z0, &z1](float t) {
		return z0 + (z1 - z0) * t;
		}, 0, 1);
}

ComplexCurve ComplexCurve::circle(Complex z0, float r)
{
	return ComplexCurve([&z0, &r](float t) {
		return z0 + exp(I * t) * r;
		}, 0, TAU, TAU, .001f);
}

ComplexCurve ComplexCurve::arc(Complex center, Complex z0, Complex z1)
{
	float r = abs(z0 - center);
	float phi0 = (z0 - center).arg();
	float phi1 = (z1 - center).arg();
	return ComplexCurve([&center, &r, &phi0, &phi1](float t) {
		return center + exp(I * t) * r;
		}, phi0, phi1);
}

PlaneSmoothEndomorphism::PlaneSmoothEndomorphism() {
	_f = make_shared<std::function<vec2(vec2)>>([](vec2 v) {return v; });
	_df = make_shared<std::function<mat2(vec2)>>([](vec2 v) {return mat2(1, 0, 0, 1); });
}

PlaneSmoothEndomorphism::PlaneSmoothEndomorphism(shared_ptr<std::function<vec2(vec2)>> f, shared_ptr<std::function<mat2(vec2)>> df)
{
	_f = f;
	_df = df;
}

vec2 PlaneSmoothEndomorphism::operator()(vec2 x) const {
	return (*_f)(x);
}

vec2 PlaneSmoothEndomorphism::df(vec2 x, vec2 v) const {
	return (*_df)(x) * v;
}
mat2 PlaneSmoothEndomorphism::df(vec2 x) const {
	return (*_df)(x);
}


PlaneAutomorphism::PlaneAutomorphism(shared_ptr<std::function<vec2(vec2)>> f, shared_ptr<std::function<mat2(vec2)>> df,
	shared_ptr<std::function<vec2(vec2)>> f_inv) : PlaneSmoothEndomorphism(f, df) {
	_f_inv = f_inv;
}

vec2 PlaneAutomorphism::f_inv(vec2 x) const {
	return (*_f_inv)(x);
}

PlaneAutomorphism PlaneAutomorphism::operator~() const {
	auto df_cpy = _df;
	auto f_inv_cpy = _f_inv;
	auto df_inv = [df_cpy, f_inv_cpy](vec2 x) { return inverse((*df_cpy)((*f_inv_cpy)(x))); };
	return PlaneAutomorphism( f_inv_cpy, make_shared<std::function<mat2(vec2)>>(df_inv), _f);
}

Meromorphism::Meromorphism() {
	_f = make_shared<endC>([](Complex z) {return z; });
	_df =  make_shared<endC>([](Complex z) {return ONE; });
}

Meromorphism::Meromorphism(shared_ptr<endC> f, shared_ptr<endC> df) {
	_f = f;
	_df = df;
}


Meromorphism::operator PlaneSmoothEndomorphism() const {
	auto new_f = _f;
	auto new_df = _df;
	return PlaneSmoothEndomorphism(make_shared<std::function<vec2(vec2)>>([new_f](vec2 v) {return vec2((*new_f)(Complex(v))); }),
		make_shared<std::function<mat2(vec2) >>([new_df](vec2 x) {return mat2((*new_df)(Complex(x)).re(), (*new_df)(Complex(x)).im(), -(*new_df)(Complex(x)).im(), (*new_df)(Complex(x)).re()); }));
}

Meromorphism Meromorphism::compose(Meromorphism g) const {
	auto new_f = _f;
	auto new_df = _df;
	auto g_f = g._f;
	auto g_df = g._df;
	return Meromorphism(make_shared<endC>([new_f, g_f](Complex z) {return (*new_f)((*g_f)(z)); }),
		make_shared<endC>([new_f, new_df, g_f, g_df](Complex z) {return (*new_df)((*g_f)(z)) * (*g_df)(z); }));
}

Meromorphism Meromorphism::operator+(Meromorphism g) const {
	auto new_f = _f;
	auto new_df = _df;
	auto g_f = g._f;
	auto g_df = g._df;
	return Meromorphism(make_shared<endC>([new_f, g_f](Complex z) {return (*new_f)(z) + (*g_f)(z); }),
		make_shared<endC>([new_df, g_df](Complex z) {return (*new_df)(z) + (*g_df)(z); }));
}

Meromorphism Meromorphism::operator*(Meromorphism g) const {
	auto new_f = _f;
	auto new_df = _df;
	auto g_f = g._f;
	auto g_df = g._df;
	return Meromorphism(make_shared<endC>([new_f, g_f](Complex z) {return (*new_f)(z) * (*g_f)(z); }),
		make_shared<endC>([new_f, new_df, g_f, g_df](Complex z) {return (*new_df)(z) * (*g_f)(z) + (*new_f)(z) * (*g_df)(z); }));
}

Meromorphism Meromorphism::operator-() const {
	auto new_f = _f;
	auto new_df = _df;
	return Meromorphism(make_shared<endC>([new_f](Complex z) {return -(*new_f)(z); }),
		make_shared<endC>([new_df](Complex z) {return -(*new_df)(z); }));

}

Meromorphism Meromorphism::operator-(Meromorphism g) const {
	auto new_f = _f;
	auto new_df = _df;
	auto g_f = g._f;
	auto g_df = g._df;
	return Meromorphism(make_shared<endC>([new_f, g_f](Complex z) {return (*new_f)(z) - (*g_f)(z); }),
		make_shared<endC>([new_df, g_df](Complex z) {return (*new_df)(z) - (*g_df)(z); }));
}

Meromorphism Meromorphism::operator/(Meromorphism g) const {
	auto new_f = _f;
	auto new_df = _df;
	auto g_f = g._f;
	auto g_df = g._df;
	return Meromorphism(make_shared<endC>([new_f, g_f](Complex z) {return (*new_f)(z) / (*g_f)(z); }),
		make_shared<endC>([new_f, new_df, g_f, g_df](Complex z) {return ((*new_df)(z) * (*g_f)(z) - (*new_f)(z) * (*g_df)(z)) / ((*g_f)(z) * (*g_f)(z)); }));
}


Complex Meromorphism::operator()(Complex z) const {
	return (*_f)(z);
}

Complex Meromorphism::df(Complex z) const {
	return  (*_df)(z);
}


Biholomorphism::Biholomorphism() {
	_f = make_shared<endC>([](Complex z) {return z; });
	_df = make_shared<endC>([](Complex z) {return ONE; });
	_f_inv = make_shared<endC>([](Complex z) {return z; });
}

Biholomorphism::Biholomorphism(shared_ptr<endC> f, shared_ptr<endC> df, shared_ptr<endC> f_inv) : Meromorphism(f, df) {
	_f_inv = f_inv;
}

Biholomorphism Biholomorphism::mobius(Matrix<Complex, 2> mobius) {
	auto _f = make_shared<endC>([mobius](Complex z) {return mobius.mobius(z); });
	auto _df = make_shared<endC>([mobius](Complex z) {return mobius.mobius_derivative(z); });
	auto _f_inv = make_shared<endC>([mobius](Complex z) {return (~mobius).mobius(z); });
	return Biholomorphism(_f, _df, _f_inv);
}

Biholomorphism Biholomorphism::linear(Complex a, Complex b) {
	return mobius(Matrix<Complex, 2>(a, b, 0, 1));
}

Biholomorphism Biholomorphism::_LOG() {
	auto _f = make_shared<endC>([](Complex z) {return log(z); });
	auto _df = make_shared<endC>([](Complex z) {return ONE / z; });
	auto _f_inv = make_shared<endC>([](Complex z) {return exp(z); });
	return Biholomorphism(_f, _df, _f_inv);
}

Biholomorphism Biholomorphism::_EXP() {
	auto _f = make_shared<endC>([](Complex z) {return exp(z); });
	auto _df = make_shared<endC>([](Complex z) {return exp(z); });
	auto _f_inv = make_shared<endC>([](Complex z) {return log(z); });
	return Biholomorphism(_f, _df, _f_inv);
}

Biholomorphism Biholomorphism::power(float a) {
	auto _f = make_shared<endC>([a](Complex z) {return z.pow(a); });
	auto _df = make_shared<endC>([a](Complex z) {return a!=0 ? z.pow(a - 1)*a : ZERO; });
	auto _f_inv = make_shared<endC>([a](Complex z) {return z.pow(1 / a); });
	return Biholomorphism(_f, _df, _f_inv);
}



SpaceEndomorphism::SpaceEndomorphism(std::function<vec3(vec3)> f, std::function<mat3(vec3)> df) {
	_f = f;
	_df = df;
}

SpaceEndomorphism::SpaceEndomorphism(std::function<vec3(vec3)> f, float epsilon) {
	_f = f;
    _df = [f, epsilon](vec3 x) {
        return mat3(
            (f(x + vec3(epsilon, 0, 0)) - f(x - vec3(epsilon, 0, 0))) / (2 * epsilon),
            (f(x + vec3(0, epsilon, 0)) - f(x - vec3(0, epsilon, 0))) / (2 * epsilon),
            (f(x + vec3(0, 0, epsilon)) - f(x - vec3(0, 0, epsilon))) / (2 * epsilon)
        );
    };
}

vec3 SpaceEndomorphism::directional_derivative(vec3 x, vec3 v) const {
	return _df(x) * v;
}

vec3 SpaceEndomorphism::operator()(vec3 v) const {
	return _f(v);
}

SpaceEndomorphism SpaceEndomorphism::compose(SpaceEndomorphism g) const {
	return SpaceEndomorphism([f=_f, g](vec3 v) {return f(g(v)); },
			[d=_df, g](vec3 x) {return d(g(x)) * g.df(x); });
}

SpaceEndomorphism SpaceEndomorphism::linear(mat3 A) {
	return SpaceEndomorphism([A](vec3 v) {return A * v; },
		[A](vec3 x) {return A; });
}

SpaceEndomorphism SpaceEndomorphism::translation(vec3 v) {
	return affine(mat3(1), v);
}

SpaceEndomorphism SpaceEndomorphism::scaling(float x, float y, float z) {
	return linear(mat3(x, 0, 0, 0, y, 0, 0, 0, z));
}

SpaceEndomorphism SpaceEndomorphism::affine(mat3 A, vec3 v) {
	return SpaceEndomorphism([A, v](vec3 x) {return A * x + v; },
		[A](vec3 x) {return A; });
}

SpaceAutomorphism SpaceAutomorphism::linear(mat3 A) {
    return SpaceAutomorphism([A](vec3 v) {return A * v; },
                            [A](vec3 v) {return inverse(A) * v; },
                             [A](vec3 x) {return A; });
}
SpaceAutomorphism SpaceAutomorphism::translation(vec3 v) {
    return SpaceAutomorphism([v](vec3 x) {return x + v; },
                            [v](vec3 x) {return x - v; },
                             [](vec3 x) {return mat3(1); });
}
SpaceAutomorphism SpaceAutomorphism::scaling(float x, float y, float z) {
    return SpaceAutomorphism([x, y, z](vec3 v) {return vec3(v.x*x, v.y*y, v.z*z); },
                            [x, y, z](vec3 v) {return vec3(v.x/x, v.y/y, v.z/z); },
                             [x, y, z](vec3 v) {return mat3(x, 0, 0, 0, y, 0, 0, 0, z); });
}


SpaceAutomorphism SpaceAutomorphism::affine(mat3 A, vec3 v) {
    return SpaceAutomorphism([A, v](vec3 x) { return A * x + v; }, [A, v](vec3 x) { return inverse(A) * (x - v); },
                             [A](vec3 x) { return A; });
}
SpaceAutomorphism SpaceAutomorphism::rotation(float angle) {
    return linear(rotationMatrix3(angle));
}


SpaceAutomorphism SpaceAutomorphism::rotation(vec3 axis, float angle) {
    return linear(rotationMatrix3(axis, angle));
}


VectorFieldR3::VectorFieldR3() {
	_field = [](vec3 v) {return vec3(0, 0, 0); };
}

VectorFieldR3::VectorFieldR3(std::function<vec3(vec3)> field) {
	_field = field;
}

VectorFieldR3::VectorFieldR3(SpaceEndomorphism f) {
	_field = f;
}

VectorFieldR3::VectorFieldR3(VectorFieldR2 field) {
	_field = [X=field](vec3 v) {return vec3(X(vec2(v.x, v.y)), 0); };
}

VectorFieldR3 VectorFieldR3::operator*(mat3 A) const {
	return VectorFieldR3([X=_field, A](vec3 v) {return A * X(v); });
}

VectorFieldR3 VectorFieldR3::radial(vec3 scale) {
	return VectorFieldR3([scale](vec3 v) {return vec3(v.x*scale.x, v.y*scale.y, v.z*scale.z)/norm(v); });
}

vec3 VectorFieldR3::operator()(vec3 v) const {
	return _field(v);
}

VectorFieldR3 VectorFieldR3::operator+(VectorFieldR3 g) const {
	return VectorFieldR3([f=_field, g](vec3 v) {return f(v) + g(v); });
}

VectorFieldR3 VectorFieldR3::operator*(float a) const {
	return VectorFieldR3([f=_field, a](vec3 v) {return f(v) * a; });
}

VectorFieldR3 VectorFieldR3::operator-(VectorFieldR3 g) const {
	return VectorFieldR3([f=_field, g](vec3 v) {return f(v) - g(v); });
}

VectorFieldR3 VectorFieldR3::operator*(SmoothRealFunctionR3 g) const {
	return VectorFieldR3([f=_field, g](vec3 v) {return f(v) * g(v); });
}

VectorFieldR3 VectorFieldR3::constant(vec3 v) {
	return VectorFieldR3([v](vec3 x) {return v; });
}

VectorFieldR3 VectorFieldR3::linear(mat3 A) {
	return VectorFieldR3([A](vec3 v) {return A * v; });
}

VectorFieldR3 VectorFieldR3::rotational() {
	return VectorFieldR3([](vec3 v) {return vec3(-v.y, v.x, 0); });
}

VectorFieldR3 VectorFieldR3::wirlpool() {
	return VectorFieldR3([](vec3 v) {return vec3(-v.y, v.x, 0)/norm(vec2(v)); });
}

VectorFieldR3 VectorFieldR3::gradient(SmoothRealFunctionR3 f) {
	return f.gradient();
}


SpaceAutomorphism::SpaceAutomorphism(std::function<vec3(vec3)> f, std::function<vec3(vec3)> f_inv, std::function<mat3(vec3)> df) : SpaceEndomorphism(f, df) {
	this->_f_inv = f_inv;
}

SpaceAutomorphism::SpaceAutomorphism(std::function<vec3(vec3)> f, std::function<vec3(vec3)> f_inv, float epsilon) : SpaceEndomorphism(f, epsilon) {
	this->_f_inv = f_inv;
}

vec3 SpaceAutomorphism::inv(vec3 v) const {
	return _f_inv(v);
}

SpaceAutomorphism SpaceAutomorphism::operator~() const {
	return SpaceAutomorphism(_f_inv, _f, [d=_df](vec3 x) {return inverse(d(x)); });
}

SpaceAutomorphism SpaceAutomorphism::compose(SpaceAutomorphism g) const {
    return SpaceAutomorphism([f = _f, g](vec3 v) { return f(g(v)); }, [f_inv = _f_inv, g](vec3 v) { return g.inv(f_inv(v)); },
                             [d = _df, g](vec3 x) { return d(g(x)) * g.df(x); });
}


Complex Biholomorphism::f_inv(Complex z) const {
	return (*_f_inv)(z);
}

Biholomorphism Biholomorphism::operator~() const {
	auto df_cpy = _df;
	auto f_inv_cpy = _f_inv;
	return Biholomorphism(_f_inv, make_shared<endC>([df_cpy, f_inv_cpy](Complex z) {return ONE / (*df_cpy)((*f_inv_cpy)(z)); }), _f);
}

Biholomorphism::operator PlaneAutomorphism() const {
	auto new_f = _f;
	auto new_df = _df;
	auto new_f_inv = _f_inv;
	return PlaneAutomorphism(make_shared<std::function<vec2(vec2)>>([new_f](vec2 v) {return vec2((*new_f)(Complex(v))); }),
		make_shared<std::function<mat2(vec2) >>([new_df](vec2 x) {return mat2((*new_df)(Complex(x)).re(), (*new_df)(Complex(x)).im(), -(*new_df)(Complex(x)).im(), (*new_df)(Complex(x)).re()); }),
		make_shared<std::function<vec2(vec2)>>([new_f_inv](vec2 v) {return vec2((*new_f_inv)(Complex(v))); }));
}

Biholomorphism Biholomorphism::compose(Biholomorphism g) const {
	auto new_f = _f;
	auto new_df = _df;
	auto new_inv = _f_inv;
	auto g_f = g._f;
	auto g_df = g._df;
	auto g_inv = g._f_inv;
	return Biholomorphism(make_shared<endC>([new_f, g_f](Complex z) {return (*new_f)((*g_f)(z)); }),
		make_shared<endC>([new_f, new_df, g_f, g_df](Complex z) {return (*new_df)((*g_f)(z)) * (*g_df)(z); }),
		make_shared<endC>([new_inv, g_inv](Complex z) {return (*g_inv)((*new_inv)(z)); }));
}

//ComplexCurve::operator SmoothParametricPlaneCurve()
//{
//	auto new_f = *_f;
//	return SmoothParametricPlaneCurve([new_f](float t) {
//		return vec2(new_f(t));
//	}, t0, t1, period, epsilon);
//}