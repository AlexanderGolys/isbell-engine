#include "func.hpp"
#include <functional>
#include <vector>
#include <iostream>
#include <memory>


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
std::function<glm::vec2(float)> derivativeOperator(std::function<glm::vec2(float)> f, float epsilon) {
	auto fx = [f](float x) {return f(x).x; };
	auto fy = [f](float x) {return f(x).y; };
	return [fx, fy, epsilon](float x) {return vec2(derivativeOperator(fx, epsilon)(x), derivativeOperator(fy, epsilon)(x)); };
}
std::function<glm::vec3(float)> derivativeOperator(std::function<glm::vec3(float)> f, float epsilon) {
	auto fx = [f](float x) {return f(x).x; };
	auto fy = [f](float x) {return f(x).y; };
	auto fz = [f](float x) {return f(x).z; };
	return [fx, fy, fz, epsilon](float x) {return vec3(derivativeOperator(fx, epsilon)(x), derivativeOperator(fy, epsilon)(x), derivativeOperator(fz, epsilon)(x)); };
}



SmoothRealFunctionR3::SmoothRealFunctionR3() {
	this->_f = [](vec3 v) {return 0; };
	this->_df = [](vec3 v) {return vec3(0, 0, 0); };
}


SmoothRealFunctionR3::SmoothRealFunctionR3(std::function<float(glm::vec3)> f, std::function<glm::vec3(glm::vec3)> df) {
	this->_f = f;
	this->_df = df;
}

SmoothRealFunctionR3::SmoothRealFunctionR3(std::function<float(glm::vec3)> f, float epsilon) {
	this->_f = f;
	this->_df = [f, epsilon](vec3 v) {
		return  (f(v + vec3(epsilon, 0, 0)) - f(v - vec3(epsilon, 0, 0))) / (2 * epsilon) * vec3(1, 0, 0) +
			(f(v + vec3(0, epsilon, 0)) - f(v - vec3(0, epsilon, 0))) / (2 * epsilon) * vec3(0, 1, 0) +
			(f(v + vec3(0, 0, epsilon)) - f(v - vec3(0, 0, epsilon))) / (2 * epsilon) * vec3(0, 0, 1);
	};
}

float SmoothRealFunctionR3::operator()(glm::vec3 v) const {
	return _f(v);
}

glm::vec3 SmoothRealFunctionR3::df(glm::vec3 v) const {
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

SmoothRealFunctionR3 SmoothRealFunctionR3::linear(glm::vec3 v) {
	return SmoothRealFunctionR3([v](vec3 x) {return dot(x, v); },
		[v](vec3 x) {return v; });
}

SmoothRealFunctionR3 SmoothRealFunctionR3::projection(int i) {
	return SmoothRealFunctionR3([i](vec3 x) {return x[i]; },
		[i](vec3 x) {return vec3(0, 0, 0); });
}


SmoothRealFunctionR3 SmoothRealFunctionR3::constant(float a) {
	return SmoothRealFunctionR3([a](vec3 x) {return a; },
		[a](vec3 x) {return vec3(0, 0, 0); });
}

ParametricCurve::ParametricCurve()
{
	this->_f = [](float t) {return vec3(0, 0, 0); };
}

ParametricCurve::ParametricCurve(std::function<vec3(float)> curve)
{
	this->_f = curve;
}

vec3 ParametricCurve::operator()(float t) const {
    return this->_f(t);
}

SmoothParametricCurve::SmoothParametricCurve()
{
	this->_f = [this](float t) {return vec3(0, 0, 0); };
	this->_df = [this](float t) {return vec3(0, 0, 0); };
	this->_ddf = [this](float t) {return vec3(0, 0, 0); };
	this->regularity = SMOOTH;
}

SmoothParametricCurve::SmoothParametricCurve(std::function<glm::vec3(float)> f, std::function<glm::vec3(float)> df,
	std::function<glm::vec3(float)> ddf, float epsilon) {
	this->_f = f;
	this->_df = df;
	this->_ddf = ddf;
	this->regularity = SMOOTH;
	eps = epsilon;
}



SmoothParametricCurve::SmoothParametricCurve(std::function<vec3(float)> f, std::function<vec3(float)> df, float epsilon)

{
	this->_f = f;
	this->_df = df;
	this->_ddf = derivativeOperator(_df, epsilon);
	this->eps = epsilon;
	this->regularity = SMOOTH;
}

SmoothParametricCurve::SmoothParametricCurve(std::function<vec3(float)> f, float epsilon)
{
	this->_f = f;
	this->_df = derivativeOperator(_f, epsilon);
	this->_ddf = derivativeOperator(_df, epsilon);
	this->eps = epsilon;
	this->regularity = SMOOTH;
}

SmoothParametricCurve::SmoothParametricCurve(std::function<glm::vec3(float)> f,
	std::function<std::function<glm::vec3(float)>(int)> derivativeOperator, float epsilon) {
	this->_f = f;
	this->_df = derivativeOperator(1);
	this->_ddf = derivativeOperator(2);
	this->eps = epsilon;
	this->regularity = SMOOTH;
	this->_der_higher = derivativeOperator;
}

SmoothParametricCurve::SmoothParametricCurve(std::function<glm::vec3(float)> f,
	std::vector<std::function<glm::vec3(float)>> derivatives, float epsilon) {
	this->_f = f;
	if (derivatives.size() > 0)
		this->_df = derivatives[0];
	else
		this->_df = derivativeOperator(_f, epsilon);

	if (derivatives.size() > 1)
		this->_ddf = derivatives[1];
	else
		this->_ddf = derivativeOperator(_df, epsilon);

	for (int i = 2; i < derivatives.size(); i++) {
		auto d = derivatives[i];
		auto _der_higher_prev = this->_der_higher;
		this->_der_higher = [d, _der_higher_prev, i](int n) {return n==i+1 ? d : _der_higher_prev(n); };
	}
	this->eps = epsilon;
}

AffinePlane SmoothParametricCurve::osculatingPlane(float t) const {
	return AffinePlane(binormal(t), dot(binormal(t), _f(t)));
}

SmoothParametricCurve SmoothParametricCurve::constCurve(glm::vec3 v) {
	return SmoothParametricCurve([v](float t) {return v; },
		[](float t) {return vec3(0, 0, 0); },
		[](float t) {return vec3(0, 0, 0); });
}

float SmoothParametricCurve::length(float t0, float t1, int n) const {
	float res = 0;
	float dt = (t1 - t0) / n;
	for (int i = 0; i < n; i++)
	{
		res += norm(_f(lerp(t0, t1, dt*i)) - _f(lerp(t0, t1, dt*(i + 1))));
	}
	return res;
}

SmoothParametricCurve SmoothParametricCurve::operator*(SpaceEndomorphism g) const {
	auto new_f = _f;
	auto new_df = _df;
	auto new_ddf = _ddf;

	return SmoothParametricCurve([new_f, g](float t) {return g(new_f(t)); },
		[new_f, new_df, g](float t) {return g.df(new_f(t)) * new_df(t); },
		[new_f, new_df, new_ddf, g](float t) {return g.df(new_f(t)) * new_ddf(t) * g.df(new_f(t)); });
}

void SmoothParametricCurve::operator*=(SpaceEndomorphism g) {
	_f = [this, g](float t) {return g(this->_f(t)); };
	_df = [this, g](float t) {return g.df(this->_f(t)) * this->_df(t); };
	_ddf = [this, g](float t) {return g.df(this->_f(t)) * this->_ddf(t) * g.df(this->_f(t)); };
}

glm::mat3 SmoothParametricCurve::FrenetFrame(float t) const {
	return mat3(tangent(t), normal(t), binormal(t));

}

float SmoothParametricCurve::curvature(float t) const {
	return norm(cross(_df(t), _ddf(t))) / pow(norm(_df(t)), 3);
}

float SmoothParametricCurve::torsion(float t) const {
	return dot(cross(df(t), ddf(t)), higher_derivative(t, 3)) / norm2(cross(df(t), ddf(t)));
}

glm::vec3 SmoothParametricCurve::curvature_vector(float t) const {
	return normal(t)*curvature(t);
}

AffineLine::AffineLine(glm::vec3 p0, glm::vec3 v) :
SmoothParametricCurve([p0, v](float t) {return p0 + t * v; },
					  [v](float t) {return v; },
					  [](float t) {return vec3(0, 0, 0); }) {
	this->p0 = p0;
	this->v = v;
}

AffineLine AffineLine::spanOfPts(glm::vec3 p0, glm::vec3 p1) {
	return AffineLine(p0, p1 - p0);
}


SmoothRealFunctionR3 AffineLine::distanceField() const {
	return SmoothRealFunctionR3([this](vec3 p) {return distance(p); });
}

glm::vec3 AffineLine::orthogonalProjection(glm::vec3 p) const {
	return p0 + dot(p - p0, v) * v / dot(v, v);
}

glm::vec3 AffineLine::pivot() const {
	return p0;
}

glm::vec3 AffineLine::direction() const {
	return v;
}

float AffineLine::distance(AffineLine &l) const {
	return dot((p0 - l.pivot()), cross(v, l.direction())) / norm(cross(v, l.direction()));
}

AffineLine AffineLine::operator+(glm::vec3 v) const {
	return AffineLine(p0 + v, this->v);
}

SmoothParametricSurface::SmoothParametricSurface(std::function<glm::vec3(float, float)> f,
	std::function<glm::vec3(float, float)> df_t, std::function<glm::vec3(float, float)> df_u, float t0, float t1,
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

SmoothParametricSurface::SmoothParametricSurface(std::function<glm::vec3(float, float)> f, float t0, float t1, float u0,
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

glm::vec3 SmoothParametricSurface::operator()(float t, float s) const {
	return _f(t, s);
}

SmoothImplicitSurface::SmoothImplicitSurface(SmoothRealFunctionR3 F) {
	_F = F;
}

float SmoothImplicitSurface::operator()(glm::vec3 p) const {
	return _F(p);
}

AffinePlane::AffinePlane(glm::vec3 n, float d) :
SmoothImplicitSurface(SmoothRealFunctionR3(
	[n, d](vec3 p) {return dot(p, n) - d; }, [n](vec3 p) {return n; })) {
	this->n = normalise(n);
	this->d = d;
	auto tangents = orthogonalComplementBasis(n);
	this->v1 = tangents.first;
	this->v2 = tangents.second;
	this->pivot = this->n * d;
}

AffinePlane::AffinePlane(glm::vec3 pivot, glm::vec3 v1, glm::vec3 v2):
SmoothImplicitSurface(SmoothRealFunctionR3(
	[v1, v2, pivot](vec3 p) {return dot(p, normalise(cross(v1, v2))) - dot(pivot, normalise(cross(v1, v2))); },
	[v1, v2](vec3 p) {return normalise(cross(v1, v2)); })) {
	this->pivot = pivot;
	this->v1 = v1;
	this->v2 = v2;
	this->n = normalise(cross(v1, v2));
	this->d = dot(pivot, n);
}


AffinePlane AffinePlane::spanOfPts(glm::vec3 p0, glm::vec3 p1, glm::vec3 p2) {
	return AffinePlane(cross(p1-p0, p2-p0), dot(p0, cross(p1-p0, p2-p0)));
}

AffineLine AffinePlane::intersection(AffinePlane &p) const {
	vec3 v = cross(n, p.normal());
	float denominator = norm2(n)*norm2(p.normal()) - dot(n, p.normal())*dot(n, p.normal());
	vec3 pivot = n*(d*norm2(p.normal()) - p.getD()*dot(n, p.normal()))/denominator +
				p.normal()*(p.getD()*norm2(n) - d*dot(n, p.normal()))/denominator;
	return AffineLine(pivot, v);
}

float AffinePlane::distance(glm::vec3 p) const {
	return (*this)(p);
}

glm::vec3 AffinePlane::orthogonalProjection(glm::vec3 p) const {
	return p - (*this)(p) * n;
}

bool AffinePlane::contains(glm::vec3 p, float eps) const {
	return abs(distance(p)) < eps;
}

glm::vec3 AffinePlane::normal() const {
	return n;
}

float AffinePlane::getD() const {
	return d;
}

std::pair<glm::vec3, float> AffinePlane::equationCoefs() const {
	return std::pair<glm::vec3, float>(n, d);
}

glm::vec3 AffinePlane::intersection(AffineLine &l) const {
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

glm::mat3 AffinePlane::pivotAndBasis() const {
	return glm::mat3(pivot, v1, v2);
}

glm::vec2 AffinePlane::localCoordinates(glm::vec3 p) const {
	return glm::vec2(dot(p - pivot, v1), dot(p - pivot, v2));
}

float AffineLine::distance(glm::vec3 p) const {
	return norm(cross(p - p0, v)) / norm(v);
}

bool AffineLine::contains(glm::vec3 p, float eps) const {
	return distance(p) < eps;
}


ParametricPlanarCurve::ParametricPlanarCurve(std::function<vec2(float)> curve, float t0, float t1, float period, float epsilon)
{
	this->f =  std::make_unique<std::function<vec2(float)>>(curve);
	this->cyclic = period > 0;
	this->epsilon = epsilon;
	this->t0 = t0;
	this->t1 = t1;

	this->df = std::make_unique<std::function<vec2(float)>>([this](float t) {
		return ((*f)(t + this->epsilon) - (*f)(t)) / this->epsilon;
	});

	this->ddf = std::make_unique<std::function<vec2(float)>>([this](float t) {
		return ((*f)(t + this->epsilon) - (*f)(t)*2.f + (*f)(t - this->epsilon)) / (this->epsilon * this->epsilon);
	});

	this->N = std::make_unique<std::function<vec2(float)>>([this](float t) {
		return normalise<vec2>((*ddf)(t)); });

	this->period = period;
}

vec2 ParametricPlanarCurve::operator()(float t)
{
	return (*f)(t);
}

vector<vec2> ParametricPlanarCurve::sample(float t0, float t1, int n)
{
	vector<vec2> result = std::vector<vec2>();
	result.reserve(n);
	float dt = (t1 - t0) / (n-1);
	for (int i = 0; i < n; i++)
	{
		result.push_back((*f)(t0 + dt*i));
	}
	return result;
}

vector<vec2> ParametricPlanarCurve::sample(int n)
{
	return sample(t0, t1, n);
}

vector<vec3> ParametricPlanarCurve::adjacency_lines_buffer(float t0, float t1, int n, float z) const {
	vector<vec3> result = vector<vec3>();
	result.reserve(4 * n - 4);
	float dt = (t1 - t0) / n;
	for (int i = 1; i < n; i++)
	{
		vec3 p1 = vec3((*f)(t0 + (i - 1) * dt), z);
		vec3 p2 = vec3((*f)(t0 + i * dt), z);
		vec3 p0 = p1 + vec3((*N)(t0 + (i - 1) * dt), 0);
		vec3 p3 = p2 + vec3((*N)(t0 + i * dt), 0);
		result.push_back(p0);
		result.push_back(p1);
		result.push_back(p2);
		result.push_back(p3);
	}
	return result;
}

ComplexCurve::ComplexCurve(std::function<Complex(float)> curve, float t0, float t1, float period, float epsilon)
{
	this->f = std::make_unique<std::function<Complex(float)>>(curve);
	this->cyclic = period > 0;
	this->epsilon = epsilon;
	this->t0 = t0;
	this->t1 = t1;

	this->df = std::make_unique<std::function<Complex(float)>>([this](float t) {
		return ((*f)(t + this->epsilon) - (*f)(t)) / this->epsilon;
	});

	this->ddf = std::make_unique<std::function<Complex(float)>>([this](float t) {
		return ((*f)(t + this->epsilon) - (*f)(t)*2.f + (*f)(t - this->epsilon)) / (this->epsilon * this->epsilon);
	});

	this->N = std::make_unique<std::function<Complex(float)>>([this](float t) {
		return Complex(normalise<vec2>((*ddf)(t))); });

	this->period = period;
}

ComplexCurve::ComplexCurve(ParametricPlanarCurve* curve)
{
	std::function<vec2(float)> new_f = std::function<vec2(float)>(*curve->f);

	this->f = std::make_unique<std::function<Complex(float)>>([new_f](float t) {
		return Complex(new_f(t));
	});
	this->cyclic = curve->cyclic;
	this->epsilon = curve->epsilon;
	this->t0 = curve->t0;
	this->t1 = curve->t1;

	this->df = std::make_unique<std::function<Complex(float)>>([this](float t) {
		return ((*f)(t + this->epsilon) - (*f)(t)) / this->epsilon;
	});

	this->ddf = std::make_unique<std::function<Complex(float)>>([this](float t) {
		return ((*f)(t + this->epsilon) - (*f)(t)*2.f + (*f)(t - this->epsilon)) / (this->epsilon * this->epsilon);
	});

	this->N = std::make_unique<std::function<Complex(float)>>([this](float t) {
		return Complex(normalise<vec2>((*ddf)(t))); });

	this->period = curve->period;
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
	return Biholomorphism::mobius(Matrix<Complex, 2>(a, b, ZERO, ONE));
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



SpaceEndomorphism::SpaceEndomorphism(std::function<glm::vec3(glm::vec3)> f, std::function<glm::mat3(glm::vec3)> df) {
	_f = f;
	_df = df;
}

SpaceEndomorphism::SpaceEndomorphism(std::function<glm::vec3(glm::vec3)> f, float epsilon) {
	_f = f;
	auto dir_der = [f, epsilon](glm::vec3 x, glm::vec3 v) {
		return (f(x +  v*epsilon) - f(x - v*epsilon)) / (2 * epsilon);
	};
	_df = [dir_der](glm::vec3 x) {
		return mat3(dir_der(x, vec3(1, 0, 0)), dir_der(x, vec3(0, 1, 0)), dir_der(x, vec3(0, 0, 1)));
	};
}

glm::vec3 SpaceEndomorphism::directional_derivative(glm::vec3 x, glm::vec3 v) const {
	return _df(x) * v;
}

glm::vec3 SpaceEndomorphism::operator()(glm::vec3 v) const {
	return _f(v);
}

SpaceEndomorphism SpaceEndomorphism::compose(SpaceEndomorphism g) const {
	return SpaceEndomorphism([this, g](glm::vec3 v) {return this->_f(g(v)); },
			[this, g](glm::vec3 x) {return this->_df(g(x)) * g._df(x); });
}

SpaceEndomorphism SpaceEndomorphism::linear(glm::mat3 A) {
	return SpaceEndomorphism([A](glm::vec3 v) {return A * v; },
		[A](glm::vec3 x) {return A; });
}

SpaceEndomorphism SpaceEndomorphism::translation(glm::vec3 v) {
	return affine(mat3(1), v);
}

SpaceEndomorphism SpaceEndomorphism::scaling(float x, float y, float z) {
	return linear(mat3(x, 0, 0, 0, y, 0, 0, 0, z));
}

SpaceEndomorphism SpaceEndomorphism::affine(glm::mat3 A, glm::vec3 v) {
	return SpaceEndomorphism([A, v](glm::vec3 x) {return A * x + v; },
		[A](glm::vec3 x) {return A; });
}

VectorFieldR3::VectorFieldR3() {
	_field = [](vec3 v) {return vec3(0, 0, 0); };
}

VectorFieldR3::VectorFieldR3(std::function<glm::vec3(glm::vec3)> field) {
	_field = field;
}

VectorFieldR3::VectorFieldR3(SpaceEndomorphism f) {
	_field = f;
}

VectorFieldR3::VectorFieldR3(VectorFieldR2 field) {
	_field = [field](vec3 v) {return vec3(field(vec2(v.x, v.y)), 0); };
}

VectorFieldR3 VectorFieldR3::operator*(glm::mat3 A) const {
	return VectorFieldR3([this, A](vec3 v) {return A * this->_field(v); });
}

VectorFieldR3 VectorFieldR3::radial(glm::vec3 scale) {
	return VectorFieldR3([scale](vec3 v) {return vec3(v.x*scale.x, v.y*scale.y, v.z*scale.z)/norm(v); });
}

glm::vec3 VectorFieldR3::operator()(glm::vec3 v) const {
	return _field(v);
}

VectorFieldR3 VectorFieldR3::operator+(VectorFieldR3 g) const {
	return VectorFieldR3([this, g](vec3 v) {return this->_field(v) + g(v); });
}

VectorFieldR3 VectorFieldR3::operator*(float a) const {
	return VectorFieldR3([this, a](vec3 v) {return this->_field(v) * a; });
}

VectorFieldR3 VectorFieldR3::operator-(VectorFieldR3 g) const {
	return VectorFieldR3([this, g](vec3 v) {return this->_field(v) - g(v); });
}

VectorFieldR3 VectorFieldR3::operator*(SmoothRealFunctionR3 f) const {
	return VectorFieldR3([this, f](vec3 v) {return this->_field(v) * f(v); });
}

VectorFieldR3 VectorFieldR3::constant(glm::vec3 v) {
	return VectorFieldR3([v](vec3 x) {return v; });
}

VectorFieldR3 VectorFieldR3::linear(glm::mat3 A) {
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


SpaceAutomorphism::SpaceAutomorphism(std::function<glm::vec3(glm::vec3)> f, std::function<glm::vec3(glm::vec3)> f_inv, std::function<glm::mat3(glm::vec3)> df) : SpaceEndomorphism(f, df) {
	this->_f_inv = f_inv;
}

SpaceAutomorphism::SpaceAutomorphism(std::function<glm::vec3(glm::vec3)> f, std::function<glm::vec3(glm::vec3)> f_inv, float epsilon) : SpaceEndomorphism(f, epsilon) {
	this->_f_inv = f_inv;
}

glm::vec3 SpaceAutomorphism::inv(glm::vec3 v) const {
	return _f_inv(v);
}

SpaceAutomorphism SpaceAutomorphism::operator~() const {
	return SpaceAutomorphism(_f_inv, _f, [this](glm::vec3 x) {return inverse(_df(x)); });
}

SpaceAutomorphism SpaceAutomorphism::compose(SpaceAutomorphism g) const {
	return SpaceAutomorphism([this, g](glm::vec3 v) {return this->_f(g(v)); },
		[this, g](glm::vec3 v) {return g.inv(this->inv(v)); },
		[this, g](glm::vec3 x) {return this->_df(g(x)) * g._df(x); });
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

//ComplexCurve::operator ParametricPlanarCurve()
//{
//	auto new_f = *_f;
//	return ParametricPlanarCurve([new_f](float t) {
//		return vec2(new_f(t));
//	}, t0, t1, period, epsilon);
//}