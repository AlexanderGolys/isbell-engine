#pragma once
#include "utils/flows.hpp"


class SmoothParametricCurve;

class SmoothParametricPlaneCurve {
    Foo12 _f;
    Foo12 _df;
    Foo12 _ddf;
    std::function<Foo12(int)> _der_higher = [this](int n) {return n == 0 ? _f : n == 1 ? _df : n == 2 ? _ddf : derivativeOperator(_der_higher(n-1), this->eps);};
    float eps = 0.01;
    optional<float> t0 = std::nullopt;
    optional<float> t1 = std::nullopt;
    bool periodic = true;
public:
	explicit SmoothParametricPlaneCurve(const Foo12& curve, float t0=0, float t1=TAU, bool period=true, float epsilon=0.01);
    SmoothParametricPlaneCurve(Foo12 f,const Foo12& df, float t0=0, float t1=TAU, bool period=true, float epsilon=0.01);
    SmoothParametricPlaneCurve(Foo12 f, Foo12 df, Foo12 ddf, float t0 = 0, float t1 = TAU, bool period = true, float epsilon = 0.01);
    SmoothParametricPlaneCurve(const SmoothParametricPlaneCurve &other);
    SmoothParametricPlaneCurve(SmoothParametricPlaneCurve &&other) noexcept;
    SmoothParametricPlaneCurve &operator=(const SmoothParametricPlaneCurve &other);
    SmoothParametricPlaneCurve &operator=(SmoothParametricPlaneCurve &&other) noexcept;
    vec2 operator()(float t) const { return _f(t); }
    vec2 derivative(float t) const { return _df(t); }
    vec2 df(float t) const { return derivative(t); }
    vec2 second_derivative(float t) const { return _ddf(t); }
    vec2 ddf(float t) const { return second_derivative(t); }
    vec2 higher_derivative(float t, int n) const { return _der_higher(n)(t); }
    vec2 tangent(float t) const { return normalise(_df(t)); }
    vec2 normal(float t) const { return orthogonalComplement(tangent(t)); }

	SmoothParametricPlaneCurve operator+(vec2 v) const;
	SmoothParametricPlaneCurve operator-(vec2 v) const;

	vector<vec2> sample(float t0, float t1, int n) const;
	vector<vec2> sample(int n) const {return sample(t0.value_or(-1), t1.value_or(1), n);}
	vector<vec3> adjacency_lines_buffer(float t0, float t1, int n, float z=0) const;
    SmoothParametricCurve embedding(vec3 v1=e1, vec3 v2=e2, vec3 pivot=ORIGIN_R3) const;
    vec2 bounds() const { return vec2(t0.value_or(-1), t1.value_or(1)); }
	bool isPeriodic() const { return periodic; }
};

class PlanarSurface {
	HOM(vec2, vec2) f;
	CONST_PROPERTY(vec2, uBounds);
	CONST_PROPERTY(vec2, vBounds);
	CONST_PROPERTY(bool, uPeriodic);
	CONST_PROPERTY(bool, vPeriodic);
public:
	PlanarSurface(const HOM(vec2, vec2)& f, vec2 u_bounds, vec2 v_bounds, bool u_periodic=false, bool v_periodic=false);

	vec2 operator()(float u, float v) const;
	vec2 operator()(vec2 uv) const;

	bool normalisedDomain() const;
	void normaliseDomain();
};