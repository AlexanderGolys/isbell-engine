#include "func.hpp"

#include <chrono>
#include <functional>
#include <iosfwd>
#include <iostream>
#include <memory>
#include <random>
#include <utility>
#include <vector>
#include <glm/gtx/transform.hpp>

using namespace glm;
using std::vector, std::shared_ptr, std::make_shared, std::max, std::min;


Regularity operator+(Regularity a, int b);

int randomInt() {
	auto seed = std::chrono::system_clock::now().time_since_epoch().count();
	auto generator = std::default_random_engine(seed);
	auto distribution = std::uniform_int_distribution<int>(0, 1215752192);
	return distribution(generator);
}

long randomLong() {
	auto seed = std::chrono::system_clock::now().time_since_epoch().count();
	auto generator = std::default_random_engine(seed);
	auto distribution = std::uniform_int_distribution<long>(0, 1000000000000000);
	return distribution(generator);
}

float randomFloat(float a, float b) {
	auto seed = std::chrono::system_clock::now().time_since_epoch().count()*randomInt();
	auto generator = std::default_random_engine(seed);
	auto distribution = std::uniform_real_distribution<float>(a, b);
	return distribution(generator);
}

std::string randomString() {
	return std::to_string(randomLong());
}

Regularity operator-(Regularity a, int b) { // NOLINT(*-no-recursion)
    if (b < 0)
        return a + (-b);
    if (INT(a) < 0 || INT(a) > 6)
        return a;
    return static_cast<Regularity>(std::max(-1, INT(a) - b));
}

Regularity operator+(Regularity a, int b) { // NOLINT(*-no-recursion)
    if (b < 0)
        return a - (-b);
    if (INT(a) < -1 || INT(a) > 5)
        return a;
    return static_cast<Regularity>(std::min(6, INT(a) + b));
}




Foo31 partialDerivativeOperator(Foo31 f, int i, float epsilon) {
    return directionalDerivativeOperator<vec3, float>(std::move(f), vec3(i == 0, i == 1, i == 2), epsilon);
}

Foo21 partialDerivativeOperator(Foo21 f, int i, float epsilon) {
    return directionalDerivativeOperator<vec2, float>(std::move(f), vec2(i == 0, i == 1), epsilon);
}

Foo33 derivativeOperator(const Foo31 &f, float epsilon) {
    return [f, epsilon](vec3 v) { return
        vec3(partialDerivativeOperator(f, 0, epsilon)(v), partialDerivativeOperator(f, 1, epsilon)(v), partialDerivativeOperator(f, 2, epsilon)(v)); };
}

Foo22 derivativeOperator(const Foo21 &f, float epsilon) {
    return [f, epsilon](vec2 v) { return
        vec2(partialDerivativeOperator(f, 0, epsilon)(v), partialDerivativeOperator(f, 1, epsilon)(v)); };
}

Foo13 derivativeOperator(const Foo13  &f, float epsilon) {
    Fooo f1 = [f](float x) { return f(x).x; };
    Fooo f2 = [f](float x) { return f(x).y; };
    Fooo f3 = [f](float x) { return f(x).z; };
    return [f1, f2, f3, epsilon](float x) { return vec3(derivativeOperator(f1, epsilon)(x), derivativeOperator(f2, epsilon)(x), derivativeOperator(f3, epsilon)(x)); };
}

Foo12 derivativeOperator(const Foo12  &f, float epsilon) {
    Fooo f1 = [f](float x) { return f(x).x; };
    Fooo f2 = [f](float x) { return f(x).y; };
    return [f1, f2, epsilon](float x) { return vec2(derivativeOperator(f1, epsilon)(x), derivativeOperator(f2, epsilon)(x));};
}

RealFunctionR3::RealFunctionR3(RealFunctionR3 &&other) noexcept: _f(std::move(other._f)),
                                                                 _df(std::move(other._df)),
                                                                 eps(other.eps),
                                                                 regularity(other.regularity) {}

RealFunctionR3 & RealFunctionR3::operator=(const RealFunctionR3 &other) {
    if (this == &other)
        return *this;
    _f = other._f;
    _df = other._df;
    eps = other.eps;
    regularity = other.regularity;
    return *this;
}

RealFunctionR3 & RealFunctionR3::operator=(RealFunctionR3 &&other) noexcept {
    if (this == &other)
        return *this;
    _f = std::move(other._f);
    _df = std::move(other._df);
    eps = other.eps;
    regularity = other.regularity;
    return *this;
}

float RealFunctionR3::operator()(vec3 v) const {
	return _f(v);
}

vec3 RealFunctionR3::df(vec3 v) const {
	return _df(v);
}


RealFunctionR3 RealFunctionR3::operator*(float a) const {
	return RealFunctionR3([this, a](vec3 v) {return this->_f(v) * a; },
		[this, a](vec3 v) {return this->_df(v) * a; });
}

RealFunctionR3 RealFunctionR3::operator+(float a) const {
	return RealFunctionR3([this, a](vec3 v) {return this->_f(v) + a; },
		[this](vec3 v) {return this->_df(v); });
}

RealFunctionR3 RealFunctionR3::operator-(float a) const {
	return RealFunctionR3([this, a](vec3 v) {return this->_f(v) - a; },
		[this](vec3 v) {return this->_df(v); });
}

RealFunctionR3 RealFunctionR3::operator/(float a) const {
    return RealFunctionR3([this, a](vec3 v) { return this->_f(v) / a; }, [this, a](vec3 v) { return this->_df(v) / a; });
}

RealFunctionR3 RealFunctionR3::operator+(const RealFunctionR3 &g) const {
    return RealFunctionR3([f=_f, g_=g](vec3 x) { return f(x) + g_(x); },
                                [df=_df, g_=g](vec3 x) { return df(x) + g_.df(x); });
}


RealFunctionR3 RealFunctionR3::operator*(const RealFunctionR3 &g) const {
    return RealFunctionR3([f=_f, g_=g](vec3 x) { return f(x) * g_(x); },
                                [f=_f, df=_df, g_=g](vec3 x) { return f(x) * g_.df(x) + df(x) * g_(x); });
}

RealFunctionR3 RealFunctionR3::operator/(const RealFunctionR3 &g) const {
    return RealFunctionR3([f=_f, g_=g](vec3 x) { return f(x) / g_(x); },
                                [f=_f, df=_df, g_=g](vec3 x) { return (f(x) * g_.df(x) - df(x) * g_(x)) / (g_(x) * g_(x)); });
}

RealFunctionR3 RealFunctionR3::linear(vec3 v) {
	return RealFunctionR3([v](vec3 x) {return dot(x, v); },
		[v](vec3 x) {return v; });
}

RealFunctionR3 RealFunctionR3::projection(int i) {
    return RealFunctionR3([i](vec3 x) { return x[i]; }, [i](vec3 x) { return vec3(0, 0, 0); });
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


RealFunctionR3 RealFunctionR3::constant(float a) {
	return RealFunctionR3([a](vec3 x) {return a; },
		[a](vec3 x) {return vec3(0, 0, 0); });
}




RealFunctionR1::RealFunctionR1(const RealFunctionR1 &other): _f(other._f),
															 _df(other._df),
															 _ddf(other._ddf),
															 eps(other.eps) {}

RealFunctionR1::RealFunctionR1(RealFunctionR1 &&other) noexcept: _f(std::move(other._f)),
																 _df(std::move(other._df)),
																 _ddf(std::move(other._ddf)),
																 eps(other.eps) {}

RealFunctionR1 & RealFunctionR1::operator=(const RealFunctionR1 &other) {
	if (this == &other)
		return *this;
	_f   = other._f;
	_df  = other._df;
	_ddf = other._ddf;
	eps  = other.eps;
	return *this;
}

RealFunctionR1 & RealFunctionR1::operator=(RealFunctionR1 &&other) noexcept {
	if (this == &other)
		return *this;
	_f   = std::move(other._f);
	_df  = std::move(other._df);
	_ddf = std::move(other._ddf);
	eps  = other.eps;
	return *this;
}

RealFunctionR1 RealFunctionR1::operator+(const RealFunctionR1 &g) const { return RealFunctionR1([f=_f, g](float x) { return f(x) + g(x); },
																								[df=_df, dg=g._df](float x) { return df(x) + dg(x); },
																								[ddf=_ddf, ddg=g._ddf](float x) { return ddf(x) + ddg(x); }); }

RealFunctionR1 RealFunctionR1::operator*(float a) const { return RealFunctionR1([f=_f, a](float x) { return f(x) * a; },
																				[df=_df, a](float x) { return df(x) * a; },
																				[ddf=_ddf, a](float x) { return ddf(x) * a; }); }

RealFunctionR1 RealFunctionR1::operator+(float a) const { return RealFunctionR1([this, a](float x) {return this->_f(x) + a; },
																				[this](float x) {return this->_df(x); },
																				[this](float x) {return this->_ddf(x); }); }

RealFunctionR1 RealFunctionR1::operator*(const RealFunctionR1 &g_) const { return RealFunctionR1(
		[f=_f, g=g_._f](float x) { return f(x)*g(x); },
		[f=_f, g=g_._f, df=_df, dg=g_._df](float x) { return f(x)*dg(x) + g(x)*df(x); },
		[f=_f, g=g_._f, df=_df, dg=g_._df, ddf=_ddf, ddg=g_._ddf](float x) { return f(x)*ddg(x) + 2*df(x)*dg(x) + g(x)*ddf(x); }, eps); }

RealFunctionR1 RealFunctionR1::operator/(const RealFunctionR1 &g_) const { return RealFunctionR1(
		[f=_f, g=g_._f](float x) { return f(x)/g(x); },
		[f=_f, g=g_._f, df=_df, dg=g_._df](float x) { return (f(x)*dg(x) - g(x)*df(x)) / (g(x)*g(x)); },
		[f=_f, g=g_._f, df=_df, dg=g_._df, ddf=_ddf, ddg=g_._ddf](float x) {
			return (ddf(x)*g(x)*g(x) - 2*df(x)*dg(x)*g(x) + f(x)*ddg(x)*g(x)*g(x) - f(x)*g(x)*g(x)*g(x)*g(x))/(g(x)*g(x)*g(x)*g(x)); }, eps); }

RealFunctionR1 RealFunctionR1::operator&(const RealFunctionR1 &g_) const {
	return RealFunctionR1([f=_f, g=g_._f](float x) { return f(g(x)); },
						  [df=_df, dg=g_._df, g=g_._f](float x) { return df(g(x)) * dg(x); },
						  [ddf=_ddf, ddg=g_._ddf, dg=g_._df, g=g_._f, df=_df](float x) { return ddf(g(x)) * dg(x) * dg(x) + df(g(x)) * ddg(x); });
}

RealFunctionR1 RealFunctionR1::monomial(int n) { if (n == 0) return constant(1); return RealFunctionR1([n](float x) { return ::pow(x, n); }, [n](float x) { return n * ::pow(x, n - 1); }, [n](float x) { return n * (n - 1) * ::pow(x, n - 2); }); }

RealFunctionR1 RealFunctionR1::polynomial(std::vector<float> coeffs) {
	if (coeffs.size() == 1) return constant(coeffs[0]);
	std::vector<float> c_lower = rangeFrom(coeffs, 1);
	int degree                 = coeffs.size() - 1;
	return monomial(degree) * coeffs[0] + polynomial(c_lower);
}

RealLineAutomorphism RealLineAutomorphism::operator&(const RealLineAutomorphism &g) const {
	return RealLineAutomorphism(static_cast<RealFunctionR1>(*this) & static_cast<RealFunctionR1>(g), g.inverse & inverse);
}

RealLineAutomorphism RealLineAutomorphism::pow(int n) {
	if (n % 2 == 0) throw IllegalArgumentError("Only odd degree monomials are invertible");
	return RealLineAutomorphism(RealFunctionR1::pow(n), RealFunctionR1::pow(1.f/n));}

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

RealFunctionR1 smoothenReLu(float c0, float x_change_to_id) {
	return RealFunctionR1([c0, x_change_to_id](float t) {
		if (t>=x_change_to_id) return t;
		float a = 2.0*c0 - x_change_to_id;
		float b = 2.0*x_change_to_id - 3.0*c0;
		float x = t/x_change_to_id;
		return (a*x + b)*x*x + c0;
	});
}

RealFunctionR1 expImpulse(float peak) {
	return RealFunctionR1([peak](float t) {
		return t/peak*std::exp(1-t/peak);
	});
}

RealFunctionR1 cubicShroom(float center, float margin) {
	return RealFunctionR1([center, margin](float t) {
			if (abs(t-center) > margin) return 0.f;
			return 1.f - sq(t-center)*(3-2*abs(t-center)/margin)/sq(margin);
		});
}

RealFunctionR1 powerShroom(float begin, float end, float zeroExponent, float oneExponent) {
	if (zeroExponent <= 0) throw IllegalArgumentError("zeroExponent must be positive.");
	if (oneExponent <= 0) throw IllegalArgumentError("oneExponent must be positive.");
	return RealFunctionR1([begin, end, zeroExponent, oneExponent](float t) {
		float x = (t-begin)/(end-begin);
		if (x < 0) return 0.f;
		if (x > 1) return 0.f;
		float c = pow(zeroExponent+oneExponent, zeroExponent+oneExponent)/(pow(zeroExponent, zeroExponent) * pow(oneExponent, oneExponent));
		return c * pow(x, zeroExponent) * pow(1.f-x, oneExponent);
	});
}

RealFunctionR1 toneMap(float k) {
	return RealFunctionR1([k](float t) {
		return (k*t)/(1+k*t);
	});
}

RealFunctionR1 rationalInfiniteShroom(float steepFactor, float center) {
	if (steepFactor <= 0) throw IllegalArgumentError("steepFactor must be positive.");
	return RealFunctionR1([steepFactor, center](float t) {
			return 1.f/(steepFactor*sq(t-center) + 1);
		});
}

RealFunctionR1 expStep(int exponentDegree, float center) {
	if (exponentDegree < 1) throw IllegalArgumentError("Exponent should be positive (and at least 2).");
	if (exponentDegree < 2) throw IllegalArgumentError("Step with exponent 1 is not differentiable at zero.");
	return RealFunctionR1([exponentDegree, center](float t) {
			if (t < center) return 1.f;
			return (float)exp(-1.f*exp(exponentDegree)*pow(t, exponentDegree));
		});
}

RealFunctionR1 flatFunctionStep(float center) {
	return RealFunctionR1([ center](float t) {
			return t > center ? exp(-1/(t-center)) : 0;
		});
}

RealFunctionR1 flatAtRoofShroom(float center) {
	return RealFunctionR1([ center](float t) {
			return 1- exp(-1/sq(t-center));
		});
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


SpaceEndomorphism SpaceEndomorphism::compose(const SpaceEndomorphism &g) const {
	return SpaceEndomorphism([f=_f, g](vec3 v) {return f(g(v)); },
			[d=_df, g](vec3 x) {return d(g(x)) * g.df(x); });
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
	vec3 n = normalize(axis);
	float nx = n.x;
	float ny = n.y;
	float nz = n.z;
    return linear(rotationMatrix3(n, angle));
}


VectorFieldR3::VectorFieldR3() {
    _X = [](vec3 v) { return vec3(0, 0, 0); };
}
VectorFieldR3::VectorFieldR3(Foo33 X, float eps) : _X(std::move(X)), eps(eps) {
    RealFunctionR3 Fx = RealFunctionR3([X=_X](vec3 x) { return X(x).x;}, eps);
    RealFunctionR3 Fy = RealFunctionR3([X=_X](vec3 x) { return X(x).y;}, eps);
    RealFunctionR3 Fz = RealFunctionR3([X=_X](vec3 x) { return X(x).z;}, eps);
    _dX = [Fx, Fy, Fz](vec3 v) { return mat3(Fx.df(v), Fy.df(v), Fz.df(v)); };
}

VectorFieldR3::VectorFieldR3(RealFunctionR3 Fx, RealFunctionR3 Fy, RealFunctionR3 Fz, float epsilon) : eps(epsilon) {
    _X = [fx = Fx, fy = Fy, fz = Fz](vec3 v) { return vec3(fx(v), fy(v), fz(v)); };
    _dX = [fx = Fx, fy = Fy, fz = Fz](vec3 v) { return mat3(fx.df(v), fy.df(v), fz.df(v)); };
}


VectorFieldR3::VectorFieldR3(VectorFieldR2 f) : VectorFieldR3([_f=f](vec3 v) {return vec3(_f(vec2(v.x, v.y)), 0); }, 0.01) {}



VectorFieldR3 VectorFieldR3::constant(vec3 v) { return VectorFieldR3([v](vec3 x) {return v; }, 0.01); }


VectorFieldR3 VectorFieldR3::operator+(const VectorFieldR3 &Y) const {
    return VectorFieldR3([f=_X, g=Y._X](vec3 v) {return f(v) + g(v); },  [df=_dX, dg=Y._dX](vec3 v) {return df(v) + dg(v); }, eps);
}

VectorFieldR3 VectorFieldR3::operator*(float a) const {
    return VectorFieldR3([f=_X, a](vec3 v) {return f(v) * a; }, [df=_dX, a](vec3 v) {return df(v) * a; }, eps);
}

VectorFieldR3 VectorFieldR3::operator*(const RealFunctionR3 &f) const {
    return VectorFieldR3([X=_X, g=f](vec3 v) {return X(v) * g(v); }, [dX=_dX, g=f](vec3 v) {return dX(v) * g(v); }, eps);
}

VectorFieldR3 VectorFieldR3::linear(mat3 A) {
    return VectorFieldR3([A](vec3 v) {return A * v; }, [A](vec3 v) {return A; }, 0.01);
}

VectorFieldR3 VectorFieldR3::radial(vec3 scale) {
    return VectorFieldR3([scale](vec3 v) {return vec3(v.x*scale.x, v.y*scale.y, v.z*scale.z); }, [scale](vec3 v) {return mat3(scale.x, 0, 0, 0, scale.y, 0, 0, 0, scale.z); }, 0.01);
}


RealFunctionR3 VectorFieldR3::divergence() const {
    return RealFunctionR3([f = _X](vec3 v) {
        return dot(f(v), vec3(1, 1, 1));
    });
}

VectorFieldR3 VectorFieldR3::curl() const {
    RealFunctionR3 F1 = F_x();
    RealFunctionR3 F2 = F_y();
    RealFunctionR3 F3 = F_z();
    return VectorFieldR3([F1, F2, F3](vec3 v) {
        return vec3(F3.dy(v) - F2.dz(v), F1.dz(v) - F3.dx(v), F2.dx(v) - F1.dy(v));
    }, 0.01);
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
