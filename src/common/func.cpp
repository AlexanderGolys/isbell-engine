#include "func.hpp"
#include "mat.hpp"

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


std::string polyGroupIDtoString(std::variant<int, std::string> id) { return std::holds_alternative<int>(id) ? std::to_string(std::get<int>(id)) : std::get<std::string>(id); }

std::variant<int, std::string> prefix(const std::variant<int, std::string> &id, const std::string &prefix) { return PolyGroupID(prefix + polyGroupIDtoString(id)); }

std::variant<int, std::string> make_unique_id(const std::variant<int, std::string> &id) { return prefix(id, randomStringNumeric()); }

std::variant<int, std::string> bdGroup(int n) { return PolyGroupID("bd" + std::to_string(n));}

std::variant<int, std::string> curveGroup(int n) { return PolyGroupID("curva" + std::to_string(n));}

std::variant<int, std::string> randomCurvaID() { return PolyGroupID("curva" + randomStringNumeric());}

std::variant<int, std::string> randomID() { return PolyGroupID(randomStringNumeric());}

std::variant<int, std::string> randomBdID() { return PolyGroupID("bd" + randomStringNumeric()+ randomStringNumeric()+ randomStringNumeric());}

bool idBelongsToCurve(const std::variant<int, std::string> &id) { return std::holds_alternative<std::string>(id) && std::get<std::string>(id).starts_with("curva");}

bool idBelongsToBd(const std::variant<int, std::string> &id) { return std::holds_alternative<std::string>(id) && std::get<std::string>(id).starts_with("bd");}

Regularity max(Regularity a, Regularity b) {return static_cast<Regularity>(std::max(INT(a), INT(b))); }

Regularity min(Regularity a, Regularity b) {return static_cast<Regularity>(std::min(INT(a), INT(b))); }

Regularity operator--(Regularity a) { return a - 1; }

Regularity operator++(Regularity a) { return a + 1; }

template<typename codomain>
std::function<codomain(float)> derivativeOperator(std::function<codomain(float)> f, float epsilon) {
	return [f, epsilon](float x) {return (f(x + epsilon) - f(x - epsilon)) / (2 * epsilon); };
}

template<typename domain, typename codomain>
std::function<codomain(domain)> directionalDerivativeOperator(std::function<codomain(domain)> f, domain v, float epsilon) {
	return [f, epsilon, v](domain x) { return (f(x + v * epsilon) - f(x - v * epsilon)) / (2 * epsilon); };
}


VectorFieldR2::VectorFieldR2(): field([](vec2 v) { return vec2(0, 0); }) {}

VectorFieldR2::VectorFieldR2(const std::function<vec2(vec2)> &field): field(field) {}

vec2 VectorFieldR2::operator()(vec2 v) const { return field(v); }

vec2 VectorFieldR2::normal(vec2 x) const { return normalise(vec2(-field(x).y, field(x).x)); }

float VectorFieldR2::speed(vec2 x) const { return length(field(x)); }

RealFunctionR3::RealFunctionR3(): _f([](vec3 v) { return 0; }), _df([](vec3 v) { return vec3(0); }), regularity(Regularity::ANALYTIC) {}

RealFunctionR3::RealFunctionR3(const RealFunctionR3 &other) = default;

RealFunctionR3::RealFunctionR3(std::function<float(vec3)> f, std::function<vec3(vec3)> df, float eps, Regularity regularity): _f(f), _df(df), eps(eps), regularity(regularity) {}

RealFunctionR3::RealFunctionR3(std::function<float(vec3)> f, float epsilon, Regularity regularity):  _f(f), _df(derivativeOperator(f, epsilon)), eps(epsilon), regularity(regularity) {}

RealFunctionR3::RealFunctionR3(std::function<float(float, float, float)> f, std::function<vec3(float, float, float)> df, float eps, Regularity regularity): _f(wrap(f)), _df(wrap(df)), eps(eps), regularity(regularity) {}

RealFunctionR3::RealFunctionR3(std::function<float(float, float, float)> f, float epsilon, Regularity regularity):  _f(wrap(f)), _df(derivativeOperator(wrap(f), epsilon)), eps(epsilon), regularity(regularity) {}

float RealFunctionR3::operator()(float x, float y, float z) const { return _f(vec3(x, y, z)); }

vec3 RealFunctionR3::df(float x, float y, float z) const { return _df(vec3(x, y, z)); }

RealFunctionR3 RealFunctionR3::operator-() const { return *this * -1; }

RealFunctionR3 RealFunctionR3::operator-(const RealFunctionR3 &g) const { return *this + (-g); }

RealFunctionR3 operator*(float a, const RealFunctionR3 &f) { return f * a; }

RealFunctionR3 operator+(float a, const RealFunctionR3 &f) { return f + a; }

RealFunctionR3 operator-(float a, const RealFunctionR3 &f) { return -f + a; }

RealFunctionR3 operator/(float a, const RealFunctionR3 &f) { return RealFunctionR3::constant(a)/f; }

RealFunctionR3 RealFunctionR3::operator~() const { return 1/(*this); }

float RealFunctionR3::getEps() const { return eps; }

float RealFunctionR3::dx(vec3 x) const { return _df(x).x; }

float RealFunctionR3::dy(vec3 x) const { return _df(x).y; }

float RealFunctionR3::dz(vec3 x) const { return _df(x).z; }

SpaceEndomorphism::~SpaceEndomorphism() = default;

SpaceEndomorphism::SpaceEndomorphism(const SpaceEndomorphism &other): _f(other._f), _df(other._df) {}

SpaceEndomorphism::SpaceEndomorphism(SpaceEndomorphism &&other) noexcept: _f(std::move(other._f)), _df(std::move(other._df)) {}

SpaceEndomorphism::SpaceEndomorphism(std::function<vec3(vec3)> f, std::function<mat3(vec3)> df, float eps): _f(std::move(f)), _df(std::move(df)), eps(eps) {}

SpaceEndomorphism::SpaceEndomorphism(mat3 A): _f([A](vec3 x) { return A * x; }), _df([A](vec3 x) { return A; }) {}

SpaceEndomorphism::SpaceEndomorphism(mat4 A): _f([A](vec3 x) { return vec3(A * vec4(x, 1)); }), _df([A](vec3 x) { return mat3(A); }) {}

vec3 SpaceEndomorphism::directional_derivative(vec3 x, vec3 v) const { return _df(x) * v; }

vec3 SpaceEndomorphism::dfdv(vec3 x, vec3 v) const { return directional_derivative(x, v); }

mat3 SpaceEndomorphism::df(vec3 x) const { return _df(x); }

vec3 SpaceEndomorphism::operator()(vec3 x) const { return _f(x); }

SpaceEndomorphism SpaceEndomorphism::operator&(const SpaceEndomorphism &g) const { return compose(g); }

float SpaceEndomorphism::getEps() const { return eps; }

SpaceEndomorphism SpaceEndomorphism::linear(const mat3 &A) { return SpaceEndomorphism(A); }

PlaneEndomorphism::~PlaneEndomorphism() = default;

PlaneEndomorphism::PlaneEndomorphism(const PlaneEndomorphism &other): _f(other._f), _df(other._df) , eps(other.eps) {}

PlaneEndomorphism::PlaneEndomorphism(PlaneEndomorphism &&other) noexcept: _f(std::move(other._f)), _df(std::move(other._df)) {}

PlaneEndomorphism::PlaneEndomorphism(std::function<vec2(vec2)> f, std::function<glm::mat2(glm::vec2)> df, float eps): _f(std::move(f)), _df(std::move(df)), eps(eps) {}

PlaneEndomorphism::PlaneEndomorphism(mat2 A): _f([A](vec2 x) { return A * x; }), _df([A](vec2 x) { return A; }) {}

vec2 PlaneEndomorphism::directional_derivative(vec2 x, vec2 v) const { return _df(x) * v; }

vec2 PlaneEndomorphism::dfdv(vec2 x, vec2 v) const { return directional_derivative(x, v); }

mat2 PlaneEndomorphism::df(vec2 x) const { return _df(x); }

mat2 PlaneEndomorphism::df(float x, float y) const { return _df(vec2(x, y)); }

vec2 PlaneEndomorphism::operator()(vec2 x) const { return _f(x); }

vec2 PlaneEndomorphism::operator()(float x, float y) const { return _f(vec2(x, y)); }

PlaneEndomorphism PlaneEndomorphism::operator&(const PlaneEndomorphism &g) const { return compose(g); }

PlaneEndomorphism PlaneEndomorphism::linear(mat2 A) { return PlaneEndomorphism(A); }

PlaneEndomorphism::PlaneEndomorphism(std::function<vec2(vec2)> f, float epsilon) {
	_f = f;
    _df = [f, epsilon](vec2 x) {
        return mat2(
            (f(x + vec2(epsilon, 0)) - f(x - vec2(epsilon, 0))) / (2 * epsilon),
            (f(x + vec2(0, epsilon)) - f(x - vec2(0, epsilon))) / (2 * epsilon)
        );
    };
	eps = epsilon;
}

PlaneEndomorphism PlaneEndomorphism::compose(const PlaneEndomorphism &g) const {
	return PlaneEndomorphism([f=_f, g=g._f](vec2 x) { return f(g(x)); }, [df=_df, dg=g._df, g=g._f](vec2 x) { return df(g(x)) * dg(x); }, eps);
}

VectorFieldR3::VectorFieldR3(std::function<vec3(vec3)> X, std::function<mat3(vec3)> dX, float eps): _X(X), _dX(std::move(dX)), eps(eps) {}

VectorFieldR3 VectorFieldR3::operator-() const { return *this * -1; }

VectorFieldR3 VectorFieldR3::operator-(const VectorFieldR3 &Y) const { return *this + (-Y); }

RealFunctionR3 VectorFieldR3::F_x() const { return RealFunctionR3([this](vec3 x) { return _X(x).x; }, [this](vec3 x) { return _dX(x)[0]; }); }

RealFunctionR3 VectorFieldR3::F_y() const { return RealFunctionR3([this](vec3 x) { return _X(x).y; }, [this](vec3 x) { return _dX(x)[1]; }); }

RealFunctionR3 VectorFieldR3::F_z() const { return RealFunctionR3([this](vec3 x) { return _X(x).z; }, [this](vec3 x) { return _dX(x)[2]; }); }

std::array<RealFunctionR3, 3> VectorFieldR3::components() const { return {F_x(), F_y(), F_z()}; }

glm::vec3 VectorFieldR3::operator()(glm::vec3 v) const { return _X(v); }

VectorFieldR3 operator*(const mat3 &A, const VectorFieldR3 &X) {
	return VectorFieldR3([f=X._X, A](vec3 v) {return A * f(v); }, [df=X._dX, A](vec3 v) {return A * df(v); }, X.eps);
}

vec3 VectorFieldR3::moveAlong(vec3 v, float dt) const { return v + _X(v) * dt; }

VectorFieldR3Dynamic::VectorFieldR3Dynamic(std::function<vec3(float, vec3)> X, float epsilon): _X(X), eps(epsilon) {}

vec3 VectorFieldR3Dynamic::operator()(float t, vec3 x) const { return _X(t, x); }

VectorFieldR3Dynamic VectorFieldR3Dynamic::operator+(const VectorFieldR3Dynamic &Y) const { return VectorFieldR3Dynamic([__X=_X, _Y=Y._X](float t, vec3 x) { return  __X(t, x) + _Y(t, x); }, eps); }

VectorFieldR3Dynamic VectorFieldR3Dynamic::operator*(float a) const { return VectorFieldR3Dynamic([__X=_X, a](float t, vec3 x) { return  __X(t, x) * a; }, eps); }

SpaceAutomorphism SpaceAutomorphism::inv() const { return ~(*this); }

SpaceAutomorphism SpaceAutomorphism::operator&(const SpaceAutomorphism &g) const { return compose(g); }

SpaceAutomorphism SpaceAutomorphism::applyWithBasis(mat3 A) const { return linear(A) & *this & linear(inverse(A)); }

SpaceAutomorphism SpaceAutomorphism::applyWithBasis(vec3 v1, vec3 v2, vec3 v3) const { return applyWithBasis(mat3(v1, v2, v3)); }

SpaceAutomorphism SpaceAutomorphism::applyWithShift(vec3 v) const { return translation(v) & *this & translation(-v); }

SpaceAutomorphism SpaceAutomorphism::scaling(vec3 factors) { return scaling(factors.x, factors.y, factors.z); }

SpaceAutomorphism SpaceAutomorphism::scaling(float x) { return scaling(x, x, x); }

SpaceAutomorphism SpaceAutomorphism::scaling(float x, float y, float z, vec3 center) { return translation(center) & scaling(x, y, z) & translation(-center); }

SpaceAutomorphism SpaceAutomorphism::scaling(float x, vec3 center) { return scaling(x, x, x, center); }

SpaceAutomorphism SpaceAutomorphism::scaling(vec3 factors, vec3 center) {  return scaling(factors.x, factors.y, factors.z, center); }

SpaceAutomorphism SpaceAutomorphism::rotation(vec3 axis, float angle, vec3 center) { return rotation(axis, angle).applyWithShift(center); }

SpaceAutomorphism SpaceAutomorphism::deltaRotation(vec3 v1, vec3 v2) { return linear(rotationBetween(v1, v2)); }

SpaceAutomorphism SpaceAutomorphism::deltaRotation(vec3 v1, vec3 v2, vec3 center) { return translation(center) & deltaRotation(v1, v2) & translation(-center); }

RealFunctionR1::RealFunctionR1(std::function<float(float)> f, std::function<float(float)> df, std::function<float(float)> ddf, float epsilon): _f(f), _df(df), _ddf(ddf), eps(epsilon) {}

RealFunctionR1::RealFunctionR1(std::function<float(float)> f, std::function<float(float)> df, float epsilon):  RealFunctionR1(f, df, derivativeOperator(df, epsilon), epsilon) {}

RealFunctionR1::RealFunctionR1(std::function<float(float)> f, float epsilon): RealFunctionR1(f, derivativeOperator(f, epsilon), epsilon) {}

float RealFunctionR1::operator()(float x) const { return _f(x); }

float RealFunctionR1::df(float x) const { return _df(x); }

float RealFunctionR1::ddf(float x) const { return _ddf(x); }

RealFunctionR1 RealFunctionR1::df() const {return RealFunctionR1(_df, _ddf, eps);}

RealFunctionR1 RealFunctionR1::operator-(float a) const { return *this + (-a); }

RealFunctionR1 RealFunctionR1::operator-() const { return *this * -1; }

RealFunctionR1 RealFunctionR1::operator-(const RealFunctionR1 &g) const { return *this + (-g); }

RealFunctionR1 RealFunctionR1::operator/(float a) const { return *this * (1/a); }

RealFunctionR1 operator*(float a, const RealFunctionR1 &f) { return f * a; }

RealFunctionR1 operator+(float a, const RealFunctionR1 &f) { return f + a; }

RealFunctionR1 operator-(float a, const RealFunctionR1 &f) { return -f + a; }

RealFunctionR1 operator/(float a, const RealFunctionR1 &f) { return RealFunctionR1::constant(a)/f; }

RealFunctionR1 RealFunctionR1::constant(float a) { return RealFunctionR1([a](float x) { return a; }, [a](float x) { return 0; }, [a](float x) { return 0; }); }

RealFunctionR1 RealFunctionR1::x() { return RealFunctionR1([](float x) { return x; }, [](float x) { return 1; }, [](float x) { return 0; }); }

RealFunctionR1 RealFunctionR1::one() { return constant(1); }

RealFunctionR1 RealFunctionR1::zero() { return constant(0); }

RealFunctionR1 RealFunctionR1::linear(float a, float b) { return x() * a + b; }

RealFunctionR1 RealFunctionR1::quadratic(float a, float b, float c) { return x() * x() * a + x() * b + c; }

RealFunctionR1 RealFunctionR1::sin() { return RealFunctionR1([](float x) { return std::sin(x); }, [](float x) { return std::cos(x); }, [](float x) { return -std::sin(x); }); }

RealFunctionR1 RealFunctionR1::cos() { return RealFunctionR1([](float x) { return std::cos(x); }, [](float x) { return -std::sin(x); }, [](float x) { return -std::cos(x); }); }

RealFunctionR1 RealFunctionR1::exp() { return RealFunctionR1([](float x) { return std::exp(x); }, [](float x) { return std::exp(x); }, [](float x) { return std::exp(x); }); }

RealFunctionR1 RealFunctionR1::log() { return RealFunctionR1([](float x) { return std::log(x); }, [](float x) { return 1/x; }, [](float x) { return -1/(x*x); }); }

RealFunctionR1 RealFunctionR1::sqrt() { return RealFunctionR1([](float x) { return std::sqrt(x); }, [](float x) { return 1/(2*std::sqrt(x)); }, [](float x) { return -1/(4*x*std::sqrt(x)); }); }

RealFunctionR1 RealFunctionR1::pow(float a) { return RealFunctionR1([a](float x) { return std::pow(x, a); }, [a](float x) { return a * std::pow(x, a - 1); }, [a](float x) { return a * (a - 1) * std::pow(x, a - 2); }); }

RealFunctionR1 RealFunctionR1::pow(int a) { return RealFunctionR1([a](float x) { return std::pow(x, a); }, [a](float x) { return a * std::pow(x, a - 1); }, [a](float x) { return a * (a - 1) * std::pow(x, a - 2); }); }

RealLineAutomorphism::RealLineAutomorphism(RealFunctionR1 f, RealFunctionR1 inv): RealFunctionR1(f), inverse(inv) {}

RealLineAutomorphism::RealLineAutomorphism(std::function<float(float)> f, std::function<float(float)> f_inv, float epsilon): RealFunctionR1(f, epsilon), inverse(f_inv, epsilon) {}

RealLineAutomorphism::RealLineAutomorphism(std::function<float(float)> f, std::function<float(float)> df, std::function<float(float)> f_inv, float epsilon): RealFunctionR1(f, df, epsilon), inverse(f_inv, epsilon) {}

RealLineAutomorphism::RealLineAutomorphism(std::function<float(float)> f, std::function<float(float)> df, std::function<float(float)> f_inv, std::function<float(float)> df_inv, float epsilon): RealFunctionR1(f, df, epsilon), inverse(f_inv, df_inv, epsilon) {}

RealLineAutomorphism::RealLineAutomorphism(std::function<float(float)> f, std::function<float(float)> df, std::function<float(float)> ddf, std::function<float(float)> f_inv, std::function<float(float)> df_inv, std::function<float(float)> ddf_inv, float epsilon): RealFunctionR1(f, df, ddf, epsilon), inverse(f_inv, df_inv, ddf_inv, epsilon) {}

float RealLineAutomorphism::inv(float x) const {return inverse(x);}

RealLineAutomorphism RealLineAutomorphism::inv() const {return RealLineAutomorphism(inverse, *this);}

RealLineAutomorphism RealLineAutomorphism::operator~() const {return inv();}

RealLineAutomorphism RealLineAutomorphism::operator*(float c) const {return RealLineAutomorphism(static_cast<RealFunctionR1>(*this)*c, inverse/c);}

RealLineAutomorphism RealLineAutomorphism::operator+(float c) const {return RealLineAutomorphism(static_cast<RealFunctionR1>(*this)+c, inverse-c);}

RealLineAutomorphism RealLineAutomorphism::operator/(float c) const {return (*this)*(1/c);}

RealLineAutomorphism RealLineAutomorphism::operator-(float c) const {return (*this)+(-c);}

RealLineAutomorphism RealLineAutomorphism::operator-() const {return RealLineAutomorphism(-static_cast<RealFunctionR1>(*this), -inverse);}

RealLineAutomorphism operator*(float c, const RealLineAutomorphism &f) {return f*c;}

RealLineAutomorphism operator+(float c, const RealLineAutomorphism &f) {return f+c;}

RealLineAutomorphism operator-(float c, const RealLineAutomorphism &f) {return -f+c;}

RealLineAutomorphism RealLineAutomorphism::Id() {return RealLineAutomorphism(RealFunctionR1::x(), RealFunctionR1::x());}

RealLineAutomorphism RealLineAutomorphism::x() {return Id();}

vec2 PlaneSmoothEndomorphism::operator()(float t, float u) const { return (*this)(vec2(t, u)); }

vec2 PlaneSmoothEndomorphism::df(float t, float u, vec2 v) const { return df(vec2(t, u), v); }

mat2 PlaneSmoothEndomorphism::df(float t, float u) const { return df(vec2(t, u)); }

vec2 PlaneAutomorphism::inv(vec2 x) const { return f_inv(x); }

vec2 PlaneAutomorphism::inv(float t, float u) const { return f_inv(vec2(t, u)); }

PlaneAutomorphism PlaneAutomorphism::inv() const { return ~(*this); }


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

RealFunctionR3 RealFunctionR1::operator&(const RealFunctionR3 &g_) const { return RealFunctionR3([f=*this, g=g_](vec3 v) { return f(g(v)); }, eps); }

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
