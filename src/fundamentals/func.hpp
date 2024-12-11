#pragma once

#include <functional>
#include <iostream>
#include <memory>
#include <optional>
#include <string>
#include <utility>
#include <variant>

//#include <src/geometry/smoothParametric.hpp>

#include "mat.hpp"



std::string randomString();
int randomInt();
long randomLong();
float randomFloat(float a=0, float b=1);


inline std::string polyGroupIDtoString(PolyGroupID id) { return std::holds_alternative<int>(id) ? std::to_string(std::get<int>(id)) : std::get<std::string>(id); }
inline PolyGroupID prefix(const PolyGroupID &id, const std::string& prefix) { return PolyGroupID(prefix + polyGroupIDtoString(id)); }

inline PolyGroupID make_unique_id(const PolyGroupID &id) { return prefix(id, randomString()); }
inline PolyGroupID bdGroup(int n) { return PolyGroupID("bd" + std::to_string(n));}
inline PolyGroupID curveGroup(int n) { return PolyGroupID("curva" + std::to_string(n));}
inline PolyGroupID randomCurvaID() { return PolyGroupID("curva" + randomString());}
inline PolyGroupID randomID() { return PolyGroupID(randomString());}

inline PolyGroupID randomBdID() { return PolyGroupID("bd" + randomString()+ randomString()+ randomString());}
inline bool idBelongsToCurve(const PolyGroupID &id) { return std::holds_alternative<std::string>(id) && std::get<std::string>(id).starts_with("curva");}
inline bool idBelongsToBd(const PolyGroupID &id) { return std::holds_alternative<std::string>(id) && std::get<std::string>(id).starts_with("bd");}


class VectorFieldR3;





enum class Regularity {
	PATHOLOGICAL=-99,
    UNKNOWN = -2,
	DISCONTINUOUS=-1,
	C0 = 0,
	C1 = 1,
	C2 = 2,
	C3 = 3,
	C4 = 4,
    C5 = 5,
	POSSIBLY_FINITE_REGULARITY_AT_LEAST_C6 = 6,
	SMOOTH = 69,
    ANALYTIC = 70,
	CONFORMAL= 71,
	MEROMORPHIC = 420,
	HOLOMORPHIC = 2137
};

inline Regularity max(Regularity a, Regularity b) {return static_cast<Regularity>(std::max(INT(a), INT(b))); }
inline Regularity min(Regularity a, Regularity b) {return static_cast<Regularity>(std::min(INT(a), INT(b))); }
Regularity operator-(Regularity a, int b);
Regularity operator+(Regularity a, int b);
inline Regularity operator--(Regularity a) { return a - 1; }
inline Regularity operator++(Regularity a) { return a + 1; }


// universal morphism of product
// template <typename X, typename Y>
// prod(X, Y) x(X x, Y y) {
//     return prod(X, Y)({x, y});
// }

#define Hom(X, Y) Morphism<X, Y>
#define Iso(X, Y) Isomorphism<X, Y>
#define End(X) Endomorphism<X>

// template <typename X, typename Y>
// Hom(prod(X, Y), X) pi1 = Morphism([](prod(X, Y) xy) { return xy.first; } );
//
// template <typename X, typename Y>
// Hom(prod(X, Y), Y) pi2 = Morphism([](prod(X, Y) xy) { return xy.second; } );
//
// template <typename X, typename Y, typename Z>
// std::function<Z(X, Y)> f2(Hom(prod(X, Y), Z) f) {
//     return [f](X x, Y y) { return f(p(x, y)); };
// }

// template <typename X, typename Y, typename A>
// auto productUP = auto ([](prod(Hom(A, X), Hom(A, Y)) fg) { return auto ([fg](A a) { return x(fg&pi1(a), fg&pi2(a)); }); });


// template <typename X, typename Y, typename Z>
// Hom(X, Hom(Y, Z)) curry(const Morphism<prod(X, Y), Z> &f) {
//     return auto ([f_=f](X x) { return auto ([f_, x](Y y) { return f_(x, y); }); });
// }
//
// template <typename X, typename Y, typename Z>
// Morphism<prod(X, Y), Z> uncurry(const Morphism<X, Morphism<Y, Z>> &f) {
//     return auto ([f_=f](prod(X, Y) xy) { return f_(xy.first)(xy.second); });
// }




template<typename codomain>
std::function<codomain(float)> derivativeOperator(std::function<codomain(float)> f, float epsilon) {
    return [f, epsilon](float x) {return (f(x + epsilon) - f(x - epsilon)) / (2 * epsilon); };
}

template<typename  domain, typename  codomain>
std::function<codomain(domain)> directionalDerivativeOperator(std::function<codomain(domain)> f, domain v, float epsilon) {
    return [f, epsilon, v](domain x) { return (f(x + v * epsilon) - f(x - v * epsilon)) / (2 * epsilon); };
}

Foo31 partialDerivativeOperator(Foo31 f, int i, float epsilon);
Foo21 partialDerivativeOperator(Foo21 f, int i, float epsilon);

Foo33 derivativeOperator(const Foo31 &f, float epsilon);
Foo22 derivativeOperator(const Foo21 &f, float epsilon);

Foo13 derivativeOperator(const Foo13 &f, float epsilon);
Foo12 derivativeOperator(const Foo12 &f, float epsilon);

template<typename Y>
HOM(vec3, Y) wrap(const TRIHOM(float, float, float, Y) &f_) { return [f=f_](vec3 v) { return f(v.x, v.y, v.z); }; }
template<typename Y>
TRIHOM(float, float, float, Y) unwrap(const HOM(vec3, Y) &f_) { return [f=f_](float x, float y, float z) { return f(vec3(x, y, z)); }; }

class VectorFieldR2 {
public:
	Foo22 field;
	VectorFieldR2() : field([](vec2 v) { return vec2(0, 0); }) {}
	explicit VectorFieldR2(const Foo22 &field) : field(field) {}
	vec2 operator()(vec2 v) const { return field(v); }
	vec2 normal(vec2 x) const { return normalise(vec2(-field(x).y, field(x).x)); }
	float speed(vec2 x) const { return length(field(x)); }
};


class RealFunctionR3 {
	Foo31 _f;
	Foo33 _df;
    float eps = 0.01;
    Regularity regularity;
public:
	RealFunctionR3() : _f([](vec3 v) { return 0; }), _df([](vec3 v) { return vec3(0); }), regularity(Regularity::ANALYTIC) {}
    RealFunctionR3(const RealFunctionR3 &other) = default;
    RealFunctionR3(RealFunctionR3 &&other) noexcept;
    RealFunctionR3 & operator=(const RealFunctionR3 &other);
    RealFunctionR3 & operator=(RealFunctionR3 &&other) noexcept;

    RealFunctionR3(Foo31 f,Foo33 df, float eps=.01, Regularity regularity = Regularity::SMOOTH)
    : _f(f), _df(df), eps(eps), regularity(regularity){};

	explicit RealFunctionR3(Foo31 f, float epsilon=0.01, Regularity regularity = Regularity::SMOOTH)
    :  _f(f), _df(derivativeOperator(f, epsilon)), eps(epsilon), regularity(regularity){};

	RealFunctionR3(Foo1111 f,Foo1113 df, float eps=.01, Regularity regularity = Regularity::SMOOTH)
	: _f(wrap(f)), _df(wrap(df)), eps(eps), regularity(regularity){};

	explicit RealFunctionR3(Foo1111 f, float epsilon=0.01, Regularity regularity = Regularity::SMOOTH)
	:  _f(wrap(f)), _df(derivativeOperator(wrap(f), epsilon)), eps(epsilon), regularity(regularity){};

	float operator()(vec3 v) const;
	vec3 df(vec3 v) const;
	float operator()(float x, float y, float z) const { return _f(vec3(x, y, z)); }
	vec3 df(float x, float y, float z) const { return _df(vec3(x, y, z)); }

    RealFunctionR3 operator*(float a) const;
    RealFunctionR3 operator+(float a) const;
    RealFunctionR3 operator-(float a) const;
    RealFunctionR3 operator/(float a) const;
    RealFunctionR3 operator-() const { return *this * -1; }
	RealFunctionR3 operator+(const RealFunctionR3 &g) const;
	RealFunctionR3 operator-(const RealFunctionR3 &g) const { return *this + (-g); }
	RealFunctionR3 operator*(const RealFunctionR3 &g) const;
	RealFunctionR3 operator/(const RealFunctionR3 &g) const;
    friend RealFunctionR3 operator*(float a, const RealFunctionR3 &f) { return f * a; }
    friend RealFunctionR3 operator+(float a, const RealFunctionR3 &f) { return f + a; }
    friend RealFunctionR3 operator-(float a, const RealFunctionR3 &f) { return -f + a; }
    friend RealFunctionR3 operator/(float a, const RealFunctionR3 &f) { return constant(a)/f; }
    RealFunctionR3 operator~() const { return 1/(*this); }
	float getEps() const { return eps; }

    VectorFieldR3 gradient() const;
    RealFunctionR3 Laplacian() const;
    float dx (vec3 x) const { return _df(x).x; }
    float dy (vec3 x) const { return _df(x).y; }
    float dz (vec3 x) const { return _df(x).z; }

    static RealFunctionR3 linear(vec3 v);
    static RealFunctionR3 constant(float a);
    static RealFunctionR3 projection(int i);
};


// R3 -> R3
class SpaceEndomorphism {
protected:
	Foo33 _f;
	Foo3Foo33 _df;
    float eps = 0.01;
	virtual SpaceEndomorphism compose(const SpaceEndomorphism &g) const;


public:
	virtual ~SpaceEndomorphism() = default;

	SpaceEndomorphism(const SpaceEndomorphism &other) : _f(other._f), _df(other._df){}
  SpaceEndomorphism(SpaceEndomorphism &&other) noexcept : _f(std::move(other._f)), _df(std::move(other._df)) {}
  SpaceEndomorphism(Foo33 f, Foo3Foo33 df, float eps=.01) : _f(std::move(f)), _df(std::move(df)), eps(eps) {};
  explicit SpaceEndomorphism(Foo33 f, float epsilon=0.01);
  SpaceEndomorphism &operator=(const SpaceEndomorphism &other);
  SpaceEndomorphism &operator=(SpaceEndomorphism &&other) noexcept;
  explicit SpaceEndomorphism(mat3 A) : _f([A](vec3 x) { return A * x; }), _df([A](vec3 x) { return A; }) {}
  explicit SpaceEndomorphism(mat4 A) : _f([A](vec3 x) { return vec3(A * vec4(x, 1)); }), _df([A](vec3 x) { return mat3(A); }) {}

	vec3 directional_derivative(vec3 x, vec3 v) const { return _df(x) * v; }
    vec3 dfdv(vec3 x, vec3 v) const { return directional_derivative(x, v); }
	mat3 df(vec3 x) const { return _df(x); }
	vec3 operator()(vec3 x) const { return _f(x); }
    virtual SpaceEndomorphism operator&(const SpaceEndomorphism &g) const { return compose(g); }

	float getEps() const { return eps; }

	static SpaceEndomorphism linear(const mat3 &A) { return SpaceEndomorphism(A); }
	static SpaceEndomorphism translation(vec3 v);
	static SpaceEndomorphism scaling(float x, float y, float z);
	static SpaceEndomorphism affine(mat3 A, vec3 v);
};


// R3 -> R3
class PlaneEndomorphism {
protected:
	Foo22 _f;
	Foo2Foo22 _df;
	float eps = 0.01;
	virtual PlaneEndomorphism compose(const PlaneEndomorphism &g) const;

public:
	virtual ~PlaneEndomorphism() = default;

	PlaneEndomorphism(const PlaneEndomorphism &other) : _f(other._f), _df(other._df) {}
	PlaneEndomorphism(PlaneEndomorphism &&other) noexcept : _f(std::move(other._f)), _df(std::move(other._df)) {}
	PlaneEndomorphism(Foo22 f, Foo2Foo22 df, float eps=.01) : _f(std::move(f)), _df(std::move(df)), eps(eps) {}
	explicit PlaneEndomorphism(Foo22 f, float epsilon=0.01);
	PlaneEndomorphism &operator=(const PlaneEndomorphism &other);
	PlaneEndomorphism &operator=(PlaneEndomorphism &&other) noexcept;
	explicit PlaneEndomorphism(mat2 A) : _f([A](vec2 x) { return A * x; }), _df([A](vec2 x) { return A; }) {}

	vec2 directional_derivative(vec2 x, vec2 v) const { return _df(x) * v; }
	vec2 dfdv(vec2 x, vec2 v) const { return directional_derivative(x, v); }
	mat2 df(vec2 x) const { return _df(x); }
	mat2 df(float x, float y) const { return _df(vec2(x, y)); }
	vec2 operator()(vec2 x) const { return _f(x); }
	vec2 operator()(float x, float y) const { return _f(vec2(x, y)); }
	virtual PlaneEndomorphism operator&(const PlaneEndomorphism &g) const { return compose(g); }

	static PlaneEndomorphism linear(mat2 A) { return PlaneEndomorphism(A); }
	static PlaneEndomorphism translation(vec2 v);
	static PlaneEndomorphism scaling(float x, float y);
	static PlaneEndomorphism affine(mat2 A, vec2 v);
};

// R3 -> R3
class VectorFieldR3 {
	Foo33 _X;
    Foo3Foo33 _dX;
    float eps = 0.01;
public:
	VectorFieldR3();
    VectorFieldR3(Foo33 X, Foo3Foo33 dX, float eps=.01) : _X(X), _dX(std::move(dX)), eps(eps) {}
    explicit VectorFieldR3(Foo33 X, float eps=.01);
    VectorFieldR3(RealFunctionR3 Fx , RealFunctionR3 Fy, RealFunctionR3 Fz, float epsilon=0.01);
	explicit VectorFieldR3(Foo33 field);
	explicit VectorFieldR3(VectorFieldR2 f);
	VectorFieldR3 operator+(const VectorFieldR3 &Y) const;
	VectorFieldR3 operator*(float a) const;
    VectorFieldR3 operator-() const { return *this * -1; }
	VectorFieldR3 operator-(const VectorFieldR3 &Y) const { return *this + (-Y); }
	VectorFieldR3 operator*(const RealFunctionR3 &f) const;
    RealFunctionR3 F_x() const { return RealFunctionR3([this](vec3 x) { return _X(x).x; }, [this](vec3 x) { return _dX(x)[0]; }); }
    RealFunctionR3 F_y() const { return RealFunctionR3([this](vec3 x) { return _X(x).y; }, [this](vec3 x) { return _dX(x)[1]; }); }
    RealFunctionR3 F_z() const { return RealFunctionR3([this](vec3 x) { return _X(x).z; }, [this](vec3 x) { return _dX(x)[2]; }); }
    std::array<RealFunctionR3, 3> components() const { return {F_x(), F_y(), F_z()}; }
    R3 operator()(R3 v) const { return _X(v); }

    friend VectorFieldR3 operator*(const mat3 &A, const VectorFieldR3 &X) {
      return VectorFieldR3([f=X._X, A](vec3 v) {return A * f(v); }, [df=X._dX, A](vec3 v) {return A * df(v); }, X.eps);
    }

	static VectorFieldR3 constant(vec3 v);
	static VectorFieldR3 linear(mat3 A);
	static VectorFieldR3 radial(vec3 scale);

    RealFunctionR3 divergence() const;
    VectorFieldR3 curl() const;

    vec3 moveAlong(vec3 v, float dt=1) const { return v + _X(v) * dt; }
};

class VectorFieldR3Dynamic {
	BIHOM(float, vec3, vec3) _X;
	float eps = 0.01;
	public:
	explicit VectorFieldR3Dynamic(BIHOM(float, vec3, vec3) X, float epsilon=0.01) : _X(X), eps(epsilon) {}
	vec3 operator()(float t, vec3 x) const { return _X(t, x); }
	VectorFieldR3Dynamic operator+(const VectorFieldR3Dynamic &Y) const { return VectorFieldR3Dynamic([__X=_X, _Y=Y._X](float t, vec3 x) { return  __X(t, x) + _Y(t, x); }, eps); }
	VectorFieldR3Dynamic operator*(float a) const { return VectorFieldR3Dynamic([__X=_X, a](float t, vec3 x) { return  __X(t, x) * a; }, eps); }
};


inline VectorFieldR3 RealFunctionR3::gradient() const { return VectorFieldR3(_df, eps); }
inline RealFunctionR3 RealFunctionR3::Laplacian() const { return gradient().divergence(); }


// aut(R3)
class SpaceAutomorphism : public SpaceEndomorphism {
	Foo33 _f_inv;
public:
	SpaceAutomorphism(Foo33 f, Foo33 f_inv, std::function<mat3(vec3)> df);
	SpaceAutomorphism(Foo33 f, Foo33 f_inv, float epsilon=0.01);

	vec3 inv(vec3 v) const;
	SpaceAutomorphism operator~() const;
	SpaceAutomorphism inv() const { return ~(*this); }
	SpaceAutomorphism compose(SpaceAutomorphism g) const;
    SpaceAutomorphism operator&(const SpaceAutomorphism &g) const { return compose(g); }
    SpaceAutomorphism applyWithBasis(mat3 A) const { return linear(A) & *this & linear(inverse(A)); }
    SpaceAutomorphism applyWithBasis(vec3 v1, vec3 v2, vec3 v3) const { return applyWithBasis(mat3(v1, v2, v3)); }
    SpaceAutomorphism applyWithShift(vec3 v) const { return translation(v) & *this & translation(-v); }

    static SpaceAutomorphism linear(mat3 A);
    static SpaceAutomorphism translation(vec3 v);
    static SpaceAutomorphism scaling(float x, float y, float z);
    static SpaceAutomorphism scaling(vec3 factors) { return scaling(factors.x, factors.y, factors.z); }

    static SpaceAutomorphism scaling(float x) { return scaling(x, x, x); }
    static SpaceAutomorphism scaling(float x, float y, float z, vec3 center) { return translation(center) & scaling(x, y, z) & translation(-center); }
    static SpaceAutomorphism scaling(float x, vec3 center) { return scaling(x, x, x, center); }
    static SpaceAutomorphism scaling(vec3 factors, vec3 center) {  return scaling(factors.x, factors.y, factors.z, center); }

    static SpaceAutomorphism affine(mat3 A, vec3 v);
    static SpaceAutomorphism rotation(float angle);
    static SpaceAutomorphism rotation(vec3 axis, float angle);
	static SpaceAutomorphism rotation(vec3 axis, float angle, vec3 center) { return rotation(axis, angle).applyWithShift(center); }
    static SpaceAutomorphism deltaRotation(vec3 v1, vec3 v2) { return linear(rotationBetween(v1, v2)); }
    static SpaceAutomorphism deltaRotation(vec3 v1, vec3 v2, vec3 center) { return translation(center) & deltaRotation(v1, v2) & translation(-center); }

};

class RealFunctionR1 {
	Fooo _f;
	Fooo _df;
	Fooo _ddf;
	float eps = 0.01;
public:
	RealFunctionR1(Fooo f, Fooo df, Fooo ddf, float epsilon=0.01) : _f(f), _df(df), _ddf(ddf), eps(epsilon) {}
	RealFunctionR1(Fooo f, Fooo df, float epsilon=0.01) :  RealFunctionR1(f, df, derivativeOperator(df, epsilon), epsilon) {}
	explicit RealFunctionR1(Fooo f, float epsilon=0.01) : RealFunctionR1(f, derivativeOperator(f, epsilon), epsilon) {}

	RealFunctionR1(const RealFunctionR1 &other);
	RealFunctionR1(RealFunctionR1 &&other) noexcept;
	RealFunctionR1 & operator=(const RealFunctionR1 &other);
	RealFunctionR1 & operator=(RealFunctionR1 &&other) noexcept;

	float operator()(float x) const { return _f(x); }
	float df(float x) const { return _df(x); }
	float ddf(float x) const { return _ddf(x); }
	RealFunctionR1 df() const {return RealFunctionR1(_df, _ddf, eps);}

	RealFunctionR1 operator+(const RealFunctionR1 &g) const;
	RealFunctionR1 operator*(float a) const;
	RealFunctionR1 operator+(float a) const;
	RealFunctionR1 operator-(float a) const { return *this + (-a); }
	RealFunctionR1 operator*(const RealFunctionR1 &g_) const;
	RealFunctionR1 operator-() const { return *this * -1; }
	RealFunctionR1 operator-(const RealFunctionR1 &g) const { return *this + (-g); }
	RealFunctionR1 operator/(const RealFunctionR1 &g_) const;
	RealFunctionR1 operator/(float a) const { return *this * (1/a); }
	friend RealFunctionR1 operator*(float a, const RealFunctionR1 &f) { return f * a; }
	friend RealFunctionR1 operator+(float a, const RealFunctionR1 &f) { return f + a; }
	friend RealFunctionR1 operator-(float a, const RealFunctionR1 &f) { return -f + a; }
	friend RealFunctionR1 operator/(float a, const RealFunctionR1 &f) { return constant(a)/f; }

	RealFunctionR1 operator & (const RealFunctionR1 &g_) const;

	static RealFunctionR1 constant(float a) { return RealFunctionR1([a](float x) { return a; }, [a](float x) { return 0; }, [a](float x) { return 0; }); }
	static RealFunctionR1 x() { return RealFunctionR1([](float x) { return x; }, [](float x) { return 1; }, [](float x) { return 0; }); }
	static RealFunctionR1 one() { return constant(1); }
	static RealFunctionR1 zero() { return constant(0); }
	static RealFunctionR1 linear(float a, float b) { return x() * a + b; }
	static RealFunctionR1 quadratic(float a, float b, float c) { return x() * x() * a + x() * b + c; }
	static RealFunctionR1 monomial(int n);
	static RealFunctionR1 polynomial(std::vector<float> coeffs);
	static RealFunctionR1 sin() { return RealFunctionR1([](float x) { return std::sin(x); }, [](float x) { return std::cos(x); }, [](float x) { return -std::sin(x); }); }
	static RealFunctionR1 cos() { return RealFunctionR1([](float x) { return std::cos(x); }, [](float x) { return -std::sin(x); }, [](float x) { return -std::cos(x); }); }
	static RealFunctionR1 exp() { return RealFunctionR1([](float x) { return std::exp(x); }, [](float x) { return std::exp(x); }, [](float x) { return std::exp(x); }); }
	static RealFunctionR1 log() { return RealFunctionR1([](float x) { return std::log(x); }, [](float x) { return 1/x; }, [](float x) { return -1/(x*x); }); }
	static RealFunctionR1 sqrt() { return RealFunctionR1([](float x) { return std::sqrt(x); }, [](float x) { return 1/(2*std::sqrt(x)); }, [](float x) { return -1/(4*x*std::sqrt(x)); }); }
	static RealFunctionR1 pow(float a) { return RealFunctionR1([a](float x) { return std::pow(x, a); }, [a](float x) { return a * std::pow(x, a - 1); }, [a](float x) { return a * (a - 1) * std::pow(x, a - 2); }); }
	static RealFunctionR1 pow(int a) { return RealFunctionR1([a](float x) { return std::pow(x, a); }, [a](float x) { return a * std::pow(x, a - 1); }, [a](float x) { return a * (a - 1) * std::pow(x, a - 2); }); }
};

class RealLineAutomorphism : public RealFunctionR1 {
	RealFunctionR1 inverse;
public:
	RealLineAutomorphism(RealFunctionR1 f, RealFunctionR1 inv) : RealFunctionR1(f), inverse(inv) {}
	RealLineAutomorphism(Fooo f, Fooo f_inv, float epsilon=0.01) : RealFunctionR1(f, epsilon), inverse(f_inv, epsilon) {}
	RealLineAutomorphism(Fooo f, Fooo df, Fooo f_inv, float epsilon=0.01) : RealFunctionR1(f, df, epsilon), inverse(f_inv, epsilon) {}
	RealLineAutomorphism(Fooo f, Fooo df, Fooo f_inv, Fooo df_inv, float epsilon=0.01) : RealFunctionR1(f, df, epsilon), inverse(f_inv, df_inv, epsilon) {}
	RealLineAutomorphism(Fooo f, Fooo df, Fooo ddf, Fooo f_inv, Fooo df_inv, Fooo ddf_inv, float epsilon=0.01) : RealFunctionR1(f, df, ddf, epsilon), inverse(f_inv, df_inv, ddf_inv, epsilon) {}

	float inv(float x) const {return inverse(x);}
	RealLineAutomorphism inv() const {return RealLineAutomorphism(inverse, *this);}
	RealLineAutomorphism operator~() const {return inv();}

	RealLineAutomorphism operator& (const RealLineAutomorphism &g) const;
	RealLineAutomorphism operator*(float c) const {return RealLineAutomorphism(static_cast<RealFunctionR1>(*this)*c, inverse/c);}
	RealLineAutomorphism operator+(float c) const {return RealLineAutomorphism(static_cast<RealFunctionR1>(*this)+c, inverse-c);}
	RealLineAutomorphism operator/(float c) const {return (*this)*(1/c);}
	RealLineAutomorphism operator-(float c) const {return (*this)+(-c);}
	RealLineAutomorphism operator-() const {return RealLineAutomorphism(-static_cast<RealFunctionR1>(*this), -inverse);}
	friend RealLineAutomorphism operator*(float c, const RealLineAutomorphism &f) {return f*c;}
	friend RealLineAutomorphism operator+(float c, const RealLineAutomorphism &f) {return f+c;}
	friend RealLineAutomorphism operator-(float c, const RealLineAutomorphism &f) {return -f+c;}

	static RealLineAutomorphism Id() {return RealLineAutomorphism(RealFunctionR1::x(), RealFunctionR1::x());}
	static RealLineAutomorphism x() {return Id();}
	static RealLineAutomorphism pow(int n);
};




class PlaneSmoothEndomorphism {
protected:
	std::shared_ptr<Foo22> _f;
	std::shared_ptr<std::function<mat2(vec2)>> _df;
public:
	PlaneSmoothEndomorphism();
	PlaneSmoothEndomorphism(std::shared_ptr<Foo22> f, std::shared_ptr<std::function<mat2(vec2)>> _df);
	vec2 operator()(vec2 x) const;
	vec2 operator()(float t, float u) const { return (*this)(vec2(t, u)); }
	vec2 df(vec2 x, vec2 v) const;
	mat2 df(vec2 x) const;
	vec2 df(float t, float u, vec2 v) const { return df(vec2(t, u), v); }
	mat2 df(float t, float u) const { return df(vec2(t, u)); }
};

class PlaneAutomorphism : public PlaneSmoothEndomorphism {
protected:
	std::shared_ptr<Foo22> _f_inv;
public:
	PlaneAutomorphism(std::shared_ptr<Foo22> f, std::shared_ptr<std::function<mat2(vec2)>> _df, std::shared_ptr<Foo22> f_inv);
	vec2 f_inv(vec2 x) const;
	vec2 inv(vec2 x) const { return f_inv(x); }
	vec2 inv(float t, float u) const { return f_inv(vec2(t, u)); }
	PlaneAutomorphism operator~() const;
	PlaneAutomorphism inv() const { return ~(*this); }
};

template<typename Dom, typename Cod>
vector<Cod> map(vector<Dom> v, const HOM(const Dom&, Cod) &f) {
	vector<Cod> res;
	res.reserve(v.size());
	for (auto x : v) res.push_back(f(x));
	return res;
}




template<typename Dom, typename Cod, int k>
std::array<Cod, k> map(std::array<Dom, k> v, const HOM(const Dom&, Cod) &f) {
	std::array<Cod, k> res;
	for (int i = 0; i < k; i++) res[i] = f(v[i]);
	return res;
}

template<TotallyOrderedAbelianSemigroup Dom, typename Cod>
vector<Cod> mapArange(Dom a, Dom b, Dom step, const HOM(const Dom&, Cod) &f) {
	return map(arange(a, b, step), f);
}

template<typename Cod>
vector<Cod> mapRange(int a, int b, int step, const HOM(int, Cod) &f) {
	return map(range(a, b, step), f);
}

template<VectorSpaceConcept<float> Dom, typename Cod>
vector<Cod> mapLinspace(Dom a, Dom b, int steps, const HOM(const Dom&, Cod) &f) {
	return map(linspace(a, b, steps), f);
}




template<typename Dom, typename Cod, int k>
std::array<Cod, k> mapByVal(std::array<Dom, k> v, const HOM(Dom, Cod) &f) {
	std::array<Cod, k> res;
	for (int i = 0; i < k; i++) res[i] = f(v[i]);
	return res;
}


template<typename container, typename A>
A combine(container arr, BIHOM(A, A, A) binaryOperator) {
	A res = arr[0];
	for (int i = 1; i < arr.size(); i++) res = binaryOperator(res, arr[i]);
	return res;
}

template<typename  A=int, typename container=vector<A>>
A sum(container arr) {
	A res = arr[0];
	for (int i = 1; i < arr.size(); i++) res += arr[i];
	return res;
}

template<typename container, Semigroup A>
A mult(container arr) {return combine(arr, [](A a, A b) { return a * b; }); }


template<typename A, int k, typename homspace>
std::array<A, k> arrayComprehension(const homspace &f) {
	std::array<A, k> res = {f(0)};
	for (int i = 1; i < k; i++) res[i] = f(i);
	return res;
}

template<typename A, typename homspace>
vector<A> vectorComprehension(const homspace &f, int n) {
	vector res = {f(0)};
	res.reserve(n);
	for (int i = 1; i < n; i++) res[i] = f(i);
	return res;
}

inline vec3 stereoProjection(vec4 v) {
	return vec3(v.x / (1-v.w), v.y / (1-v.w), v.z / (1-v.w));
}


RealFunctionR1 smoothenReLu(float c0, float x_change_to_id);
RealFunctionR1 expImpulse(float peak);
RealFunctionR1 cubicShroom(float center, float margin);
RealFunctionR1 rationalInfiniteShroom(float steepFactor, float center);
RealFunctionR1 powerShroom(float begin, float end, float zeroExponent, float oneExponent);
RealFunctionR1 toneMap(float k);
