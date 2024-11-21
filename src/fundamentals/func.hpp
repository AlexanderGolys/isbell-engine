#pragma once

#include "mat.hpp"
#include <functional>
#include <memory>
#include <iostream>
#include <optional>
#include <utility>
#include <variant>
#include <string>






inline std::string polyGroupIDtoString(PolyGroupID id) { return std::holds_alternative<int>(id) ? std::to_string(std::get<int>(id)) : std::get<std::string>(id); }
inline PolyGroupID prefix(const PolyGroupID &id, std::string prefix) { return PolyGroupID(prefix + polyGroupIDtoString(id)); }
inline std::string randomString() { int a = rand(); return std::to_string(a); }
inline PolyGroupID make_unique_id(const PolyGroupID &id) { return prefix(id, randomString()); }
inline PolyGroupID bdGroup(int n) { return PolyGroupID("bd" + std::to_string(n));}
inline PolyGroupID curveGroup(int n) { return PolyGroupID("curva" + std::to_string(n));}
inline PolyGroupID randomCurvaID() { return PolyGroupID("curva" + randomString());}
inline PolyGroupID randomBdID() { return PolyGroupID("bd" + randomString());}
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



template<typename domain, typename codomain=float>
class Morphism {
protected:
    std::function<codomain(domain)> _f;
public:
    Morphism(std::function<codomain(domain)> f) : _f(f) {} // NOLINT(*-explicit-constructor)
    codomain operator()(domain x) const { return _f(x); }
    Morphism operator+(const Morphism &g) const requires AbelianSemigroup<codomain> { return Morphism([f_=_f, g_=g._f](domain x) { return f_(x) + g_(x); }); }
    // friend Morphism operator+(codomain a, const Morphism &f) requires AbelianSemigroup<codomain> { return Morphism([f_ = _f, a](domain x) { return f_(x) + a; }); }

     Morphism operator-(const Morphism &g) const requires AbelianGroupConcept<codomain> { return Morphism([this, g](domain x) { return (*this)(x) - g(x); }); }
     Morphism operator*(const Morphism &g) const requires GroupConcept<codomain> { return Morphism([this, g](domain x) { return (*this)(x) * g(x); }); }
     Morphism operator/(const Morphism &g) const requires DivisionRing<codomain> { return Morphism([this, g](domain x) { return (*this)(x) / g(x); }); }
     Morphism operator*(codomain a) const requires GroupConcept<codomain> { return Morphism([this, a](domain x) { return (*this)(x) * a; }); }
    // friend Morphism operator*(codomain a, const Morphism &f) requires GroupConcept<codomain> { return Morphism([f, a](domain x) { return a*f(x); }); }

    template<DivisionRing K>  Morphism operator/(K a) const requires VectorSpaceConcept<codomain, K>      { return Morphism([this, a](domain x) { return (*this)(x) / a; }); }
    template<Rng R>           Morphism operator*(R a) const requires ModuleConcept<codomain, R>           { return Morphism([this, a](domain x) { return (*this)(x) * a; }); }
};

template<typename domain, typename codomain>
class Isomorphism : public Morphism<domain, codomain> {
    Morphism<codomain, domain> _inverse;
public:
    Isomorphism(std::function<codomain(domain)> f, std::function<domain(codomain)> g) : Morphism<domain, codomain>(f), _inverse(g) {}
    Isomorphism inverseMorphism() const { return Isomorphism(_inverse, this->_f); }
    Isomorphism operator~() const { return inverseMorphism(); }
    domain inv(codomain y) const { return _inverse(y); }
};

template<typename domain>
class Endomorphism : public Morphism<domain, domain> {
public:
    using Morphism<domain, domain>::Morphism;
    static Endomorphism id() { return Endomorphism([](domain x) { return x; }); }
    Endomorphism pow(int p) const {
        if (p < 0) return ~(*this).pow(-p);
        if (p == 0) return Endomorphism::id();
        if (p == 1) return *this;
        if (p % 2 == 0) return this->pow(p / 2) * this->pow(p / 2);
        return (*this) * this->pow(p / 2) * this->pow(p / 2); }
    Endomorphism operator^(int p) const { return pow(p); }
};



// composition, used with operator f&g
template <typename X, typename Y, typename Z>
Morphism<X, Z> compose(const Morphism<Y, Z> &f, const Morphism<X, Y> &g) {
    return Morphism<X, Z>([f_=f, g_=g](X x) { return f_(g_(x)); });
}
template <typename X, typename Y, typename Z>
Morphism<X, Z> operator&(const Morphism<Y, Z> &f, const Morphism<X, Y> &g) {
    return Morphism<X, Z>([f_=f, g_=g](X x) { return f_(g_(x)); });
}

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

Foo31 partialDerivativeOperator(Foo31 f, int i, float epsilon=0.01);
Foo21 partialDerivativeOperator(Foo21 f, int i, float epsilon=0.01);

Foo33 derivativeOperator(const Foo31 &f, float epsilon=0.01);
Foo22 derivativeOperator(const Foo21 &f, float epsilon=0.01);

Foo13 derivativeOperator(const Foo13 &f, float epsilon=0.01);
Foo12 derivativeOperator(const Foo12 &f, float epsilon=0.01);

class VectorFieldR2 {
public:
	Foo22 field;
	VectorFieldR2() : field([](glm::vec2 v) { return glm::vec2(0, 0); }) {}
	explicit VectorFieldR2(const Foo22 &field) : field(field) {}
	glm::vec2 operator()(glm::vec2 v) const { return field(v); }
};


class RealFunctionR3 {
	Foo31 _f;
	Foo33 _df;
    float eps = 0.01;
    Regularity regularity;
public:
	RealFunctionR3() : _f([](glm::vec3 v) { return 0; }), _df([](glm::vec3 v) { return glm::vec3(0); }), regularity(Regularity::ANALYTIC) {}
    RealFunctionR3(const RealFunctionR3 &other) = default;
    RealFunctionR3(RealFunctionR3 &&other) noexcept;
    RealFunctionR3 & operator=(const RealFunctionR3 &other);
    RealFunctionR3 & operator=(RealFunctionR3 &&other) noexcept;

    RealFunctionR3(Foo31 f,Foo33 df, float eps=.01, Regularity regularity = Regularity::SMOOTH)
    : _f(std::move(f)), _df(std::move(df)), eps(eps), regularity(regularity){};

	explicit RealFunctionR3(Foo31 f, float epsilon=0.01, Regularity regularity = Regularity::SMOOTH)
    :  _f(std::move(f)), _df(derivativeOperator(f, eps)), eps(eps), regularity(regularity){};

	float operator()(glm::vec3 v) const;
	glm::vec3 df(glm::vec3 v) const;

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

    VectorFieldR3 gradient() const;
    RealFunctionR3 Laplacian() const;
    float dx (glm::vec3 x) const { return _df(x).x; }
    float dy (glm::vec3 x) const { return _df(x).y; }
    float dz (glm::vec3 x) const { return _df(x).z; }

    static RealFunctionR3 linear(glm::vec3 v);
    static RealFunctionR3 constant(float a);
    static RealFunctionR3 projection(int i);
};


// R3 -> R3
class SpaceEndomorphism {
protected:
	Foo33 _f;
	Foo3Foo33 _df;
    float eps = 0.01;
public:
  SpaceEndomorphism(const SpaceEndomorphism &other) : _f(other._f), _df(other._df){}
  SpaceEndomorphism(SpaceEndomorphism &&other) noexcept : _f(std::move(other._f)), _df(std::move(other._df)) {}
  SpaceEndomorphism(Foo33 f, Foo3Foo33 df, float eps=.01) : _f(std::move(f)), _df(std::move(df)), eps(eps) {};
  explicit SpaceEndomorphism(Foo33 f, float epsilon=0.01);
  SpaceEndomorphism &operator=(const SpaceEndomorphism &other);
  SpaceEndomorphism &operator=(SpaceEndomorphism &&other) noexcept;
  explicit SpaceEndomorphism(glm::mat3 A) : _f([A](glm::vec3 x) { return A * x; }), _df([A](glm::vec3 x) { return A; }) {}
  explicit SpaceEndomorphism(glm::mat4 A) : _f([A](glm::vec3 x) { return glm::vec3(A * glm::vec4(x, 1)); }), _df([A](glm::vec3 x) { return glm::mat3(A); }) {}

	glm::vec3 directional_derivative(glm::vec3 x, glm::vec3 v) const { return _df(x) * v; }
    glm::vec3 dfdv(glm::vec3 x, glm::vec3 v) const { return directional_derivative(x, v); }
	glm::mat3 df(glm::vec3 x) const { return _df(x); }
	glm::vec3 operator()(glm::vec3 x) const { return _f(x); }
	virtual SpaceEndomorphism compose(const SpaceEndomorphism &g) const;
    virtual SpaceEndomorphism operator&(const SpaceEndomorphism &g) const { return compose(g); }

	static SpaceEndomorphism linear(glm::mat3 A) { return SpaceEndomorphism(A); }
	static SpaceEndomorphism translation(glm::vec3 v);
	static SpaceEndomorphism scaling(float x, float y, float z);
	static SpaceEndomorphism affine(glm::mat3 A, glm::vec3 v);
};

// R3 -> R3
class VectorFieldR3 {
	Foo33 _X;
    Foo3Foo33 _dX;
    float eps = 0.01;
public:
	VectorFieldR3();
    VectorFieldR3(Foo33 X, Foo3Foo33 dX, float eps=.01) : _X(std::move(X)), _dX(std::move(dX)), eps(eps) {}
    explicit VectorFieldR3(Foo33 X, float eps=.01);
    VectorFieldR3(RealFunctionR3 Fx , RealFunctionR3 Fy, RealFunctionR3 Fz, float epsilon=0.01);
	explicit VectorFieldR3(Foo33 field);
	explicit VectorFieldR3(VectorFieldR2 f);
	VectorFieldR3 operator+(const VectorFieldR3 &Y) const;
	VectorFieldR3 operator*(float a) const;
    VectorFieldR3 operator-() const { return *this * -1; }
	VectorFieldR3 operator-(const VectorFieldR3 &Y) const { return *this + (-Y); }
	VectorFieldR3 operator*(const RealFunctionR3 &f) const;
    RealFunctionR3 F_x() const { return RealFunctionR3([this](glm::vec3 x) { return _X(x).x; }, [this](glm::vec3 x) { return _dX(x)[0]; }); }
    RealFunctionR3 F_y() const { return RealFunctionR3([this](glm::vec3 x) { return _X(x).y; }, [this](glm::vec3 x) { return _dX(x)[1]; }); }
    RealFunctionR3 F_z() const { return RealFunctionR3([this](glm::vec3 x) { return _X(x).z; }, [this](glm::vec3 x) { return _dX(x)[2]; }); }
    std::array<RealFunctionR3, 3> components() const { return {F_x(), F_y(), F_z()}; }
    R3 operator()(R3 v) const { return _X(v); }

    friend VectorFieldR3 operator*(const glm::mat3 &A, const VectorFieldR3 &X) {
      return VectorFieldR3([f=X._X, A](glm::vec3 v) {return A * f(v); }, [df=X._dX, A](glm::vec3 v) {return A * df(v); }, X.eps);
    }

	static VectorFieldR3 constant(glm::vec3 v);
	static VectorFieldR3 linear(glm::mat3 A);
	static VectorFieldR3 radial(glm::vec3 scale);

    RealFunctionR3 divergence() const;
    VectorFieldR3 curl() const;

    glm::vec3 moveAlong(glm::vec3 v, float dt=1) const { return v + _X(v) * dt; }
};


inline VectorFieldR3 RealFunctionR3::gradient() const { return VectorFieldR3(_df, eps); }
inline RealFunctionR3 RealFunctionR3::Laplacian() const { return gradient().divergence(); }


// aut(R3)
class SpaceAutomorphism : public SpaceEndomorphism {
	Foo33 _f_inv;
public:
	SpaceAutomorphism(Foo33 f, Foo33 f_inv, std::function<glm::mat3(glm::vec3)> df);
	SpaceAutomorphism(Foo33 f, Foo33 f_inv, float epsilon=0.01);

	glm::vec3 inv(glm::vec3 v) const;
	SpaceAutomorphism operator~() const;
	SpaceAutomorphism inv() const { return ~(*this); }
	SpaceAutomorphism compose(SpaceAutomorphism g) const;
    SpaceAutomorphism operator&(const SpaceAutomorphism &g) const { return compose(g); }
    SpaceAutomorphism applyWithBasis(glm::mat3 A) const { return linear(A) & *this & linear(glm::inverse(A)); }
    SpaceAutomorphism applyWithBasis(glm::vec3 v1, glm::vec3 v2, glm::vec3 v3) const { return applyWithBasis(glm::mat3(v1, v2, v3)); }
    SpaceAutomorphism applyWithShift(glm::vec3 v) const { return translation(v) & *this & translation(-v); }

    static SpaceAutomorphism linear(glm::mat3 A);
    static SpaceAutomorphism translation(glm::vec3 v);
    static SpaceAutomorphism scaling(float x, float y, float z);
    static SpaceAutomorphism scaling(glm::vec3 factors) { return scaling(factors.x, factors.y, factors.z); }

    static SpaceAutomorphism scaling(float x) { return scaling(x, x, x); }
    static SpaceAutomorphism scaling(float x, float y, float z, glm::vec3 center) { return translation(center) & scaling(x, y, z) & translation(-center); }
    static SpaceAutomorphism scaling(float x, glm::vec3 center) { return scaling(x, x, x, center); }
    static SpaceAutomorphism scaling(glm::vec3 factors, glm::vec3 center) {  return scaling(factors.x, factors.y, factors.z, center); }

    static SpaceAutomorphism affine(glm::mat3 A, glm::vec3 v);
    static SpaceAutomorphism rotation(float angle);
    static SpaceAutomorphism rotation(glm::vec3 axis, float angle);
    static SpaceAutomorphism deltaRotation(glm::vec3 v1, glm::vec3 v2) { return linear(rotationBetween(v1, v2)); }
    static SpaceAutomorphism deltaRotation(glm::vec3 v1, glm::vec3 v2, glm::vec3 center) { return translation(center) & deltaRotation(v1, v2) & translation(-center); }

};



class PlaneSmoothEndomorphism {
protected:
	std::shared_ptr<Foo22> _f;
	std::shared_ptr<std::function<glm::mat2(glm::vec2)>> _df;
public:
	PlaneSmoothEndomorphism();
	PlaneSmoothEndomorphism(std::shared_ptr<Foo22> f, std::shared_ptr<std::function<glm::mat2(glm::vec2)>> _df);
	glm::vec2 operator()(glm::vec2 x) const;
	glm::vec2 df(glm::vec2 x, glm::vec2 v) const;
	glm::mat2 df(glm::vec2 x) const;
};

class PlaneAutomorphism : public PlaneSmoothEndomorphism {
protected:
	std::shared_ptr<Foo22> _f_inv;
public:
	PlaneAutomorphism(std::shared_ptr<Foo22> f, std::shared_ptr<std::function<glm::mat2(glm::vec2)>> _df, std::shared_ptr<Foo22> f_inv);
	glm::vec2 f_inv(glm::vec2 x) const;
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

template<typename Dom, typename Cod, int k>
std::array<Cod, k> mapMove(std::array<Dom, k> v, const HOM( Dom&&, Cod) &f) {
	std::array<Cod, k> res;
	for (int i = 0; i < k; i++) res[i] = f(v[i]);
	return res;
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

template<typename  A=float, typename container=vector<A>>
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
