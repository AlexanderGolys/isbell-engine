#pragma once

#include <functional>
#include <iostream>
#include <memory>
#include <optional>
#include <string>
#include <utility>
#include <variant>


#include "filesUtils.hpp"


std::string polyGroupIDtoString(PolyGroupID id);
PolyGroupID prefix(const PolyGroupID &id, const std::string& prefix);

PolyGroupID make_unique_id(const PolyGroupID &id);
PolyGroupID bdGroup(int n);
PolyGroupID curveGroup(int n);
PolyGroupID randomCurvaID();
PolyGroupID randomID();

PolyGroupID randomBdID();
bool idBelongsToCurve(const PolyGroupID &id);
bool idBelongsToBd(const PolyGroupID &id);


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

Regularity max(Regularity a, Regularity b);
Regularity min(Regularity a, Regularity b);
Regularity operator-(Regularity a, int b);
Regularity operator+(Regularity a, int b);
Regularity operator--(Regularity a);
Regularity operator++(Regularity a);


template<typename codomain>
std::function<codomain(float)> derivativeOperator(std::function<codomain(float)> f, float epsilon);

template<typename  domain, typename  codomain>
std::function<codomain(domain)> directionalDerivativeOperator(std::function<codomain(domain)> f, domain v, float epsilon);

Foo31 partialDerivativeOperator(Foo31 f, int i, float epsilon);
Foo21 partialDerivativeOperator(Foo21 f, int i, float epsilon);

Foo33 derivativeOperator(const Foo31 &f, float epsilon);
Foo22 derivativeOperator(const Foo21 &f, float epsilon);

Foo13 derivativeOperator(const Foo13 &f, float epsilon);
Foo12 derivativeOperator(const Foo12 &f, float epsilon);




class VectorFieldR2 {
public:
	Foo22 field;
	VectorFieldR2();
	explicit VectorFieldR2(const Foo22 &field);
	vec2 operator()(vec2 v) const;
	vec2 normal(vec2 x) const;
	float speed(vec2 x) const;
};


class RealFunctionR3 {
	Foo31 _f;
	Foo33 _df;
    float eps = 0.01;
    Regularity regularity;
public:
	RealFunctionR3();
    RealFunctionR3(const RealFunctionR3 &other);
    RealFunctionR3(RealFunctionR3 &&other) noexcept;
    RealFunctionR3 & operator=(const RealFunctionR3 &other);
    RealFunctionR3 & operator=(RealFunctionR3 &&other) noexcept;

    RealFunctionR3(Foo31 f,Foo33 df, float eps=.01, Regularity regularity = Regularity::SMOOTH);;
	explicit RealFunctionR3(Foo31 f, float epsilon=0.01, Regularity regularity = Regularity::SMOOTH);;
	RealFunctionR3(Foo1111 f,Foo1113 df, float eps=.01, Regularity regularity = Regularity::SMOOTH);;
	explicit RealFunctionR3(Foo1111 f, float epsilon=0.01, Regularity regularity = Regularity::SMOOTH);;

	float operator()(vec3 v) const;
	vec3 df(vec3 v) const;
	float operator()(float x, float y, float z) const;
	vec3 df(float x, float y, float z) const;

	RealFunctionR3 operator*(float a) const;
    RealFunctionR3 operator+(float a) const;
    RealFunctionR3 operator-(float a) const;
    RealFunctionR3 operator/(float a) const;
    RealFunctionR3 operator-() const;
	RealFunctionR3 operator+(const RealFunctionR3 &g) const;
	RealFunctionR3 operator-(const RealFunctionR3 &g) const;
	RealFunctionR3 operator*(const RealFunctionR3 &g) const;
	RealFunctionR3 operator/(const RealFunctionR3 &g) const;
    friend RealFunctionR3 operator*(float a, const RealFunctionR3 &f);
	friend RealFunctionR3 operator+(float a, const RealFunctionR3 &f);
	friend RealFunctionR3 operator-(float a, const RealFunctionR3 &f);
	friend RealFunctionR3 operator/(float a, const RealFunctionR3 &f);
	RealFunctionR3 operator~() const;
	float getEps() const;

	VectorFieldR3 gradient() const;
    RealFunctionR3 Laplacian() const;
    float dx (vec3 x) const;
	float dy (vec3 x) const;
	float dz (vec3 x) const;

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
	virtual ~SpaceEndomorphism();

	SpaceEndomorphism(const SpaceEndomorphism &other);
  SpaceEndomorphism(SpaceEndomorphism &&other) noexcept;
  SpaceEndomorphism(Foo33 f, Foo3Foo33 df, float eps=.01);;
  explicit SpaceEndomorphism(Foo33 f, float epsilon=0.01);
  SpaceEndomorphism &operator=(const SpaceEndomorphism &other);
  SpaceEndomorphism &operator=(SpaceEndomorphism &&other) noexcept;
  explicit SpaceEndomorphism(mat3 A);
  explicit SpaceEndomorphism(mat4 A);

	vec3 directional_derivative(vec3 x, vec3 v) const;
	vec3 dfdv(vec3 x, vec3 v) const;
	mat3 df(vec3 x) const;
	vec3 operator()(vec3 x) const;
	virtual SpaceEndomorphism operator&(const SpaceEndomorphism &g) const;

	float getEps() const;

	static SpaceEndomorphism linear(const mat3 &A);
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
	virtual ~PlaneEndomorphism();

//	PlaneEndomorphism() = default;
	PlaneEndomorphism(const PlaneEndomorphism &other);
	PlaneEndomorphism(PlaneEndomorphism &&other) noexcept;
	PlaneEndomorphism(Foo22 f, Foo2Foo22 df, float eps=.01);
	explicit PlaneEndomorphism(Foo22 f, float epsilon=0.01);
	PlaneEndomorphism &operator=(const PlaneEndomorphism &other);
	PlaneEndomorphism &operator=(PlaneEndomorphism &&other) noexcept;
	explicit PlaneEndomorphism(mat2 A);

	vec2 directional_derivative(vec2 x, vec2 v) const;
	vec2 dfdv(vec2 x, vec2 v) const;
	mat2 df(vec2 x) const;
	mat2 df(float x, float y) const;
	vec2 operator()(vec2 x) const;
	vec2 operator()(float x, float y) const;
	virtual PlaneEndomorphism operator&(const PlaneEndomorphism &g) const;

	static PlaneEndomorphism linear(mat2 A);
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
    VectorFieldR3(Foo33 X, Foo3Foo33 dX, float eps=.01);
    explicit VectorFieldR3(Foo33 X, float eps=.01);
    VectorFieldR3(RealFunctionR3 Fx , RealFunctionR3 Fy, RealFunctionR3 Fz, float epsilon=0.01);
	explicit VectorFieldR3(Foo33 field);
	explicit VectorFieldR3(VectorFieldR2 f);
	VectorFieldR3 operator+(const VectorFieldR3 &Y) const;
	VectorFieldR3 operator*(float a) const;
    VectorFieldR3 operator-() const;
	VectorFieldR3 operator-(const VectorFieldR3 &Y) const;
	VectorFieldR3 operator*(const RealFunctionR3 &f) const;
    RealFunctionR3 F_x() const;
	RealFunctionR3 F_y() const;
	RealFunctionR3 F_z() const;
	std::array<RealFunctionR3, 3> components() const;
	R3 operator()(R3 v) const;

	friend VectorFieldR3 operator*(const mat3 &A, const VectorFieldR3 &X);

	static VectorFieldR3 constant(vec3 v);
	static VectorFieldR3 linear(mat3 A);
	static VectorFieldR3 radial(vec3 scale);

    RealFunctionR3 divergence() const;
    VectorFieldR3 curl() const;

    vec3 moveAlong(vec3 v, float dt=1) const;
};

class VectorFieldR3Dynamic {
	BIHOM(float, vec3, vec3) _X;
	float eps = 0.01;
	public:
	explicit VectorFieldR3Dynamic(BIHOM(float, vec3, vec3) X, float epsilon=0.01);
	vec3 operator()(float t, vec3 x) const;
	VectorFieldR3Dynamic operator+(const VectorFieldR3Dynamic &Y) const;
	VectorFieldR3Dynamic operator*(float a) const;
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
	SpaceAutomorphism inv() const;
	SpaceAutomorphism compose(SpaceAutomorphism g) const;
    SpaceAutomorphism operator&(const SpaceAutomorphism &g) const;
	SpaceAutomorphism applyWithBasis(mat3 A) const;
	SpaceAutomorphism applyWithBasis(vec3 v1, vec3 v2, vec3 v3) const;
	SpaceAutomorphism applyWithShift(vec3 v) const;

	static SpaceAutomorphism linear(mat3 A);
    static SpaceAutomorphism translation(vec3 v);
    static SpaceAutomorphism scaling(float x, float y, float z);
    static SpaceAutomorphism scaling(vec3 factors);

	static SpaceAutomorphism scaling(float x);
	static SpaceAutomorphism scaling(float x, float y, float z, vec3 center);
	static SpaceAutomorphism scaling(float x, vec3 center);
	static SpaceAutomorphism scaling(vec3 factors, vec3 center);

	static SpaceAutomorphism affine(mat3 A, vec3 v);
    static SpaceAutomorphism rotation(float angle);
    static SpaceAutomorphism rotation(vec3 axis, float angle);
	static SpaceAutomorphism rotation(vec3 axis, float angle, vec3 center);
	static SpaceAutomorphism deltaRotation(vec3 v1, vec3 v2);
	static SpaceAutomorphism deltaRotation(vec3 v1, vec3 v2, vec3 center);

};

class RealFunctionR1 {
	Fooo _f;
	Fooo _df;
	Fooo _ddf;
	float eps = 0.01;
public:
	RealFunctionR1(Fooo f, Fooo df, Fooo ddf, float epsilon=0.01);
	RealFunctionR1(Fooo f, Fooo df, float epsilon=0.01);
	explicit RealFunctionR1(Fooo f, float epsilon=0.01);

	RealFunctionR1(const RealFunctionR1 &other);
	RealFunctionR1(RealFunctionR1 &&other) noexcept;
	RealFunctionR1 & operator=(const RealFunctionR1 &other);
	RealFunctionR1 & operator=(RealFunctionR1 &&other) noexcept;

	float operator()(float x) const;
	float df(float x) const;
	float ddf(float x) const;
	RealFunctionR1 df() const;

	RealFunctionR1 operator+(const RealFunctionR1 &g) const;
	RealFunctionR1 operator*(float a) const;
	RealFunctionR1 operator+(float a) const;
	RealFunctionR1 operator-(float a) const;
	RealFunctionR1 operator*(const RealFunctionR1 &g_) const;
	RealFunctionR1 operator-() const;
	RealFunctionR1 operator-(const RealFunctionR1 &g) const;
	RealFunctionR1 operator/(const RealFunctionR1 &g_) const;
	RealFunctionR1 operator/(float a) const;
	friend RealFunctionR1 operator*(float a, const RealFunctionR1 &f);
	friend RealFunctionR1 operator+(float a, const RealFunctionR1 &f);
	friend RealFunctionR1 operator-(float a, const RealFunctionR1 &f);
	friend RealFunctionR1 operator/(float a, const RealFunctionR1 &f);

	RealFunctionR1 operator & (const RealFunctionR1 &g_) const;
	RealFunctionR3 operator & (const RealFunctionR3 &g_) const;

	static RealFunctionR1 constant(float a);
	static RealFunctionR1 x();
	static RealFunctionR1 one();
	static RealFunctionR1 zero();
	static RealFunctionR1 linear(float a, float b);
	static RealFunctionR1 quadratic(float a, float b, float c);
	static RealFunctionR1 monomial(int n);
	static RealFunctionR1 polynomial(std::vector<float> coeffs);
	static RealFunctionR1 sin();
	static RealFunctionR1 cos();
	static RealFunctionR1 exp();
	static RealFunctionR1 log();
	static RealFunctionR1 sqrt();
	static RealFunctionR1 pow(float a);
	static RealFunctionR1 pow(int a);
};

class RealLineAutomorphism : public RealFunctionR1 {
	RealFunctionR1 inverse;
public:
	RealLineAutomorphism(RealFunctionR1 f, RealFunctionR1 inv);
	RealLineAutomorphism(Fooo f, Fooo f_inv, float epsilon=0.01);
	RealLineAutomorphism(Fooo f, Fooo df, Fooo f_inv, float epsilon=0.01);
	RealLineAutomorphism(Fooo f, Fooo df, Fooo f_inv, Fooo df_inv, float epsilon=0.01);
	RealLineAutomorphism(Fooo f, Fooo df, Fooo ddf, Fooo f_inv, Fooo df_inv, Fooo ddf_inv, float epsilon=0.01);

	float inv(float x) const;
	RealLineAutomorphism inv() const;
	RealLineAutomorphism operator~() const;

	RealLineAutomorphism operator& (const RealLineAutomorphism &g) const;
	RealLineAutomorphism operator*(float c) const;
	RealLineAutomorphism operator+(float c) const;
	RealLineAutomorphism operator/(float c) const;
	RealLineAutomorphism operator-(float c) const;
	RealLineAutomorphism operator-() const;
	friend RealLineAutomorphism operator*(float c, const RealLineAutomorphism &f);
	friend RealLineAutomorphism operator+(float c, const RealLineAutomorphism &f);
	friend RealLineAutomorphism operator-(float c, const RealLineAutomorphism &f);

	static RealLineAutomorphism Id();
	static RealLineAutomorphism x();
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
	vec2 operator()(float t, float u) const;
	vec2 df(vec2 x, vec2 v) const;
	mat2 df(vec2 x) const;
	vec2 df(float t, float u, vec2 v) const;
	mat2 df(float t, float u) const;
};

class PlaneAutomorphism : public PlaneSmoothEndomorphism {
protected:
	std::shared_ptr<Foo22> _f_inv;
public:
	PlaneAutomorphism(std::shared_ptr<Foo22> f, std::shared_ptr<std::function<mat2(vec2)>> _df, std::shared_ptr<Foo22> f_inv);
	vec2 f_inv(vec2 x) const;
	vec2 inv(vec2 x) const;
	vec2 inv(float t, float u) const;
	PlaneAutomorphism operator~() const;
	PlaneAutomorphism inv() const;
};

template<typename A, typename B, typename C>
    HOM(A, C) compose(const HOM(B, C) &f, const HOM(A, B) &g) {
		return [f, g](A x) { return f(g(x)); };
	}
