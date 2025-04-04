#pragma once

#include <cmath>
#include <functional>
#include <iostream>
#include <memory>
#include <optional>
#include <string>
#include <utility>
#include <variant>


#include "filesUtils.hpp"
#include "monads.hpp"
#include "mat.hpp"
#include "modules.hpp"


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

class RealFunctionR3;
typedef RealFunctionR3 SteadyScalarField;


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
Regularity operator+(Regularity a, Regularity b);
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


class VectorField;
class RealFunction;

class RealFunctionR2 {
	Foo21 _f;
	Foo22 _df;
    float eps = 0.01;
    Regularity regularity;
public:
	RealFunctionR2();
    RealFunctionR2(const RealFunctionR2 &other);
    RealFunctionR2(RealFunctionR2 &&other) noexcept;
    RealFunctionR2 & operator=(const RealFunctionR2 &other);
	RealFunctionR2 & operator=(RealFunctionR2 &&other) noexcept;
	RealFunctionR2(Foo21 f, Foo22 df, float eps=.01, Regularity regularity = Regularity::SMOOTH);
	RealFunctionR2(Foo21 f, Foo21 dfdx, Foo21 dfdy, float eps=.01, Regularity regularity = Regularity::SMOOTH);
	explicit RealFunctionR2(Foo21 f, float epsilon=0.01, Regularity regularity = Regularity::SMOOTH);
	explicit RealFunctionR2(Foo111 f, float epsilon=0.01, Regularity regularity = Regularity::SMOOTH) : RealFunctionR2([f](vec2 x) { return f(x.x, x.y); }, epsilon, regularity) {}
	explicit RealFunctionR2(float constant, float epsilon=0.01, Regularity regularity = Regularity::SMOOTH) : RealFunctionR2([constant](vec2 x) { return constant; }, epsilon, regularity) {}

	float operator()(vec2 v) const;
	vec2 df(vec2 v) const;
	float operator()(float x, float y) const;
	vec2 df(float x, float y) const;


	RealFunctionR2 operator-() const;
	RealFunctionR2 operator*(float a) const;
	RealFunctionR2 operator+(float a) const;
	RealFunctionR2 operator-(float a) const;
	RealFunctionR2 operator/(float a) const;
	RealFunctionR2 operator+(const RealFunctionR2 &g) const;
	RealFunctionR2 operator-(const RealFunctionR2 &g) const;
	RealFunctionR2 operator*(const RealFunctionR2 &g) const;
	RealFunctionR2 operator/(const RealFunctionR2 &g) const;
	friend RealFunctionR2 operator*(float a, const RealFunctionR2 &f);
	friend RealFunctionR2 operator+(float a, const RealFunctionR2 &f);
	friend RealFunctionR2 operator-(float a, const RealFunctionR2 &f);
	friend RealFunctionR2 operator/(float a, const RealFunctionR2 &f);
	RealFunctionR2 operator~() const;
	RealFunctionR2 pow(float a) const;
	RealFunctionR2 sqrt() const;
	float getEps() const;

	float dx(vec2 x) const;
	float dx(float x, float y) const { return dx(vec2(x, y)); }
	RealFunctionR2 dx() const { return RealFunctionR2([F=*this](vec2 x) { return F.dx(x); }, eps, regularity-1); }
	float dy(vec2 x) const;
	float dy(float x, float y) const { return dy(vec2(x, y)); }
	RealFunctionR2 dy() const { return RealFunctionR2([this](vec2 x) { return dy(x); }, eps, regularity-1); }
	RealFunctionR2 dxi_yj(int i, int j) const;
	float dxi_yj(int i, int j, vec2 x) const { return dxi_yj(i, j)(x); }
	float dxi_yj(int i, int j, float x, float y) const { return dxi_yj(i, j, vec2(x, y)); }
	float dxx(vec2 x) const { return dxi_yj(2, 0, x); }
	float dxx(float x, float y) const { return dxx(vec2(x, y)); }
	float dyy(vec2 x) const { return dxi_yj(0, 2, x); }
	float dyy(float x, float y) const { return dyy(vec2(x, y)); }
	float dxy(vec2 x) const { return dxi_yj(1, 1, x); }
	float dxy(float x, float y) const { return dxy(vec2(x, y)); }
	float Laplacian(vec2 x) const { return dxx(x) + dyy(x); }

	float integrate_rect(vec2 a, vec2 b, int precision = 100) const;
	RealFunction partially_evaulate(int variable_ind, float value) const;
	RealFunction partially_integrate_along_x(float x0, float x1, int precision = 100) const;
	RealFunction partially_integrate_along_y(float x0, float x1, int precision = 100) const;
	RealFunctionR2 convolve_x(const RealFunction &g, float L, int precision = 100) const;
	RealFunctionR2 convolve_y(const RealFunction &g, float L, int precision = 100) const;

	static RealFunctionR2 linear(vec2 v);
	static RealFunctionR2 constant(float a);
	static RealFunctionR2 projection(int i);
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
	RealFunctionR3 operator-() const;

	RealFunctionR3 operator*(float a) const;
    RealFunctionR3 operator+(float a) const;
    RealFunctionR3 operator-(float a) const;
    RealFunctionR3 operator/(float a) const;
	RealFunctionR3 operator+(const RealFunctionR3 &g) const;
	RealFunctionR3 operator-(const RealFunctionR3 &g) const;
	RealFunctionR3 operator*(const RealFunctionR3 &g) const;
	RealFunctionR3 operator/(const RealFunctionR3 &g) const;
    friend RealFunctionR3 operator*(float a, const RealFunctionR3 &f);
	friend RealFunctionR3 operator+(float a, const RealFunctionR3 &f);
	friend RealFunctionR3 operator-(float a, const RealFunctionR3 &f);
	friend RealFunctionR3 operator/(float a, const RealFunctionR3 &f);
	RealFunctionR3 operator~() const;
	RealFunctionR3 pow(float a) const;
	RealFunctionR3 sqrt() const;
	float getEps() const;

	VectorField gradient() const;
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




class ScalarField {
	BIHOM(vec3, float, float) F;
	float eps;
public:
	ScalarField();
	explicit ScalarField(BIHOM(vec3, float, float) F, float eps=.01);
	explicit ScalarField(HOM(float, SteadyScalarField) pencil);
	explicit ScalarField(SteadyScalarField steady_field);

	float operator()(vec3 x, float t) const;
	SteadyScalarField operator()(float t) const;
	SteadyScalarField fix_time(float t) const;

	ScalarField operator+(const ScalarField &Y) const;
	ScalarField operator*(float a) const;
	ScalarField operator-() const;
	ScalarField operator-(const ScalarField &Y) const;
	ScalarField operator*(const SteadyScalarField &f) const;
	ScalarField operator/(const SteadyScalarField &f) const;
	ScalarField operator*(const ScalarField &f) const;
	ScalarField operator/(const ScalarField &f) const;

	ScalarField time_derivative() const;
	ScalarField spatial_partial(int i) const;
};




	// 	Flow spatial_partial_derivative(int i) const {
// 		return Flow([F_=F, i, e=eps](vec3 x, float t){
// 			return partialDerivativeOperator([F=F_, i=i, t](vec3 y) {
// 				return F(y, t);
// 			}, i, e);
// }
// 	}

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

class CompactlySupportedRealFunction;


/**
 * @class RealFunction
 * @brief Encoding real line continuous endomorphisms
 *
 *
 *  @note precomposition syntax
 *  - both () and & operators act by precomposition on RealFunction
 *		@example
 *  @code
 *  (f & g)(x) = f(g)(x) = f(g(x))
 *  @endcode
 *  - on real numbers () act as evaluation (and thus differs from &)
 *  - on real numbers & act as precomposition with linear reparametrization
 *	@example
 *  @code
 *  (f & a)(x) = f(X_R *a)(x) = f(a *x)
 *  @endcode
 *
 *  @note Predefined elementary functions
 *	@example
 *  @code
 *  ONE_R
 *  X_R
 *  EXP_R
 *  LOG_R
 *  SIN_R
 *  COS_R
 *  SQRT_R
 *  @endcode
 *
 *   @example
 *  @code
 *  SQRT_R == pow(X_R, 0.5)
 *  SQRT_R & F == F.sqrt()
 *  F*F* F == F.pow3()
 *  SQRT & (F*F) == F
 *  F & 2 == F(X_R*2)
 *  F & 2 == F & (X_R & 2)
 *  X_R & 2 == X_R*2
 *  SQRT_R & 4 == SQRT_R*2
 *  @endcode
 */
class RealFunction {
protected:
	Fooo _f;
	Fooo _df;
	Fooo _ddf;
	float eps = 0.01;
	bool is_zero;
public:
	Regularity regularity;

	RealFunction(Fooo f, Fooo df, Fooo ddf, float epsilon=0.01, Regularity regularity = Regularity::SMOOTH);
	RealFunction(Fooo f, Fooo df, float epsilon=0.01, Regularity regularity = Regularity::SMOOTH);
	explicit RealFunction(Fooo f, float epsilon=0.01, Regularity regularity = Regularity::SMOOTH);
	explicit RealFunction(float constant, float epsilon=0.01, Regularity regularity = Regularity::SMOOTH);

	RealFunction(const RealFunction &other);
	RealFunction(RealFunction &&other) noexcept;
	RealFunction & operator=(const RealFunction &other);
	RealFunction & operator=(RealFunction &&other) noexcept;

	float operator()(float x) const;
	float df(float x) const;
	float ddf(float x) const;
	RealFunction df() const;
	float getEps() const { return eps; }


	bool isZero(float a=-10, float b=10, int prec=23);
	float sup(float a, float b, int prec) const;
	float inf(float a, float b, int prec) const;
	float argmax(float a, float b, int prec) const;
	float argmin(float a, float b, int prec) const;

	RealFunction operator+(const RealFunction &g) const;
	RealFunction operator*(float a) const;
	RealFunction operator+(float a) const;
	RealFunction operator-(float a) const;
	RealFunction operator*(const RealFunction &g_) const;
	RealFunction operator-() const;
	RealFunction operator-(const RealFunction &g) const;
	RealFunction operator/(const RealFunction &g_) const;
	RealFunction operator/(float a) const;
	RealFunction pow(float a) const;

	friend RealFunction pow(const RealFunction &f, float a);
	friend RealFunction max(const RealFunction &f, const RealFunction &g);
	friend RealFunction max(const RealFunction &f, float a);
	friend RealFunction min(const RealFunction &f, const RealFunction &g);
	friend RealFunction min(const RealFunction &f, float a);
	friend RealFunction abs(const RealFunction &f) {return max(f, -f);}

	RealFunction sqrt() const { return pow(.5f); }
	RealFunction pow2() const { return pow(2); }
	RealFunction pow3() const { return pow(3); }
	RealFunction pow4() const { return pow(4); }

	friend RealFunction operator*(float a, const RealFunction &f);
	friend RealFunction operator+(float a, const RealFunction &f);
	friend RealFunction operator-(float a, const RealFunction &f);
	friend RealFunction operator/(float a, const RealFunction &f);
	// RealFunction operator~() const { return 1.f/(*this); }

	RealFunction operator& (const RealFunction &g_) const;
	RealFunctionR3 operator& (const RealFunctionR3 &g_) const;
	RealFunction operator& (float a) const;
	RealFunction operator()(const RealFunction &g_) const;

	float integral(float a, float b, int prec) const;
	RealFunction antiderivative(float a, int prec) const;
	float L2_norm(vec2 I, int prec) const;
	float L2_product(const RealFunction &g, vec2 I, int prec) const;
	RealFunction convolve(CompactlySupportedRealFunction kernel, int prec=1000) const;
	RealFunction convolve(const RealFunction &kernel, float L, int prec=100) const;
	float repeated_integral(float a, float b, int n, int prec) const;
	RealFunction repeated_antiderivative(float a, int prec) const; // calculated using Cauchy formula for repeated integrals

	friend RealFunctionR2 separated_product(const RealFunction &f_x, const RealFunction &f_t);
	friend RealFunctionR2 separated_sum(const RealFunction &f_x, const RealFunction &f_t);


	static RealFunction constant(float a);
	static RealFunction x();
	static RealFunction one();
	static RealFunction zero();
	static RealFunction linear(float a, float b);
	static RealFunction quadratic(float a, float b, float c);
	static RealFunction monomial(float n);
	static RealFunction polynomial(std::vector<float> coeffs);
	static RealFunction sin();
	static RealFunction cos();
	static RealFunction exp();
	static RealFunction log();
	static RealFunction SQRT();
};

const auto SIN_R = RealFunction::sin();
const auto COS_R = RealFunction::cos();
const auto EXP_R = RealFunction::exp();
const auto LOG_R = RealFunction::log();
const auto SQRT_R = RealFunction::SQRT();
const auto ONE_R = RealFunction::one();
const auto X_R = RealFunction::linear(1, 0);
const auto RELU = max(X_R, 0);


class CompactlySupportedRealFunction: public RealFunction {
	vec2 support;
	vec2 support_sum(vec2 other) const;
	vec2 support_intersect(vec2 other) const;
	int default_interval_arithmetic_prec=100;

public:
	CompactlySupportedRealFunction(Fooo f, Fooo df, Fooo ddf, vec2 support, float epsilon=0.01, Regularity regularity = Regularity::SMOOTH);
	CompactlySupportedRealFunction(Fooo f, Fooo df, vec2 support, float epsilon=0.01, Regularity regularity = Regularity::SMOOTH);
	CompactlySupportedRealFunction(Fooo f, vec2 support, float epsilon=0.01, Regularity regularity = Regularity::SMOOTH);
	CompactlySupportedRealFunction(const RealFunction &f, vec2 support);
	CompactlySupportedRealFunction(const RealFunction &f, vec2 support, float epsilon);
	CompactlySupportedRealFunction(const RealFunction &f, const RealFunction &df, vec2 support, float epsilon);
	CompactlySupportedRealFunction(const RealFunction &f, const RealFunction &df, const RealFunction &ddf, vec2 support, float epsilon);

	CompactlySupportedRealFunction df() const;
	float improper_integral(RP1 a, RP1 b, int prec) const;
	float full_domain_integral(int prec) const;
	CompactlySupportedRealFunction antiderivative(float a, int prec) const;

	CompactlySupportedRealFunction operator+(const CompactlySupportedRealFunction &g) const;
	CompactlySupportedRealFunction operator*(float a) const;
	CompactlySupportedRealFunction operator/(float a) const;
	CompactlySupportedRealFunction operator*(const CompactlySupportedRealFunction &g_) const;
	CompactlySupportedRealFunction operator*(const RealFunction &g_) const;

	CompactlySupportedRealFunction operator-() const;
	CompactlySupportedRealFunction operator-(const CompactlySupportedRealFunction &g) const;
	friend CompactlySupportedRealFunction operator*(float a, const CompactlySupportedRealFunction &f);
	friend CompactlySupportedRealFunction operator*(const RealFunction &g, const CompactlySupportedRealFunction &f);



	RealFunction convolve(const RealFunction &f, int prec=1000) const {
		auto K = RealFunctionR2 ([f, kernel=_f](vec2 v) {
			return kernel(v.x-v.y) * f(v.y);
		}, eps, regularity + f.regularity);
		return K.partially_integrate_along_y(support[0], support[1], prec);
	}
};

class RealLineAutomorphism : public RealFunction {
	RealFunction inverse;
public:
	RealLineAutomorphism(RealFunction f, RealFunction inv);
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
	HOM(A, C) compose(const HOM(B, C) &f, const HOM(A, B) &g);


class HolomorphicFunction {
	HOM(Complex, Complex) _f, _df;
	HolomorphicFunction(std::function<Complex(Complex)> f, std::function<Complex(Complex)> df, float epsilon);
	HolomorphicFunction(std::function<Complex(Complex)> f, float epsilon);
	HolomorphicFunction(const PlaneSmoothEndomorphism &f, float epsilon);
	HolomorphicFunction(const RealFunctionR2 &re, const RealFunctionR2 &im, float epsilon);
	HolomorphicFunction df() const;
	HolomorphicFunction d2f() const;
	HolomorphicFunction dnf(int n) const;
	Complex df(Complex z) const;
	Complex d2f(Complex z) const;
	Complex dnf(Complex z, int n) const;
	HolomorphicFunction operator-() const;
	HolomorphicFunction pow(float a) const;
	HolomorphicFunction sqrt() const;
	HolomorphicFunction sq() const;
	HolomorphicFunction pow2() const;
	HolomorphicFunction pow3() const;
	HolomorphicFunction operator*(Complex a) const;
	HolomorphicFunction operator+(Complex a) const;
	Complex operator()(Complex z) const;
};

class ComplexValuedFunction {
	HOM(float, Complex) f;
	float eps;
public:
	explicit ComplexValuedFunction(Complex c) : f([c](float x) { return c; }), eps(0.01) {}
	explicit ComplexValuedFunction(HOM(float, Complex) f, float eps=.01);
	explicit ComplexValuedFunction(RealFunction re, float eps=.01);
	ComplexValuedFunction(RealFunction f_re, RealFunction f_im, float eps=.01);
	ComplexValuedFunction(RealFunction f_power, RealFunction f_phase, bool dump);

	explicit ComplexValuedFunction(HOM(float, float) re, float eps=.01);
	ComplexValuedFunction(HOM(float, float) f_re, HOM(float, float) f_im, float eps=.01);

	ComplexValuedFunction(const ComplexValuedFunction &other);
	ComplexValuedFunction(ComplexValuedFunction &&other) noexcept;
	ComplexValuedFunction & operator=(const ComplexValuedFunction &other);
	ComplexValuedFunction & operator=(ComplexValuedFunction &&other) noexcept;



	Complex operator()(float x) const { return f(x); }
	ComplexValuedFunction df() const { return ComplexValuedFunction([f=f, eps=eps](float x) { return (f(x+eps) - f(x-eps))/(2.f*eps); }, eps); }
	ComplexValuedFunction conj() const { return ComplexValuedFunction(re(), -im(), eps); }

	ComplexValuedFunction operator+(const ComplexValuedFunction &g) const { return ComplexValuedFunction([f=f, g=g](float x) { return f(x) + g(x); }, eps); }
	ComplexValuedFunction operator+(const RealFunction &g) const { return ComplexValuedFunction([f=f, g=g](float x) { return f(x) + g(x); }, eps); }
	ComplexValuedFunction operator+(const Complex &g) const { return ComplexValuedFunction([f=f, g=g](float x) { return f(x) + g; }, eps); }
	friend ComplexValuedFunction operator+(const RealFunction &g, const ComplexValuedFunction &f) { return f + g; }
    friend ComplexValuedFunction operator+(const Complex &g, const ComplexValuedFunction &f) { return f + g; }

	ComplexValuedFunction operator-() const { return ComplexValuedFunction([f=f](float x) { return -f(x); }, eps); }
	ComplexValuedFunction operator-(const ComplexValuedFunction &g) const { return (*this) + (-g); }
	ComplexValuedFunction operator-(const RealFunction &g) const { return (*this) + (-g); }
	ComplexValuedFunction operator-(const Complex &g) const { return (*this) + (-g); }
	friend ComplexValuedFunction operator-(const RealFunction &g, const ComplexValuedFunction &f) { return -f + g; }
	friend ComplexValuedFunction operator-(const Complex &g, const ComplexValuedFunction &f) { return -f + g; }

	ComplexValuedFunction operator*(const ComplexValuedFunction &g) const { return ComplexValuedFunction([f=f, g=g](float x) { return f(x) * g(x); }, eps); }
	ComplexValuedFunction operator*(const RealFunction &g) const { return ComplexValuedFunction([f=f, g=g](float x) { return f(x) * g(x); }, eps); }
	ComplexValuedFunction operator*(const Complex &g) const { return ComplexValuedFunction([f=f, g=g](float x) { return f(x) * g; }, eps); }
	friend ComplexValuedFunction operator*(const RealFunction &g, const ComplexValuedFunction &f) { return f * g; }
	friend ComplexValuedFunction operator*(const Complex &g, const ComplexValuedFunction &f) { return f * g; }

	ComplexValuedFunction operator/(const ComplexValuedFunction &g) const { return ComplexValuedFunction([f=f, g=g](float x) { return f(x) / g(x); }, eps); }
	ComplexValuedFunction operator/(const RealFunction &g) const { return ComplexValuedFunction([f=f, g=g](float x) { return f(x) / g(x); }, eps); }
	ComplexValuedFunction operator/(const Complex &g) const { return ComplexValuedFunction([f=f, g=g](float x) { return f(x) / g; }, eps); }
	friend ComplexValuedFunction operator/(const RealFunction &g, const ComplexValuedFunction &f) { return ComplexValuedFunction([f=f, g=g](float x) { return Complex(g(x)) / f(x); }, f.eps); }
	friend ComplexValuedFunction operator/(const Complex &g, const ComplexValuedFunction &f) { return ComplexValuedFunction([f=f, g=g](float x) { return Complex(g) / f(x); }, f.eps); }

	ComplexValuedFunction operator & (float t) const { return ComplexValuedFunction([f=f, t](float x) { return f(x*t); }, eps); }
	ComplexValuedFunction operator & (const RealFunction &g) const { return ComplexValuedFunction([f=f, g=g](float x) { return f(g(x)); }, eps); }
	ComplexValuedFunction operator() (const RealFunction &g) const { return (*this) & g; }

	RealFunction re() const { return RealFunction([f=f](float x) { return f(x).real(); }, eps); }
	RealFunction im() const { return RealFunction([f=f](float x) { return f(x).imag(); }, eps); }
	RealFunction abs() const { return RealFunction([f=f](float x) { return norm(f(x)); }, eps); }
	RealFunction arg() const { return RealFunction([f=f](float x) { return f(x).arg(); }, eps); }
	RealFunction norm2() const { return RealFunction([f=f](float x) { return ::norm2(f(x)); }, eps); }

	Complex integral(float a, float b, int prec) const { return Complex(re().integral(a, b, prec), im().integral(a, b, prec)); }
};


class ComplexFunctionR2 {
	BIHOM(float, float,  Complex) f;
	float eps;
public:
	explicit ComplexFunctionR2(Complex c) : f([c](float x, float y) { return c; }), eps(0.01) {}
	explicit ComplexFunctionR2(BIHOM(float, float,  Complex) f, float eps=.01) : f(std::move(f)), eps(eps) {}
	explicit ComplexFunctionR2(RealFunctionR2 re, float eps=.01) : ComplexFunctionR2([re](float x, float y) { return Complex(re(x, y), 0); }, eps) {}
	ComplexFunctionR2(RealFunctionR2 f_re, RealFunctionR2 f_im, float eps=.01) : ComplexFunctionR2([f_re, f_im](float x, float y) { return Complex(f_re(x, y), f_im(x, y)); }, eps) {}
	explicit ComplexFunctionR2(HOM(vec2,  Complex) f, float eps=.01) : ComplexFunctionR2([f](float x, float y) { return f(vec2(x, y)); }, eps) {}


	ComplexFunctionR2(const ComplexFunctionR2 &other);
	ComplexFunctionR2(ComplexFunctionR2 &&other) noexcept;
	ComplexFunctionR2 & operator=(const ComplexFunctionR2 &other);
	ComplexFunctionR2 & operator=(ComplexFunctionR2 &&other) noexcept;


	Complex operator()(float x, float y) const { return f(x, y); }
	Complex operator()(vec2 x) const { return f(x.x, x.y); }
	ComplexFunctionR2 dfdx() const { return ComplexFunctionR2([f=f, eps=eps](float x, float y) { return (f(x+eps, y) - f(x-eps, y))/(2.f*eps); }, eps); }
	ComplexFunctionR2 conj() const { return ComplexFunctionR2(re(), -im(), eps); }

	ComplexFunctionR2 operator+(const ComplexFunctionR2 &g) const { return ComplexFunctionR2([f=f, g=g](float x, float y) { return f(x, y) + g(x, y); }, eps); }
	ComplexFunctionR2 operator+(const RealFunctionR2 &g) const { return ComplexFunctionR2([f=f, g=g](float x, float y) { return f(x, y) + g(x, y); }, eps); }
	ComplexFunctionR2 operator+(const Complex &g) const { return ComplexFunctionR2([f=f, g=g](float x, float y) { return f(x, y) + g; }, eps); }
	friend ComplexFunctionR2 operator+(const RealFunctionR2 &g, const ComplexFunctionR2 &f) { return f + g; }
    friend ComplexFunctionR2 operator+(const Complex &g, const ComplexFunctionR2 &f) { return f + g; }

	ComplexFunctionR2 operator-() const { return ComplexFunctionR2([f=f](float x, float y) { return -f(x, y); }, eps); }
	ComplexFunctionR2 operator-(const ComplexFunctionR2 &g) const { return (*this) + (-g); }
	ComplexFunctionR2 operator-(const RealFunctionR2 &g) const { return (*this) + (-g); }
	ComplexFunctionR2 operator-(const Complex &g) const { return (*this) + (-g); }
	friend ComplexFunctionR2 operator-(const RealFunctionR2 &g, const ComplexFunctionR2 &f) { return -f + g; }
	friend ComplexFunctionR2 operator-(const Complex &g, const ComplexFunctionR2 &f) { return -f + g; }

	ComplexFunctionR2 operator*(const ComplexFunctionR2 &g) const { return ComplexFunctionR2([f=f, g=g](float x, float y) { return f(x, y) * g(x, y); }, eps); }
	ComplexFunctionR2 operator*(const RealFunctionR2 &g) const { return ComplexFunctionR2([f=f, g=g](float x, float y) { return f(x, y) * g(x, y); }, eps); }
	ComplexFunctionR2 operator*(const Complex &g) const { return ComplexFunctionR2([f=f, g=g](float x, float y) { return f(x, y) * g; }, eps); }
	friend ComplexFunctionR2 operator*(const RealFunctionR2 &g, const ComplexFunctionR2 &f) { return f * g; }
	friend ComplexFunctionR2 operator*(const Complex &g, const ComplexFunctionR2 &f) { return f * g; }

	ComplexFunctionR2 operator/(const ComplexFunctionR2 &g) const { return ComplexFunctionR2([f=f, g=g](float x, float y) { return f(x, y) / g(x, y); }, eps); }
	ComplexFunctionR2 operator/(const RealFunctionR2 &g) const { return ComplexFunctionR2([f=f, g=g](float x, float y) { return f(x, y) / g(x, y); }, eps); }
	ComplexFunctionR2 operator/(const Complex &g) const { return ComplexFunctionR2([f=f, g=g](float x, float y) { return f(x, y) / g; }, eps); }
	friend ComplexFunctionR2 operator/(const RealFunctionR2 &g, const ComplexFunctionR2 &f) { return ComplexFunctionR2([f=f, g=g](float x, float y) { return Complex(g(x, y)) / f(x, y); }, f.eps); }
	friend ComplexFunctionR2 operator/(const Complex &g, const ComplexFunctionR2 &f) { return ComplexFunctionR2([f=f, g=g](float x, float y) { return Complex(g) / f(x, y); }, f.eps); }

	RealFunctionR2 re() const { return RealFunctionR2([f=f](float x, float y) { return f(x, y).real(); }, eps); }
	RealFunctionR2 im() const { return RealFunctionR2([f=f](float x, float y) { return f(x, y).imag(); }, eps); }
	RealFunctionR2 abs() const { return RealFunctionR2([f=f](float x, float y) { return norm(f(x, y)); }, eps); }
	RealFunctionR2 arg() const { return RealFunctionR2([f=f](float x, float y) { return f(x, y).arg(); }, eps); }
	RealFunctionR2 norm2() const { return RealFunctionR2([f=f](float x, float y) { return ::norm2(f(x, y)); }, eps); }
};




const auto EXP_it = ComplexValuedFunction([](float t) { return exp(t*1.0i); }, 0.01);
const auto LOG_RC = ComplexValuedFunction([](float t) { return log(Complex(t)); }, 0.01);
const auto SQRT_RC = ComplexValuedFunction([](float t) { return sqrt(Complex(t)); }, 0.01);
const auto ONE_RC = ComplexValuedFunction(1.f);
const auto X_RC = ComplexValuedFunction([](float t) { return Complex(t); }, 0.01);
const auto SIN_RC = ComplexValuedFunction([](float t) { return sin(Complex(t)); }, 0.01);
const auto COS_RC = ComplexValuedFunction([](float t) { return cos(Complex(t)); }, 0.01);
const auto EXP_RC = ComplexValuedFunction([](float t) { return exp(Complex(t)); }, 0.01);


class DiscreteRealFunction;

class DiscreteComplexFunction {
	Vector<Complex> fn;
	vec2 domain;
public:
	DiscreteComplexFunction(const Vector<Complex> &fn, vec2 domain) : fn(fn), domain(domain) {}
	DiscreteComplexFunction(const HOM(float, Complex) &fn, vec2 domain, int sampling);

	int samples() const { return fn.size(); }
	float sampling_step() const { return (domain[1] - domain[0]) / (samples()-1); }
	Complex operator[](int i) const { return fn[i]; }
	Vector<Complex> getVector() const { return fn; }
	vec2 getDomain() const { return domain; }


	DiscreteComplexFunction(const DiscreteComplexFunction &other);
	DiscreteComplexFunction(DiscreteComplexFunction &&other) noexcept;
	DiscreteComplexFunction & operator=(const DiscreteComplexFunction &other);
	DiscreteComplexFunction & operator=(DiscreteComplexFunction &&other) noexcept;
	explicit DiscreteComplexFunction(const DiscreteRealFunction &other);

	DiscreteComplexFunction operator+(const DiscreteComplexFunction &g) const;;
	DiscreteComplexFunction operator+(Complex a) const;;
	DiscreteComplexFunction operator+(float a) const;;
	friend DiscreteComplexFunction operator+(float a, const DiscreteComplexFunction &f) { return f + a; }
	friend DiscreteComplexFunction operator+(Complex a, const DiscreteComplexFunction &f) { return f + a; }

	DiscreteComplexFunction operator-() const { return DiscreteComplexFunction(-fn, domain); }
	DiscreteComplexFunction operator-(const DiscreteComplexFunction &g) const { return (*this) + (-g); }
	DiscreteComplexFunction operator-(float a) const { return (*this) + (-a); }
	DiscreteComplexFunction operator-(Complex a) const { return (*this) + (-a); }
	friend DiscreteComplexFunction operator-(float a, const DiscreteComplexFunction &f) { return f + -a; }
	friend DiscreteComplexFunction operator-(Complex a, const DiscreteComplexFunction &f) { return f + -a; }

	DiscreteComplexFunction operator*(const DiscreteComplexFunction &g) const;
	DiscreteComplexFunction operator*(Complex a) const;
	DiscreteComplexFunction operator*(float a) const {return this->operator*(Complex(a));}
	friend DiscreteComplexFunction operator*(float a, const DiscreteComplexFunction &f) { return f * a; }
	friend DiscreteComplexFunction operator*(Complex a, const DiscreteComplexFunction &f) { return f * a; }


	DiscreteComplexFunction operator/(const DiscreteComplexFunction &g) const;
	DiscreteComplexFunction operator/(float a) const { return this->operator/(Complex(1.0/a)); }
	DiscreteComplexFunction operator/(Complex a) const { return this->operator*(1.0/a); }

	Complex operator()(float x) const;
	DiscreteRealFunction re() const;
	DiscreteRealFunction im() const;
	DiscreteRealFunction abs() const;
	DiscreteRealFunction arg() const;
	DiscreteComplexFunction conj() const;
	vector<float> args() const { return linspace(domain[0], domain[1], samples()); }


	DiscreteComplexFunction fft() const;
	DiscreteComplexFunction ifft() const;

	Complex integral() const { return sum<Complex>(fn) * sampling_step(); }
	float L2_norm() const { return norm((*this * conj()).integral()); }

	float supp_len() const { return domain[1] - domain[0]; }

	DiscreteComplexFunction shift_domain_left() const;
	DiscreteComplexFunction shift_domain_right() const;
	bool starts_from_zero() const { return isClose(0, domain[0]); }
	bool symmetric() const { return isClose(domain[0], -domain[1]); }
};



class DiscreteRealFunction {
	Vector<float> fn;
	vec2 domain;
public:
	DiscreteRealFunction() = default;
	DiscreteRealFunction(const Vector<float> &fn, vec2 domain);
	DiscreteRealFunction(const HOM(float, float) &f, vec2 domain, int sampling);

	int samples() const { return fn.size(); }
	float sampling_step() const { return (domain[1] - domain[0]) / (samples()-1); }
	float operator[](int i) const { return fn[i]; }
	Vector<float> getVector() const { return fn; }
	vec2 getDomain() const { return domain; }

	DiscreteRealFunction(const DiscreteRealFunction &other);
	DiscreteRealFunction(DiscreteRealFunction &&other) noexcept;
	DiscreteRealFunction & operator=(const DiscreteRealFunction &other);
	DiscreteRealFunction & operator=(DiscreteRealFunction &&other) noexcept;

	vector<float> args() const { return linspace(domain[0], domain[1], samples()); }

	DiscreteRealFunction operator+(const DiscreteRealFunction &g) const;
	DiscreteRealFunction operator+(float a) const;
	friend DiscreteRealFunction operator+(float a, const DiscreteRealFunction &f) { return f + a; }
	DiscreteRealFunction operator-() const { return DiscreteRealFunction(-fn, domain); }
	DiscreteRealFunction operator-(const DiscreteRealFunction &g) const { return (*this) + (-g); }
	DiscreteRealFunction operator-(float a) const { return (*this) + (-a); }
	friend DiscreteRealFunction operator-(float a, const DiscreteRealFunction &f) { return f + -a; }
	DiscreteRealFunction operator*(const DiscreteRealFunction &g) const { return DiscreteRealFunction(fn.pointwise_product(g.fn), domain); }
	DiscreteRealFunction operator*(float a) const { return DiscreteRealFunction(fn * a, domain); }
	friend DiscreteRealFunction operator*(float a, const DiscreteRealFunction &f) { return f * a; }
	DiscreteRealFunction operator/(const DiscreteRealFunction &g) const;
	DiscreteRealFunction operator/(float a) const { return this->operator*(1.0/a); }

	float operator()(float x) const;

	DiscreteRealFunction two_sided_zero_padding(int target_size) const;

	DiscreteComplexFunction fft() const {return DiscreteComplexFunction(*this).fft();}

	DiscreteRealFunction convolve(const DiscreteRealFunction &kernel) const;
	DiscreteRealFunction derivative() const;
	float integral() const { return sum(fn) * sampling_step(); }
	float L2_norm() const { return sqrt((*this * *this).integral()); }
	float integral(int b) const { return sum(fn.slice_to(b)) * sampling_step(); }
	float integral(int a, int b) const { return sum(fn.slice(a, b)) * sampling_step(); }



	float supp_len() const { return domain[1] - domain[0]; }

	DiscreteRealFunction shift_domain_left() const;
	DiscreteRealFunction shift_domain_right() const;
	bool starts_from_zero() const { return isClose(0, domain[0]); }
	bool symmetric() const { return isClose(domain[0], -domain[1]); }

};

class DiscreteRealFunctionR2;


class DiscreteComplexFunctionR2 {
	Vector<DiscreteComplexFunction> fn;
	vec2 domain;
public:
	DiscreteComplexFunctionR2(const vector<DiscreteComplexFunction> &fn, vec2 domain) : fn(fn), domain(domain) {}
	DiscreteComplexFunctionR2(const vector<HOM(float, Complex)> &f, vec2 domain, int sampling);
	DiscreteComplexFunctionR2(const DiscreteComplexFunctionR2 &other) = default;
	DiscreteComplexFunctionR2(DiscreteComplexFunctionR2 &&other) noexcept = default;
	DiscreteComplexFunctionR2 & operator=(const DiscreteComplexFunctionR2 &other) = default;
	DiscreteComplexFunctionR2 & operator=(DiscreteComplexFunctionR2 &&other) noexcept = default;

	DiscreteComplexFunction operator[](int t) const { return fn[t]; }
	int samples_t () const { return fn.size(); }
	int samples_x () const { return fn[0].samples(); }
	float sampling_step_t() const { return (domain[1] - domain[0]) / (samples_t()-1); }
	float sampling_step_x() const { return fn[0].sampling_step(); }
	vector<float> args_t() const { return linspace(domain[0], domain[1], samples_t()); }
	vector<float> args_x() const { return fn[0].args(); }

	DiscreteComplexFunctionR2 operator+(const DiscreteComplexFunctionR2 &g) const;
	DiscreteComplexFunctionR2 operator+(Complex a) const;
	DiscreteComplexFunctionR2 operator+(float a) const;
	friend DiscreteComplexFunctionR2 operator+(float a, const DiscreteComplexFunctionR2 &f) { return f + a; }
	friend DiscreteComplexFunctionR2 operator+(Complex a, const DiscreteComplexFunctionR2 &f) { return f + a; }
	DiscreteComplexFunctionR2 operator-() const { return DiscreteComplexFunctionR2((-fn).vec(), domain); }
	DiscreteComplexFunctionR2 operator-(const DiscreteComplexFunctionR2 &g) const { return (*this) + (-g); }
	DiscreteComplexFunctionR2 operator-(float a) const { return (*this) + (-a); }
	DiscreteComplexFunctionR2 operator-(Complex a) const { return (*this) + (-a); }
	friend DiscreteComplexFunctionR2 operator-(float a, const DiscreteComplexFunctionR2 &f) { return f + -a; }
	friend DiscreteComplexFunctionR2 operator-(Complex a, const DiscreteComplexFunctionR2 &f) { return f + -a; }
	DiscreteComplexFunctionR2 operator*(const DiscreteComplexFunctionR2 &g) const;
	DiscreteComplexFunctionR2 operator*(Complex a) const;
	DiscreteComplexFunctionR2 operator*(float a) const { return this->operator*(Complex(a)); }
	friend DiscreteComplexFunctionR2 operator*(float a, const DiscreteComplexFunctionR2 &f) { return f * a; }
	friend DiscreteComplexFunctionR2 operator*(Complex a, const DiscreteComplexFunctionR2 &f) { return f * a; }
	DiscreteComplexFunctionR2 operator/(const DiscreteComplexFunctionR2 &g) const;
	DiscreteComplexFunctionR2 operator/(float a) const { return this->operator*(1.0/a); }
	DiscreteComplexFunctionR2 operator/(Complex a) const { return this->operator*(1.0/a); }
	Complex operator()(float t, float x) const { return fn[t](x); }
	Complex operator()(vec2 x) const { return fn[x.x](x.y); }




	DiscreteComplexFunctionR2 integrate_t(int i_t) const;
	DiscreteComplexFunctionR2 integrate_x(int i_x) const;
	DiscreteComplexFunctionR2 transpose() const;
	Complex double_integral() const;
	Complex double_integral(int i_t) const;
	Complex double_integral(int i_t, int i_x) const;

	DiscreteComplexFunctionR2 fft_t() const {
		vector<DiscreteComplexFunction> res = vector<DiscreteComplexFunction>();
		for (auto &f : fn) {
			res.push_back(f.fft());
		}
		return DiscreteComplexFunctionR2(res, domain);
	}

	DiscreteComplexFunctionR2 ifft_t() const {
		vector<DiscreteComplexFunction> res = vector<DiscreteComplexFunction>();
		for (auto &f : fn) {
			res.push_back(f.ifft());
		}
		return DiscreteComplexFunctionR2(res, domain);
	}
	DiscreteComplexFunctionR2 fft_x() const {
		return transpose().fft_t().transpose();
	}
	DiscreteComplexFunctionR2 ifft_x() const {
		return transpose().ifft_t().transpose();
	}
	DiscreteComplexFunctionR2 fft() const {
		float Lt=domain[1]-domain[0];
		float Lx=fn[0].getDomain()[1]-fn[0].getDomain()[0];
		vector<DiscreteComplexFunction> res = vector<DiscreteComplexFunction>();

		for (int it=0; it<samples_t(); ++it) {
			auto v = Vector<Complex>(samples_x(), [fn=fn, nx=samples_x(), Lx, Lt, nt=samples_t(), it](int ix) {
				Complex c = 0;
				for (int jx=0; jx<nx; ++jx)
					for (int kt=0; kt<nt; ++kt)
						c += exp(-1.0i * (TAU/nx  * ix*jx) - 1.0i * (TAU/nt * it*kt))*fn[kt][jx];
				return c;
			});
			res.push_back(DiscreteComplexFunction(v, vec2(0, TAU/Lx)));
		}
		return DiscreteComplexFunctionR2(res, vec2(0, TAU/Lt));
	}


	DiscreteComplexFunctionR2 ifft() const {
		float wt=domain[1];
		float wx=fn[0].getDomain()[1];
		vector<DiscreteComplexFunction> res = vector<DiscreteComplexFunction>();

		for (int it=0; it<samples_t(); ++it) {
			auto v = Vector<Complex>(samples_x(), [fn=fn, nx=samples_x(), nt=samples_t(), it](int ix) {
				Complex c = 0;
				for (int jx=0; jx<nx; ++jx)
					for (int kt=0; kt<nt; ++kt)
						c += exp(1.0i * (TAU/nx  * ix*jx) + 1.0i * (TAU/nt * it*kt))*fn[kt][jx];
				return c/(nt*nx);
			});
			res.push_back(DiscreteComplexFunction(v, vec2(0, TAU/wx)));
		}
		return DiscreteComplexFunctionR2(res, vec2(0, TAU/wt));
	}

	DiscreteRealFunctionR2 re();
	DiscreteRealFunctionR2 im();
	DiscreteRealFunctionR2 abs();
	DiscreteRealFunctionR2 arg();

};


class DiscreteRealFunctionR2 {
	vector<DiscreteRealFunction> fn;
	vec2 domain;
public:
	DiscreteRealFunctionR2(const vector<DiscreteRealFunction> &fn, vec2 domain) : fn(fn), domain(domain) {}
	DiscreteRealFunctionR2(const vector<HOM(float, float)> &f, vec2 domain, int sampling) {
		fn = vector<DiscreteRealFunction>();
		for (int i=0; i<f.size(); i++) {
			fn.emplace_back(f[i], domain, sampling);
		}
	}
	DiscreteRealFunctionR2(const DiscreteRealFunctionR2 &other) = default;
	DiscreteRealFunctionR2(DiscreteRealFunctionR2 &&other) noexcept = default;
	DiscreteRealFunctionR2 & operator=(const DiscreteRealFunctionR2 &other) = default;
	DiscreteRealFunctionR2 & operator=(DiscreteRealFunctionR2 &&other) noexcept = default;

	DiscreteRealFunction operator[](int t) const;
	DiscreteRealFunction operator()(float t) const;

	DiscreteRealFunctionR2 operator+(const DiscreteRealFunctionR2 &g) const;
	DiscreteRealFunctionR2 operator+(float a) const;
	friend DiscreteRealFunctionR2 operator+(float a, const DiscreteRealFunctionR2 &f);
	DiscreteRealFunctionR2 operator-() const;
	DiscreteRealFunctionR2 operator-(const DiscreteRealFunctionR2 &g) const;
	DiscreteRealFunctionR2 operator-(float a) const;
	friend DiscreteRealFunctionR2 operator-(float a, const DiscreteRealFunctionR2 &f);
	DiscreteRealFunctionR2 operator*(float a) const;
	friend DiscreteRealFunctionR2 operator*(float a, const DiscreteRealFunctionR2 &f);
	DiscreteRealFunctionR2 operator/(float a) const;
	friend DiscreteRealFunctionR2 operator/(float a, const DiscreteRealFunctionR2 &f);
	DiscreteRealFunctionR2 operator/(const DiscreteRealFunctionR2 &g) const;



	DiscreteRealFunctionR2 operator*(const DiscreteRealFunctionR2 &g) const;

	int samples_t () const { return fn.size(); }
	int samples_x () const { return fn[0].samples(); }
	float sampling_step_t() const { return (domain[1] - domain[0]) / (samples_t()-1); }
	float sampling_step_x() const { return fn[0].sampling_step(); }
	vector<float> args_t() const { return linspace(domain[0], domain[1], samples_t()); }
	vector<float> args_x() const { return fn[0].args(); }

	DiscreteRealFunctionR2 transpose() const;
	DiscreteRealFunction integrate_t() const {return transpose().integrate_x(); }
	DiscreteRealFunction integrate_x() const {return DiscreteRealFunction(Vector<float>(fn.size(), [&](int i) { return fn[i].integral(); }), domain); }
	float double_integral() const { return integrate_x().integral(); }
	float double_integral(int i_t) const { return integrate_x().integral(i_t); }
	float double_integral(int i_t, int i_x) const;

	DiscreteComplexFunctionR2 complexify() const {
		vector<DiscreteComplexFunction> res = vector<DiscreteComplexFunction>();
		for (auto &f : fn) {
			res.push_back(DiscreteComplexFunction(f, domain, f.samples()));
		}
		return DiscreteComplexFunctionR2(res, domain);
	}

	DiscreteComplexFunctionR2 fft_t () const {
		return complexify().fft_t();
	}

	DiscreteComplexFunctionR2 fft_x () const {
		return complexify().fft_x();
	}

	DiscreteComplexFunctionR2 fft() const {
		return complexify().fft();
	}

	DiscreteRealFunctionR2 double_convolve(const DiscreteRealFunctionR2 &kernel) const {
		return (this->fft()*(kernel.fft())).ifft().re();
	}
};
