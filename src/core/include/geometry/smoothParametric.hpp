#pragma once

#include "complexGeo.hpp"


class AffinePlane;
class SmoothParametricPlaneCurve;
class SmoothParametricSurface;
class AffineLine;
class RealFunctionPS;
class Differential1FormPS;
class Differential2FormPS;

class SmoothParametricCurve {
protected:
	Foo13 _f;
	Foo13 _df;
	Foo13 _ddf;
	std::function<Foo13(int)> _der_higher = [this](int n) {
		return n == 1 ? _df : n == 2 ? _ddf : derivativeOperator(_der_higher(n - 1), this->eps);
	};
	float eps;
	float t0;
	float t1;
	bool periodic;
	PolyGroupID id;

public:
	SmoothParametricCurve(Foo13 f, Foo13 df, Foo13 ddf, PolyGroupID id = DFLT_CURV, float t0 = 0, float t1 = TAU, bool periodic = true, float epsilon = 0.01f);
	SmoothParametricCurve(Foo13 f, Foo13 df, PolyGroupID id = DFLT_CURV, float t0 = 0, float t1 = TAU, bool periodic = true, float epsilon = 0.01f);

	explicit SmoothParametricCurve(const Foo13& f, PolyGroupID id = DFLT_CURV, float t0 = 0, float t1 = TAU, bool periodic = true, float epsilon = 0.01f);

	SmoothParametricCurve(const Foo13& f, vec2 dom, PolyGroupID id = DFLT_CURV, bool periodic = false, float epsilon = 0.01f);

	SmoothParametricCurve(const RealFunction& fx, const RealFunction& fy, const RealFunction& fz, PolyGroupID id = DFLT_CURV, float t0 = 0, float t1 = TAU, bool periodic = true,
						  float epsilon = 0.01f);
	SmoothParametricCurve(Foo13 f, std::function<Foo13(int)> derivativeOperator, PolyGroupID id = DFLT_CURV, float t0 = 0, optional<float> t1 = TAU, bool periodic = true,
						  float epsilon = 0.01f);
	SmoothParametricCurve(Foo13 f, std::vector<Foo13> derivatives, PolyGroupID id = DFLT_CURV, float t0 = 0, float t1 = TAU, bool periodic = true, float epsilon = 0.01f);
	SmoothParametricCurve(const SmoothParametricCurve& other);
	SmoothParametricCurve(SmoothParametricCurve&& other) noexcept;
	SmoothParametricCurve& operator=(const SmoothParametricCurve& other);
	SmoothParametricCurve& operator=(SmoothParametricCurve&& other) noexcept;


	PolyGroupID getID() const;

	void setID(PolyGroupID id);

	void copyID(const SmoothParametricCurve& other);

	SmoothParametricCurve shift(vec3 v) const;

	SmoothParametricCurve rotate(vec3 axis, float angle, vec3 center = ORIGIN_R3) const;

	SmoothParametricCurve scale(float a, vec3 center = ORIGIN_R3) const;


	vec3 derivative(float t) const;

	vec3 df(float t) const;

	vec3 operator()(float t) const;

	vec2 bounds() const;

	float getT0() const;

	float getT1() const;

	vec3 second_derivative(float t) const;

	vec3 higher_derivative(float t, int n) const;

	vec3 ddf(float t) const;

	vec3 tangent(float t) const;
	vec3 binormal(float t) const;
	vec3 normal(float t) const;
	float length(float t0, float t1, int n) const;

	SmoothParametricCurve precompose(SpaceEndomorphism g_) const;
	void precomposeInPlace(SpaceEndomorphism g);

	mat3 FrenetFrame(float t) const;
	float curvature(float t) const;

	float curvature_radius(float t) const;

	float torsion(float t) const;
	vec3 curvature_vector(float t) const;

	float speed(float t) const;

	bool isPeriodic() const;

	float getEps() const;

	AffinePlane osculatingPlane(float t) const;

	SmoothParametricSurface surfaceOfRevolution(const AffineLine& axis) const;
	SmoothParametricSurface screwMotion(float speed, int iterations) const;
	SmoothParametricSurface cylinder(vec3 direction, float h) const;
	SmoothParametricSurface pipe(float radius, bool useFrenetFrame = true) const;
	SmoothParametricSurface pipe(float radius, float eps, bool useFrenetFrame = true) const;
	SmoothParametricSurface canal(const std::function<float(float)>& radius) const;

	static SmoothParametricCurve constCurve(vec3 v);
	static SmoothParametricCurve span(vec3 p1, vec3 p2);
};


class AffineLine : public SmoothParametricCurve {
	vec3 p0, v; // p0 + tv
public:
	AffineLine(vec3 p0, vec3 v);
	static AffineLine spanOfPts(vec3 p0, vec3 p1);
	float distance(vec3 p) const;
	RealFunctionR3 distanceField() const;
	bool contains(vec3 p, float eps = 1e-6) const;
	vec3 orthogonalProjection(vec3 p) const;
	vec3 pivot() const;
	vec3 direction() const;
	float distance(AffineLine& l) const;
	AffineLine operator+(vec3 v) const;
	SmoothParametricSurface tube(float radius, float t0, float t1) const;
};


// R2 -> R3
class SmoothParametricSurface {
	Foo113 _f;
	Foo113 _df_t;
	Foo113 _df_u;

	CONST_PROPERTY(vec2, rangeT);
	CONST_PROPERTY(vec2, rangeU);
	CONST_PROPERTY(bool, periodicT);
	CONST_PROPERTY(bool, periodicU);
	PROPERTY(float, epsilon);

	PolyGroupID id = 420;
	SmoothParametricCurve precompose(const SmoothParametricPlaneCurve& c, PolyGroupID id = 420) const;
	SmoothParametricSurface postcompose(const SpaceEndomorphism& g) const;
	SmoothParametricSurface precompose(const PlaneSmoothEndomorphism& g, vec2 t_bounds, vec2 u_bounds) const;
	SmoothParametricSurface precompose(const PlaneAutomorphism& g) const;

public:
	SmoothParametricSurface(const Foo113& f, const Foo113& df_t, const Foo113& df_u, vec2 t_range, vec2 u_range, bool t_periodic = false, bool u_periodic = false,float epsilon = .01);
	SmoothParametricSurface(const Foo113& f, vec2 t_range, vec2 u_range, bool t_periodic = false, bool u_periodic = false, float epsilon = .01);
	SmoothParametricSurface(const std::function<SmoothParametricCurve(float)>& pencil, vec2 t_range, vec2 u_range, bool t_periodic = false, bool u_periodic = false,float eps = .01);
	SmoothParametricSurface(RealFunctionR2 plot, vec2 t_range, vec2 u_range);

	vec3 operator()(float t, float s) const;
	vec3 operator()(vec2 tu) const;

	vec3 parametersNormalised(vec2 tu) const;
	vec3 parametersNormalised(float t, float u) const;

	SmoothParametricSurface operator+(const SmoothParametricSurface& S) const;
	SmoothParametricSurface operator*(float a) const;
	SmoothParametricSurface operator-(const SmoothParametricSurface& S) const;
	SmoothParametricSurface normaliseParameters() const;
	friend SmoothParametricCurve operator &(const SmoothParametricSurface& S, const SmoothParametricPlaneCurve& c);
	friend SmoothParametricSurface operator &(const SpaceEndomorphism& f, const SmoothParametricSurface& S);
	friend SmoothParametricSurface operator &(const SmoothParametricSurface& S, const PlaneAutomorphism& c);

	SmoothParametricSurface shift(vec3 v) const;
	SmoothParametricSurface rotate(vec3 axis, float angle, vec3 center = ORIGIN_R3) const;
	SmoothParametricSurface scale(float a, vec3 center = ORIGIN_R3) const;

	SmoothParametricCurve restrictToInterval(vec2 p0, vec2 p1, PolyGroupID id = 420) const;
	SmoothParametricCurve restrictToInterval(vec2 p0, vec2 p1, bool periodic, PolyGroupID id = 420) const;

	SmoothParametricCurve constT(float t0) const;
	SmoothParametricCurve constU(float u0) const;

	float t0() const;
	float t1() const;
	float u0() const;
	float u1() const;
	float periodT() const;
	float periodU() const;

	vec3 normal(float t, float s) const;
	vec3 normal(vec2 tu) const;
	void changeDomain(vec2 t_range, vec2 u_range, bool t_periodic, bool u_periodic);

	mat2 metricTensor(float t, float s) const;
	EuclideanSpace<vec2, mat2> toTangentStdBasis(float t, float s, vec3 v) const;
	vec3 embeddTangentVector(float t, float s, EuclideanSpace<vec2, mat2> v) const;
	mat3 DarbouxFrame(float t, float s) const;
	mat2x3 tangentStandardBasis(float t, float s) const;

	mat2 changeOfTangentBasisToPrincipal(float t, float s) const;
	mat2 changeOfPrincipalBasisToStandard(float t, float s) const;
	mat2x3 tangentSpacePrincipalBasis(float t, float s) const;

	mat2 firstFundamentalForm(float t, float s) const;
	mat2 secondFundamentalForm(float t, float s) const;
	float _E(float t, float u) const;
	float _F(float t, float u) const;
	float _G(float t, float u) const;
	float _L(float t, float u) const;
	float _M(float t, float u) const;
	float _N(float t, float u) const;
	vec3 d2f_tt(float t, float u) const;
	vec3 d2f_uu(float t, float u) const;
	vec3 d2f_tu(float t, float u) const;

	vec3 normalCurvature(float t, float s, vec2 v) const;
	mat2 shapeOperator(float t, float s) const;
	float meanCurvature(float t, float s) const;
	float gaussianCurvature(float t, float s) const;
	mat2x3 principalDirections(float t, float s) const;
	std::pair<float, float> principalCurvatures(float t, float s) const;
	void normaliseDomainToI2();
	vec3 Laplacian(float t, float s) const;

	float meanCurvature(vec2 tu) const;
	float globalAreaIntegral(const RealFunctionPS& f) const;
	float DirichletFunctional() const;
	float biharmonicFunctional() const;

	SmoothParametricSurface meanCurvatureFlow(float dt) const;
};

SmoothParametricSurface ruledSurfaceJoinT(const SmoothParametricCurve& c1, const SmoothParametricCurve& c2, float u0 = 0, float u1 = 1);
SmoothParametricSurface ruledSurfaceJoinU(const SmoothParametricCurve& c1, const SmoothParametricCurve& c2, float t0 = 0, float t1 = 1);

SmoothParametricSurface ruledSurfaceJoinT(const SmoothParametricCurve& c1, const SmoothParametricCurve& c2, vec2 bounds);

SmoothParametricSurface ruledSurfaceJoinU(const SmoothParametricCurve& c1, const SmoothParametricCurve& c2, vec2 bounds);

SmoothParametricSurface bilinearSurface(const vec3& p00, const vec3& p01, const vec3& p10, const vec3& p11, vec2 t_range, vec2 u_range);
SmoothParametricSurface CoonsPatch(const SmoothParametricCurve& cDown, const SmoothParametricCurve& cLeft, const SmoothParametricCurve& cUp, const SmoothParametricCurve& cRight);
SmoothParametricSurface cone(const SmoothParametricCurve& c, float h, float r);
SmoothParametricSurface polarCone(const SmoothParametricCurve& r, vec3 center);


SmoothParametricSurface CoonsPatchDisjoint(const SmoothParametricCurve& c1, const SmoothParametricCurve& c2);


class SurfaceParametricPencil {
	HOM(float, SmoothParametricSurface) pencil;

public:
	explicit SurfaceParametricPencil(const HOM(float, SmoothParametricSurface)& pencil);

	explicit SurfaceParametricPencil(const BIHOM(float, vec2, vec3)& foo, vec2 range_t = vec2(0, TAU), vec2 range_u = vec2(0, TAU), float eps = .01);

	SmoothParametricSurface operator()(float t) const;

	vec3 operator()(float t, float u, float s) const;

	vec3 operator()(float t, vec2 us) const;
};


class CurveParametricPencil {
	HOM(float, SmoothParametricCurve) pencil;

public:
	explicit CurveParametricPencil(const HOM(float, SmoothParametricCurve)& pencil);

	explicit CurveParametricPencil(BIHOM(float, float, vec3) foo, vec2 bounds = vec2(0, 1), float eps = .01);

	SmoothParametricCurve operator()(float t) const;

	vec3 operator()(float t, float u) const;
};


class ParametricSurfaceFoliation {
	HOM(float, SmoothParametricCurve) pencil_of_leaves;
	vec2 pencil_domain;
	bool pencil_periodic;
	vector<SmoothParametricCurve> special_leaves;

public:
	explicit ParametricSurfaceFoliation(HOM(float, SmoothParametricCurve) pencil, vec2 pencil_domain = vec2(0, 1), bool periodic = true,
										const vector<SmoothParametricCurve>& special_leaves = {});

	SmoothParametricCurve getLeaf(float t) const;

	SmoothParametricCurve getSpecialLeaf(int i) const;

	SmoothParametricSurface getFoliatedSurface() const;
	vector<SmoothParametricCurve> sampleLeaves(int res) const;

	vector<SmoothParametricCurve> getSpecialLeaves() const;

	vec2 getDomain() const;
};


class Differential1FormPS;

class RealFunctionPS {
	Foo111 _f;
	Foo112 _df;
	shared_ptr<SmoothParametricSurface> surface;

public:
	RealFunctionPS(const std::function<float(float, float)>& f, const shared_ptr<SmoothParametricSurface>& surface);

	RealFunctionPS(const Foo21& f, const shared_ptr<SmoothParametricSurface>& surface);

	RealFunctionPS(const Foo31& emb_pullback, const shared_ptr<SmoothParametricSurface>& surface);

	float operator()(float t, float s) const { return _f(t, s); }

	static RealFunctionPS constant(float c, const shared_ptr<SmoothParametricSurface>& surface);

	RealFunctionPS constant(float c) const;
	RealFunctionPS operator*(float a) const;
	RealFunctionPS operator/(float a) const;
	RealFunctionPS operator-() const;
	RealFunctionPS operator+(const RealFunctionPS& f) const;
	RealFunctionPS operator-(const RealFunctionPS& f) const;
	RealFunctionPS operator*(const RealFunctionPS& f) const;
	RealFunctionPS operator/(const RealFunctionPS& f) const;
	friend RealFunctionPS operator/(float a, const RealFunctionPS& g);
	friend RealFunctionPS operator+(float a, const RealFunctionPS& g);
	friend RealFunctionPS operator-(float a, const RealFunctionPS& g);
	friend RealFunctionPS operator*(float a, const RealFunctionPS& g);

	Differential1FormPS df();
};

template <typename V>
class Linear1Form2D {
	V v1, v2;
	vec2 coefs;

public:
	Linear1Form2D(vec2 omega, V basis1, V basis2);

	float operator()(V v) const;

	Linear1Form2D operator*(float a) const;

	Linear1Form2D operator/(float a) const;

	Linear1Form2D operator+(const Linear1Form2D& other) const;

	Linear1Form2D operator-(const Linear1Form2D& other) const;

	static Linear1Form2D dual(V v, V v1, V v2);

	static Linear1Form2D dx(V x, V y);

	static Linear1Form2D dy(V x, V y);

	static std::pair<Linear1Form2D, Linear1Form2D> basisForms(V v1, V v2);

	vec2 localCoefs() const;
};

template <typename V>
class Linear2Form2D {
	V v1, v2;
	float coef;

public:
	Linear2Form2D(float c, V basis1, V basis2);

	float operator()(V a, V b) const;

	Linear2Form2D operator*(float a) const;

	Linear2Form2D operator/(float a) const;

	Linear2Form2D operator+(const Linear2Form2D& other) const;

	Linear2Form2D operator-(const Linear2Form2D& other) const;

	float localCoef() const;
};


class Differential2FormPS;

class Differential1FormPS {
	std::function<Linear1Form2D<vec3>(float, float)> _omega;
	shared_ptr<SmoothParametricSurface> surface;

public:
	Differential1FormPS(const std::function<Linear1Form2D<vec3>(float, float)>& omega, const shared_ptr<SmoothParametricSurface>& surface);

	Differential1FormPS(const std::function<Linear1Form2D<vec3>(vec2)>& omega, const shared_ptr<SmoothParametricSurface>& surface);

	Differential1FormPS(const std::function<Linear1Form2D<vec3>(vec3)>& emb_pullback, const shared_ptr<SmoothParametricSurface>& surface);

	Linear1Form2D<vec3> operator()(float t, float s) const;

	Linear1Form2D<vec3> operator()(vec2 tu) const;

	float operator()(float t, float s, vec3 v) const;

	float operator()(vec2 tu, vec3 v) const;

	Differential1FormPS operator*(float a) const;

	Differential1FormPS operator/(float a) const;

	Differential1FormPS operator-() const;

	Differential1FormPS operator+(const Differential1FormPS& eta) const;

	Differential1FormPS operator-(const Differential1FormPS& eta) const;

	Differential1FormPS operator*(const RealFunctionPS& f) const;

	Differential1FormPS operator/(const RealFunctionPS& f) const;

	// Differential1FormPS<glm::vec3> d() const { return Differential1FormPS<glm::vec3>([w=_omega](float t, float s) { return w(t, s).v1; }, surface); }
};

class Differential2FormPS {
	BIHOM(float, float, Linear2Form2D<vec3>) _omega;
	shared_ptr<SmoothParametricSurface> surface;

public:
	Differential2FormPS(const BIHOM(float, float, Linear2Form2D<vec3>)& omega, const shared_ptr<SmoothParametricSurface>& surface);

	Differential2FormPS(const HOM(vec2, Linear2Form2D<vec3>)& omega, const shared_ptr<SmoothParametricSurface>& surface);

	Differential2FormPS(const HOM(vec3, Linear2Form2D<vec3>)& emb_pullback, const shared_ptr<SmoothParametricSurface>& surface);
	Linear2Form2D<vec3> operator()(float t, float s) const;
	Linear2Form2D<vec3> operator()(vec2 tu) const;

	float operator()(float t, float s, vec3 v, vec3 w) const;

	float operator()(vec2 tu, vec3 v, vec3 w) const;

	Differential2FormPS operator*(float a) const;

	Differential2FormPS operator/(float a) const;

	Differential2FormPS operator-() const;

	Differential2FormPS operator+(const Differential2FormPS& eta) const;

	Differential2FormPS operator-(const Differential2FormPS& eta) const;

	Differential2FormPS operator*(const RealFunctionPS& f) const;

	Differential2FormPS operator/(const RealFunctionPS& f) const;
};

class VectorFieldPS {
	std::function<vec3(float, float)> _f_dt;
	std::function<vec3(float, float)> _f_ds;
	shared_ptr<SmoothParametricSurface> surface;

public:
	VectorFieldPS(const std::function<vec3(float, float)>& f_dt, const std::function<vec3(float, float)>& f_du, const shared_ptr<SmoothParametricSurface>& surface);

	vec3 operator()(float t, float s) const;

	VectorFieldPS operator*(float a) const;

	VectorFieldPS operator/(float a) const;

	VectorFieldPS operator+(const VectorFieldPS& v) const;

	VectorFieldPS operator-(const VectorFieldPS& v) const;

	VectorFieldPS operator-() const;

	VectorFieldPS operator*(const RealFunctionPS& f) const;

	VectorFieldPS operator/(const RealFunctionPS& f) const;

	vec2 deform(vec2 tu, vec2 dv) const;

	vec3 shiftAmbient(vec2 tu, vec2 dv) const;

	//	vector<vec3> generate_trajectory(vec2 tu, vec2 v, int n, float h) const { return RK4([v, h, s=surface](float t, vec3 w) { return (*s)(t*v); }, 0, (*surface)(v), h).solution(n); }
};


class FunctionalPartitionOfUnity {
	std::vector<Fooo> _F_i;

public:
	explicit FunctionalPartitionOfUnity(const std::vector<Fooo>& F_i);

	Fooo operator[](int i) const;

	int size() const;
};


Fooo BernsteinPolynomial(int n, int i, float t0, float t1);

Fooo BernsteinPolynomial(int n, int i);

FunctionalPartitionOfUnity BernsteinBasis(int n);
FunctionalPartitionOfUnity BernsteinBasis(int n, float t0, float t1);

inline Fooo BSpline(int i, int k, const std::vector<float>& knots);
FunctionalPartitionOfUnity BSplineBasis(int n, int k, const std::vector<float>& knots);
vector<float> uniformKnots(int n, int k);

SmoothParametricCurve freeFormCurve(const FunctionalPartitionOfUnity& family, const std::vector<vec3>& controlPts, vec2 domain, float eps = 0.0001);

SmoothParametricCurve BezierCurve(const std::vector<vec3>& controlPoints, float t0 = 0, float t1 = 1, float eps = .001);

SmoothParametricCurve BSplineCurve(const std::vector<vec3>& controlPoints, const std::vector<float>& knots, int k, float eps = .001);

SmoothParametricCurve BSplineCurve(const std::vector<vec3>& controlPoints, float t0, float t1, int k, float eps = .001);

SmoothParametricSurface freeFormSurface(const FunctionalPartitionOfUnity& F_i, const FunctionalPartitionOfUnity& G_i, const std::vector<std::vector<vec3>>& controlPts,
										vec2 range_t, vec2 range_u, float eps = 0.01);

SmoothParametricSurface BezierSurface(const std::vector<std::vector<vec3>>& controlPoints, float t0 = 0, float t1 = 1, float u0 = 0, float u1 = 1, float eps = .01);

SmoothParametricSurface BSplineSurface(const std::vector<std::vector<vec3>>& controlPoints, const std::vector<float>& knots_t, const std::vector<float>& knots_u, int k = 3,
									   float eps = .01);

SmoothParametricSurface BSplineSurfaceUniform(const std::vector<std::vector<vec3>>& controlPoints, vec2 t_range, vec2 u_range, int k = 3, float eps = .01);
