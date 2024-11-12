# pragma once

#include "mat.hpp"
#include <functional>
#include <memory>
#include <iostream>
#include <optional>
#include <variant>
#include <string>



#define endC std::function<Complex(Complex)>
#define Fooo std::function<float(float)>
#define Foo12 std::function<glm::vec2(float)>
#define Foo13 std::function<glm::vec3(float)>
#define Foo33 std::function<glm::vec3(glm::vec3)>
#define Foo31 std::function<float(glm::vec3)>
#define Foo21 std::function<float(glm::vec2)>
#define Foo22 std::function<glm::vec2(glm::vec2)>
#define Foo32 std::function<glm::vec2(glm::vec3)>
#define Foo23 std::function<glm::vec3(glm::vec2)>
#define Foo113 std::function<glm::vec3(float, float)>
#define pencilCurv std::function<SmoothParametricCurve(float)>
#define pencilSurf std::function<SmoothParametricSurface(float)>
#define Foo3Foo33 std::function<glm::mat3(glm::vec3)>

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




enum Regularity {
	UNKNOWN = -100,
	PATHOLOGICAL=-99,
	DISCONTINUOUS=-1,
	C0 = 0,
	C1 = 1,
	C2 = 2,
	C3 = 3,
	C4 = 4,
	POSSIBLY_FINITE_REGULARITY_AT_LEAST_C5 = 5,
	SMOOTH = 69,
	CONFORMAL=70,
	MEROMORPHIC = 420,
	HOLOMORPHIC = 2137
};

typedef glm::vec2(*func_pt_R22)(glm::vec2);
typedef glm::vec2(*func_pt_R12)(float);

Fooo derivativeOperator(Fooo f, float epsilon=0.01);
Foo12 derivativeOperator(Foo12 f, float epsilon=0.01);
Foo13 derivativeOperator(Foo13 f, float epsilon=0.01);


class VectorFieldR2 {
public:
	Foo22 field;
	VectorFieldR2();
	explicit VectorFieldR2(Foo22 field);
	glm::vec2 operator()(glm::vec2 v) const;
};

class VectorFieldR3;

class SmoothRealFunctionR3 {
	Foo31 _f;
	Foo33 _df;
    float eps = 0.01;
public:
	SmoothRealFunctionR3();
	SmoothRealFunctionR3(Foo31 f, Foo33 df, float eps=.01) : _f(f), _df(df), eps(eps) {};
	SmoothRealFunctionR3(Foo31 f, float epsilon=0.01);
	float operator()(glm::vec3 v) const;
	glm::vec3 df(glm::vec3 v) const;

  static SmoothRealFunctionR3 linear(glm::vec3 v);
  static SmoothRealFunctionR3 constant(float a);
  static SmoothRealFunctionR3 projection(int i);


    SmoothRealFunctionR3 operator*(float a) const;
    SmoothRealFunctionR3 operator+(float a) const;
    SmoothRealFunctionR3 operator-(float a) const;
    SmoothRealFunctionR3 operator/(float a) const;
    SmoothRealFunctionR3 operator-() const { return *this * -1; }

	SmoothRealFunctionR3 operator+(const SmoothRealFunctionR3 &g) const;
	SmoothRealFunctionR3 operator-(const SmoothRealFunctionR3 &g) const { return *this + (-g); }
	SmoothRealFunctionR3 operator*(const SmoothRealFunctionR3 &g) const;
	SmoothRealFunctionR3 operator/(const SmoothRealFunctionR3 &g) const;
    friend SmoothRealFunctionR3 operator*(float a, const SmoothRealFunctionR3 &f) { return f * a; }
    friend SmoothRealFunctionR3 operator+(float a, const SmoothRealFunctionR3 &f) { return f + a; }
    friend SmoothRealFunctionR3 operator-(float a, const SmoothRealFunctionR3 &f) { return -f + a; }
    friend SmoothRealFunctionR3 operator/(float a, const SmoothRealFunctionR3 &f) { return constant(a)/f; }
    SmoothRealFunctionR3 operator~() const { return 1/(*this); }



    VectorFieldR3 gradient() const;
    SmoothRealFunctionR3 Laplacian() const;
    float dx (glm::vec3 x) const { return _df(x).x; }
    float dy (glm::vec3 x) const { return _df(x).y; }
    float dz (glm::vec3 x) const { return _df(x).z; }
};


// R3 -> R3
class SpaceEndomorphism {
protected:
	Foo33 _f;
	Foo3Foo33 _df;
public:
  SpaceEndomorphism(const SpaceEndomorphism &other) : _f(other._f), _df(other._df) {}
  SpaceEndomorphism(SpaceEndomorphism &&other) noexcept : _f(std::move(other._f)), _df(std::move(other._df)) {}
  SpaceEndomorphism(Foo33 f, std::function<glm::mat3(glm::vec3)> df);
  SpaceEndomorphism(Foo33 f, float epsilon=0.01);
  SpaceEndomorphism &operator=(const SpaceEndomorphism &other);
  SpaceEndomorphism &operator=(SpaceEndomorphism &&other) noexcept;
  SpaceEndomorphism(glm::mat3 A) : _f([A](glm::vec3 x) { return A * x; }), _df([A](glm::vec3 x) { return A; }) {}
  SpaceEndomorphism(glm::mat4 A) : _f([A](glm::vec3 x) { return glm::vec3(A * glm::vec4(x, 1)); }), _df([A](glm::vec3 x) { return glm::mat3(A); }) {}

	glm::vec3 directional_derivative(glm::vec3 x, glm::vec3 v) const { return _df(x) * v; }
    glm::vec3 dfdv(glm::vec3 x, glm::vec3 v) const { return directional_derivative(x, v); }
	glm::mat3 df(glm::vec3 x) const { return _df(x); }
	glm::vec3 operator()(glm::vec3 x) const { return _f(x); }
	SpaceEndomorphism compose(SpaceEndomorphism g) const;

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
    VectorFieldR3(Foo33 X, Foo3Foo33 dX, float eps=.01) : _X(X), _dX(dX), eps(eps) {}
    VectorFieldR3(Foo33 X, float eps=.01);
    VectorFieldR3(SmoothRealFunctionR3 Fx , SmoothRealFunctionR3 Fy, SmoothRealFunctionR3 Fz, float epsilon=0.01);



	explicit VectorFieldR3(Foo33 field);
	VectorFieldR3(VectorFieldR2 f);
	glm::vec3 operator()(glm::vec3 v) const { return _X(v); }
	VectorFieldR3 operator+(const VectorFieldR3 &Y) const;
	VectorFieldR3 operator*(float a) const;
    VectorFieldR3 operator-() const { return *this * -1; }
	VectorFieldR3 operator-(const VectorFieldR3 &Y) const { return *this + (-Y); }
	VectorFieldR3 operator*(const SmoothRealFunctionR3 &f) const;
    SmoothRealFunctionR3 F_x() const { return SmoothRealFunctionR3([this](glm::vec3 x) { return _X(x).x; }, [this](glm::vec3 x) { return _dX(x)[0]; }); }
    SmoothRealFunctionR3 F_y() const { return SmoothRealFunctionR3([this](glm::vec3 x) { return _X(x).y; }, [this](glm::vec3 x) { return _dX(x)[1]; }); }
    SmoothRealFunctionR3 F_z() const { return SmoothRealFunctionR3([this](glm::vec3 x) { return _X(x).z; }, [this](glm::vec3 x) { return _dX(x)[2]; }); }
    std::array<SmoothRealFunctionR3, 3> components() const { return {F_x(), F_y(), F_z()}; }

    friend VectorFieldR3 operator*(const glm::mat3 &A, const VectorFieldR3 &X) {
      return VectorFieldR3([f=X._X, A](glm::vec3 v) {return A * f(v); }, [df=X._dX, A](glm::vec3 v) {return A * df(v); }, X.eps);
    }

	static VectorFieldR3 constant(glm::vec3 v);
	static VectorFieldR3 linear(glm::mat3 A);
	static VectorFieldR3 radial(glm::vec3 scale);

    SmoothRealFunctionR3 divergence() const;
    VectorFieldR3 curl() const;

    glm::vec3 moveAlong(glm::vec3 v, float dt=1) const { return v + _X(v) * dt; }
};


inline VectorFieldR3 SmoothRealFunctionR3::gradient() const { return VectorFieldR3(_df, eps); }
inline SmoothRealFunctionR3 SmoothRealFunctionR3::Laplacian() const { return gradient().divergence(); }


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

    static SpaceAutomorphism linear(glm::mat3 A);
    static SpaceAutomorphism translation(glm::vec3 v);
    static SpaceAutomorphism scaling(float x, float y, float z);
    static SpaceAutomorphism affine(glm::mat3 A, glm::vec3 v);
    static SpaceAutomorphism rotation(float angle);
    static SpaceAutomorphism rotation(glm::vec3 axis, float angle);
};


class AffinePlane;

// R -> R3
class SmoothParametricCurve {
protected:
    Foo13 _f;
	Foo13 _df;
	Foo13 _ddf;
	std::function<Foo13(int)> _der_higher =
		[this](int n){return n == 1 ? _df : n == 2 ? _ddf : derivativeOperator(_der_higher(n-1), this->eps);};
	float eps;
    RP1 t0 = 0;
    RP1 t1 = TAU;
    bool periodic;
    PolyGroupID id;
public:
    SmoothParametricCurve(Foo13 f, Foo13 df, Foo13 ddf, PolyGroupID id=PolyGroupID(420),  RP1 t0=0, RP1 t1=TAU, bool periodic=true, float epsilon=0.01);
	SmoothParametricCurve(Foo13 f, Foo13 df,PolyGroupID id=PolyGroupID(420),  RP1 t0=0, RP1 t1=TAU, bool periodic=true, float epsilon=0.01);
    explicit SmoothParametricCurve(Foo13 f, PolyGroupID id=PolyGroupID(420),  RP1 t0=0, RP1 t1=TAU, bool periodic=true, float epsilon=0.01) : SmoothParametricCurve(f, derivativeOperator(f, epsilon), id, t0, t1, periodic, epsilon) {}

    SmoothParametricCurve(Foo13 f, std::function<Foo13(int)> derivativeOperator, PolyGroupID id=PolyGroupID(420), RP1 t0=0, RP1 t1=TAU, bool periodic=true, float epsilon=0.01);
    SmoothParametricCurve(Foo13 f, std::vector<Foo13> derivatives, PolyGroupID id=PolyGroupID(420), RP1 t0=0, RP1 t1=TAU, bool periodic=true, float epsilon=0.01);
    SmoothParametricCurve(const SmoothParametricCurve &other);
    SmoothParametricCurve(SmoothParametricCurve &&other) noexcept;
    SmoothParametricCurve &operator=(const SmoothParametricCurve &other);
    SmoothParametricCurve &operator=(SmoothParametricCurve &&other) noexcept;

    PolyGroupID getID() const { return id; }
    void setID(PolyGroupID id) { this->id = id; }
    void copyID (const SmoothParametricCurve &other) { this->id = other.id; }


	glm::vec3 derivative(float t) const { return _df(t); }
	glm::vec3 df(float t) const { return derivative(t); }
    glm::vec3 operator()(float t) const;
    glm::vec2 bounds() const { return glm::vec2(t0.value_or(-1), t1.value_or(1)); }

	glm::vec3 second_derivative(float t) const { return _ddf(t); }
	glm::vec3 higher_derivative(float t, int n) const { return _der_higher(n)(t); }
	glm::vec3 ddf(float t) const { return second_derivative(t); }
	glm::vec3 tangent(float t) const { return normalise(_df(t)); }
	glm::vec3 normal(float t) const { return cross(tangent(t), binormal(t)); }
	glm::vec3 binormal(float t) const { return normalise(cross(_df(t), _ddf(t))); }
	float length(float t0, float t1, int n) const;
	SmoothParametricCurve precompose(SpaceEndomorphism g_) const;
	void precomposeInPlace(SpaceEndomorphism g);
	glm::mat3 FrenetFrame(float t) const;
	float curvature(float t) const;
	float torsion(float t) const;
	glm::vec3 curvature_vector(float t) const;
	AffinePlane osculatingPlane(float t) const;
    float speed(float t) const { return norm(df(t));}
    bool isPeriodic() const { return periodic; }


	static SmoothParametricCurve constCurve(glm::vec3 v);
};





class AffineLine : public SmoothParametricCurve {
	glm::vec3 p0, v; // p0 + tv
public:
	AffineLine(glm::vec3 p0, glm::vec3 v);
	static AffineLine spanOfPts(glm::vec3 p0, glm::vec3 p1);
	float distance(glm::vec3 p) const;
	SmoothRealFunctionR3 distanceField() const;
	bool contains(glm::vec3 p, float eps=1e-6) const;
	glm::vec3 orthogonalProjection(glm::vec3 p) const;
	glm::vec3 pivot() const;
	glm::vec3 direction() const;
	float distance(AffineLine &l) const;
	AffineLine operator+(glm::vec3 v) const;
};








class SmoothParametricPlaneCurve {
    Foo12 _f;
    Foo12 _df;
    Foo12 _ddf;
    std::function<Foo12(int)> _der_higher = [this](int n)
            {return n == 0 ? _f : n == 1 ? _df : n == 2 ? _ddf : derivativeOperator(_der_higher(n-1), this->eps);};
    float eps = 0.01;
    RP1 t0 = std::nullopt;
    RP1 t1 = std::nullopt;
    bool periodic = true;
public:
	explicit SmoothParametricPlaneCurve(const Foo12& curve, float t0=0, float t1=TAU, bool period=true, float epsilon=0.01);
    SmoothParametricPlaneCurve(Foo12 f,const Foo12& df, float t0=0, float t1=TAU, bool period=true, float epsilon=0.01);
    SmoothParametricPlaneCurve(Foo12 f, Foo12 df, Foo12 ddf,
                               float t0 = 0, float t1 = TAU, bool period = true, float epsilon = 0.01);
    SmoothParametricPlaneCurve(const SmoothParametricPlaneCurve &other);
    SmoothParametricPlaneCurve(SmoothParametricPlaneCurve &&other) noexcept;
    SmoothParametricPlaneCurve &operator=(const SmoothParametricPlaneCurve &other);
    SmoothParametricPlaneCurve &operator=(SmoothParametricPlaneCurve &&other) noexcept;
    glm::vec2 operator()(float t) const { return _f(t); }
    glm::vec2 derivative(float t) const { return _df(t); }
    glm::vec2 df(float t) const { return derivative(t); }
    glm::vec2 second_derivative(float t) const { return _ddf(t); }
    glm::vec2 ddf(float t) const { return second_derivative(t); }
    glm::vec2 higher_derivative(float t, int n) const { return _der_higher(n)(t); }
    glm::vec2 tangent(float t) const { return normalise(_df(t)); }
    glm::vec2 normal(float t) const { return orthogonalComplement(tangent(t)); }
	std::vector<glm::vec2> sample(float t0, float t1, int n) const;
	std::vector<glm::vec2> sample(int n) const {return sample(t0.value_or(-1), t1.value_or(1), n);}
	std::vector<glm::vec3> adjacency_lines_buffer(float t0, float t1, int n, float z=0) const;
    SmoothParametricCurve embedding(glm::vec3 v1=e1, glm::vec3 v2=e2, glm::vec3 pivot=ORIGIN) const;
    glm::vec2 bounds() const { return glm::vec2(t0.value_or(-1), t1.value_or(1)); }
  bool isPeriodic() const { return periodic; }
};

class ComplexCurve // TODO: make this shit modern
{
public:
	std::unique_ptr<std::function<Complex(float)>> f;
	float epsilon;
	std::unique_ptr<std::function<Complex(float)>> df;
	std::unique_ptr<std::function<Complex(float)>> ddf;
	std::unique_ptr<std::function<Complex(float)>> N;
	bool cyclic;
	float period;
	float t0, t1;

	ComplexCurve(std::function<Complex(float)> curve, float t0, float t1, float period=0.f, float epsilon = 0.01);
	ComplexCurve(SmoothParametricPlaneCurve* curve);
	Complex operator()(float t) const;
	std::vector<Complex> sample(float t0, float t1, int n);
	std::vector<Complex> sample(int n);
	ComplexCurve disjointUnion(ComplexCurve &other);

	static ComplexCurve line(Complex z0, Complex z1);
	static ComplexCurve circle(Complex z0, float r);
	static ComplexCurve arc(Complex center, Complex z0, Complex z1);
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

class Meromorphism {
public:
	std::shared_ptr<endC> _f; // TODO: is this shared ptr even making sense
	std::shared_ptr<endC> _df;
	Meromorphism();
	Meromorphism(std::shared_ptr<endC> f, std::shared_ptr<endC> df);
	Meromorphism(std::shared_ptr<endC> f, float eps); // todo
	Complex operator()(Complex z) const;
	Complex df(Complex z) const;
	operator PlaneSmoothEndomorphism() const; // todo
	Meromorphism compose(Meromorphism g) const;
	Meromorphism operator+(Meromorphism g) const;
	Meromorphism operator*(Meromorphism g) const;
	Meromorphism operator-() const;
	Meromorphism operator-(Meromorphism g) const;
	Meromorphism operator/(Meromorphism g) const;
};


class Biholomorphism : public Meromorphism {
public:
	std::shared_ptr<endC> _f_inv;
	Biholomorphism();
	Biholomorphism(std::shared_ptr<endC> f, std::shared_ptr<endC> df, std::shared_ptr<endC> f_inv);
	Biholomorphism(std::shared_ptr<endC> f, std::shared_ptr<endC> f_inv, float eps);
	Complex f_inv(Complex z) const;
	Complex inv(Complex z) const { return f_inv(z); }
	Biholomorphism operator~() const;
	Biholomorphism operator*(Complex a) const;
	Biholomorphism operator+(Complex a) const;
	Biholomorphism operator-(Complex a) const;
	Biholomorphism operator/(Complex a) const;
	Biholomorphism inv() const { return ~(*this); }
	operator PlaneAutomorphism() const;
	Complex operator()(Complex z) const;
	Biholomorphism compose(Biholomorphism g) const;
	static Biholomorphism mobius(Matrix<Complex, 2> mobius);
	static Biholomorphism linear(Complex a, Complex b);
	static Biholomorphism _LOG();
	static Biholomorphism _EXP();
	static Biholomorphism power(float a);

};



inline Biholomorphism Biholomorphism::operator*(Complex a) const {
	return Biholomorphism::linear(a, 0).compose(*this);
}

inline Biholomorphism Biholomorphism::operator+(Complex a) const {
	return Biholomorphism::linear(ONE, a).compose(*this);
}

inline Biholomorphism Biholomorphism::operator-(Complex a) const {
	return Biholomorphism::linear(ONE, -a).compose(*this);
}

inline Biholomorphism Biholomorphism::operator/(Complex a) const {
	return Biholomorphism::linear(ONE/a, 0).compose(*this);
}

inline Complex Biholomorphism::operator()(Complex z) const {
	return (*_f)(z);
}

template<typename T>
void printVector(std::vector<T> v, std::string title="vector")
{
	std::cout << title << " [" << v.size() << "]: ";
	for (int i = 0; i < v.size(); i++)
	{
		std::cout << v[i] << ", ";
	}
	std::cout << std::endl;
}





