#pragma once

// #include <utility>

#include "hyperbolic.hpp"
#include "glm/glm.hpp"
#include "src/fundamentals/macros.hpp"
// #include "planarGeometry.hpp"

class AffinePlane;
class SmoothParametricPlaneCurve;

class SmoothParametricCurve {
protected:
    Foo13 _f;
	Foo13 _df;
	Foo13 _ddf;
	std::function<Foo13(int)> _der_higher =
		[this](int n){return n == 1 ? _df : n == 2 ? _ddf : derivativeOperator(_der_higher(n-1), this->eps);};
	float eps;
    float t0;
    float t1;
    bool periodic;
    PolyGroupID id;
public:
    SmoothParametricCurve(Foo13 f, Foo13 df, Foo13 ddf, PolyGroupID id=DFLT_CURV, float t0=0, float t1=TAU, bool periodic=true, float epsilon=0.01);
	SmoothParametricCurve(Foo13 f, Foo13 df,PolyGroupID id=DFLT_CURV,  float t0=0, float t1=TAU, bool periodic=true, float epsilon=0.01);
    explicit SmoothParametricCurve(const Foo13 &f, PolyGroupID id=DFLT_CURV,  float t0=0, float t1=TAU, bool periodic=true, float epsilon=0.01) : SmoothParametricCurve(f, derivativeOperator(f, epsilon), id, t0, t1, periodic, epsilon) {}

    SmoothParametricCurve(Foo13 f, std::function<Foo13(int)> derivativeOperator, PolyGroupID id=DFLT_CURV, float t0=0, RP1 t1=TAU, bool periodic=true, float epsilon=0.01);
    SmoothParametricCurve(Foo13 f, std::vector<Foo13> derivatives, PolyGroupID id=DFLT_CURV, float t0=0, float t1=TAU, bool periodic=true, float epsilon=0.01);
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
    glm::vec2 bounds() const { return glm::vec2(t0, t1); }
    float getT0() const { return t0; }
    float getT1() const { return t1; }

	glm::vec3 second_derivative(float t) const { return _ddf(t); }
	glm::vec3 higher_derivative(float t, int n) const { return _der_higher(n)(t); }
	glm::vec3 ddf(float t) const { return second_derivative(t); }
	glm::vec3 tangent(float t) const { return normalise(_df(t)); }
	glm::vec3 normal(float t) const { return normalise(glm::cross(tangent(t), binormal(t))); }
	glm::vec3 binormal(float t) const { return normalise(cross(_df(t), _ddf(t))); }
	float length(float t0, float t1, int n) const;
	SmoothParametricCurve precompose(SpaceEndomorphism g_) const;
	void precomposeInPlace(SpaceEndomorphism g);

	glm::mat3 FrenetFrame(float t) const;
	float curvature(float t) const;
    float curvature_radius(float t) const { return 1/curvature(t); }
	float torsion(float t) const;
	glm::vec3 curvature_vector(float t) const;
    float speed(float t) const { return norm(df(t));}

    bool isPeriodic() const { return periodic; }
    float getEps() const { return eps; }
    AffinePlane osculatingPlane(float t) const;

	static SmoothParametricCurve constCurve(glm::vec3 v);
};





class AffineLine : public SmoothParametricCurve {
	glm::vec3 p0, v; // p0 + tv
public:
	AffineLine(glm::vec3 p0, glm::vec3 v);
	static AffineLine spanOfPts(glm::vec3 p0, glm::vec3 p1);
	float distance(glm::vec3 p) const;
	RealFunctionR3 distanceField() const;
	bool contains(glm::vec3 p, float eps=1e-6) const;
	glm::vec3 orthogonalProjection(glm::vec3 p) const;
	glm::vec3 pivot() const;
	glm::vec3 direction() const;
	float distance(AffineLine &l) const;
	AffineLine operator+(glm::vec3 v) const;
};



// R2 -> R3
class SmoothParametricSurface {
  Foo113 _f;
  Foo113 _df_t;
  Foo113 _df_u;
  float t0, t1, u0, u1;
  bool t_periodic, u_periodic;
  float epsilon=.01;

public:
  SmoothParametricSurface(Foo113 f, Foo113 df_t, Foo113 df_u, glm::vec2 t_range, glm::vec2 u_range, bool t_periodic=false, bool u_periodic=false);
  SmoothParametricSurface(Foo113 f, glm::vec2 t_range, glm::vec2 u_range, bool t_periodic=false, bool u_periodic=false, float epsilon=.01);
  SmoothParametricSurface(std::function<SmoothParametricCurve(float)> pencil, glm::vec2 t_range, glm::vec2 u_range, bool t_periodic=false, bool u_periodic=false, float eps=.01);
  SmoothParametricSurface(SmoothParametricCurve canal, Fooo width, float eps=.01);
  SmoothParametricSurface(SmoothParametricCurve pipe, float width, float eps=.01) : SmoothParametricSurface(pipe, [width](float t) {return width; }, eps) {}


  glm::vec3 operator()(float t, float s) const;
  glm::vec3 operator()(glm::vec2 tu) const;
  glm::vec3 parametersNormalised(glm::vec2 tu) const { return operator()(t0 + tu.x*(t1-t0), u0 + tu.y*(u1-u0)); }
  glm::vec3 parametersNormalised(float t, float u) const { return operator()(t0 + t*(t1-t0), u0 + u*(u1-u0)); }
  SmoothParametricCurve precompose(SmoothParametricPlaneCurve c, PolyGroupID id=420) const;
  SmoothParametricCurve restrictToInterval(glm::vec2 p0, glm::vec2 p1, PolyGroupID id=420) const;

  glm::vec2 boundsT() const { return glm::vec2(t0, t1); }
  glm::vec2 boundsU() const { return glm::vec2(u0, u1); }
  float tMin() const { return t0; }
  float tMax() const { return t1; }
  float uMin() const { return u0; }
  float uMax() const { return u1; }
  bool isPeriodicT() const { return t_periodic; }
  bool isPeriodicU() const { return u_periodic; }
  float periodT() const { return t_periodic ? t1 - t0 : 0; }
  float periodU() const { return u_periodic ? u1 - u0 : 0; }

  glm::vec3 normal(float t, float s) const;
  glm::vec3 normal(glm::vec2 tu) const { return normal(tu.x, tu.y); }
  glm::mat3 DarbouxFrame(float t, float s) const;
  glm::mat2x3 tangentSpace(float t, float s) const;
  void normaliseDomainToI2();
};

// template <typename T>
// std::function<T(float, float)> unpack(const std::function<T(glm::vec2)> &f) {
//     return [f](float t, float s) { return f(glm::vec2(t, s)); };
// }

// template <typename T>
// std::function<T(glm::vec2)> pack(const std::function<T(float, float)> &f) {
//     return [f](glm::vec2 tu) { return f(tu.x, tu.y); };
// }

class Differential1FormPS;

class RealFunctionPS {
    Foo111 _f;
    Foo112 _df;
    std::shared_ptr<SmoothParametricSurface> surface;
public:
    RealFunctionPS(const std::function<float(float, float)> &f, const std::shared_ptr<SmoothParametricSurface> &surface) : _f(f), surface(surface) {}
    RealFunctionPS(const Foo21 &f, const std::shared_ptr<SmoothParametricSurface> &surface) : _f(unpack(f, glm::vec2)), _df(unpack(derivativeOperator(f) vec2)), surface(surface) {}
    RealFunctionPS(const Foo31 &emb_pullback, const std::shared_ptr<SmoothParametricSurface> &surface);
    float operator()(float t, float s) const;

    static RealFunctionPS constant(float c, const std::shared_ptr<SmoothParametricSurface> &surface) { return RealFunctionPS([c](float, float) { return c; }, surface); }
    RealFunctionPS constant(float c) const { return constant(c, surface); }

    RealFunctionPS operator*(float a) const { return RealFunctionPS([f=_f, a](float t, float s) { return f(t, s)*a; }, surface); }
    RealFunctionPS operator/(float a) const { return (*this)*(1/a); }
    RealFunctionPS operator-() const { return (*this)*(-1); }
    RealFunctionPS operator+(const RealFunctionPS& f) const { return RealFunctionPS([f1=_f, f2=f._f](float t, float s) { return f1(t, s) + f2(t, s); }, surface); }
    RealFunctionPS operator-(const RealFunctionPS& f) const { return *this + (-f); }
    RealFunctionPS operator*(const RealFunctionPS& f) const { return RealFunctionPS([f1=_f, f2=f._f](float t, float s) { return f1(t, s) * f2(t, s); }, surface); }
    RealFunctionPS operator/(const RealFunctionPS& f) const { return RealFunctionPS([f1=_f, f2=f._f](float t, float s) { return f1(t, s) / f2(t, s); }, surface); }

    friend RealFunctionPS operator/(float a, const RealFunctionPS& g) { return g.constant(a)/g; }
    friend RealFunctionPS operator+(float a, const RealFunctionPS& g) { return g.constant(a)+g; }
    friend RealFunctionPS operator-(float a, const RealFunctionPS& g) { return g.constant(a)-g; }
    friend RealFunctionPS operator*(float a, const RealFunctionPS& g) { return g.constant(a)*g; }

    Differential1FormPS df();
};

template <typename V>
class Linear1Form2D {
    V v1, v2;
    glm::vec2 coefs;
public:
    Linear1Form2D(glm::vec2 omega, V basis1, V basis2) : v1(basis1), v2(basis2), coefs(omega) {}
    float operator()(V v) const { return glm::dot(glm::vec2(glm::dot(v, v1), glm::dot(v, v2)), coefs); }
    Linear1Form2D operator*(float a) const { return Linear1Form2D(coefs*a, v1, v2 ); }
    Linear1Form2D operator/(float a) const { return Linear1Form2D(coefs/a, v1, v2 ); }
    Linear1Form2D operator+(const Linear1Form2D &other) const { return Linear1Form2D(coefs + other.coefs, v1, v2); }
    Linear1Form2D operator-(const Linear1Form2D &other) const { return Linear1Form2D(coefs - other.coefs, v1, v2); }

    static Linear1Form2D dual(V v, V v1, V v2) { return Linear1Form2D(v1, v2, glm::vec2(glm::dot(v, v1), glm::dot(v, v2))); }
    static Linear1Form2D dx(V x, V y) { return Linear1Form2D(glm::vec2(1, 0), x, y); }
    static Linear1Form2D dy(V x, V y) { return Linear1Form2D(glm::vec2(0, 1), x, y); }
    static std::pair<Linear1Form2D, Linear1Form2D> basisForms (V v1, V v2) { return {dx(v1, v2), dy(v1, v2)}; }
    glm::vec2 localCoefs() const { return coefs; }
};

template <typename V>
class Linear2Form2D {
    V v1, v2;
    float coef;
public:
    Linear2Form2D(float c, V basis1, V basis2) : v1(basis1), v2(basis2), coef(c) {}
    float operator()(V a, V b) const { return coef*glm::dot(a, v1)*glm::dot(b, v2); }
    Linear2Form2D operator*(float a) const { return Linear2Form2D(coef*a, v1, v2 ); }
    Linear2Form2D operator/(float a) const { return Linear2Form2D(coef/a, v1, v2 ); }
    Linear2Form2D operator+(const Linear2Form2D &other) const { return Linear2Form2D(coef + other.coef, v1, v2); }
    Linear2Form2D operator-(const Linear2Form2D &other) const { return Linear2Form2D(coef - other.coef, v1, v2); }

    float localCoef() const { return coef; }
};


class Differential2FormPS;

class Differential1FormPS {
    std::function<Linear1Form2D<glm::vec3>(float, float)> _omega;
    std::shared_ptr<SmoothParametricSurface> surface;
public:
    Differential1FormPS(const std::function<Linear1Form2D<glm::vec3>(float, float)> &omega, const std::shared_ptr<SmoothParametricSurface> &surface) : _omega(omega), surface(surface) {}
    Differential1FormPS(const std::function<Linear1Form2D<glm::vec3>(vec2)> &omega, const std::shared_ptr<SmoothParametricSurface> &surface) : _omega(unpack(omega, vec2)), surface(surface) {}
    Differential1FormPS(const std::function<Linear1Form2D<glm::vec3>(vec3)> &emb_pullback, const std::shared_ptr<SmoothParametricSurface> &surface);
    Linear1Form2D<glm::vec3> operator()(float t, float s) const { return _omega(t, s); }
    Linear1Form2D<glm::vec3> operator()(glm::vec2 tu) const { return _omega(tu.x, tu.y); }
    float operator()(float t, float s, glm::vec3 v) const { return _omega(t, s)(v); }
    float operator()(glm::vec2 tu, glm::vec3 v) const { return pack(_omega)(tu)(v); }
    Differential1FormPS operator*(float a) const { return Differential1FormPS([w=_omega, a](float t, float s) { return w(t, s)*a; }, surface); }
    Differential1FormPS operator/(float a) const { return (*this)*(1/a); }
    Differential1FormPS operator-() const { return (*this)*(-1); }
    Differential1FormPS operator+(const Differential1FormPS& eta) const { return Differential1FormPS([w1=_omega, w2=eta._omega](float t, float s) { return w1(t, s) + w2(t, s); }, surface); }
    Differential1FormPS operator-(const Differential1FormPS& eta) const {return *this + (-eta); }
    Differential1FormPS operator*(const RealFunctionPS& f) const { return Differential1FormPS([w=_omega, f](float t, float s) { return w(t, s)*f(t, s); }, surface); }
    Differential1FormPS operator/(const RealFunctionPS& f) const { return Differential1FormPS([w=_omega, f](float t, float s) { return w(t, s)/f(t, s); }, surface); }

    // Differential1FormPS<glm::vec3> d() const { return Differential1FormPS<glm::vec3>([w=_omega](float t, float s) { return w(t, s).v1; }, surface); }
};

class Differential2FormPS {
    std::function<Linear2Form2D<glm::vec3>(float, float)> _omega;
    std::shared_ptr<SmoothParametricSurface> surface;
public:
    Differential2FormPS(const std::function<Linear2Form2D<glm::vec3>(float, float)> &omega, const std::shared_ptr<SmoothParametricSurface> &surface) : _omega(omega), surface(surface) {}
    Differential2FormPS(const std::function<Linear2Form2D<glm::vec3>(glm::vec2)> &omega, const std::shared_ptr<SmoothParametricSurface> &surface) : _omega([omega](float t, float s) { return omega(glm::vec2(t, s)); }), surface(surface) {}
    Differential2FormPS(const std::function<Linear2Form2D<glm::vec3>(glm::vec3)> &emb_pullback, const std::shared_ptr<SmoothParametricSurface> &surface);
    Linear2Form2D<glm::vec3> operator()(float t, float s) const;
    Linear2Form2D<glm::vec3> operator()(glm::vec2 tu) const;

    float operator()(float t, float s, glm::vec3 v, glm::vec3 w) const { return _omega(t, s)(v, w); }
    float operator()(glm::vec2 tu, glm::vec3 v, glm::vec3 w) const { return [om=_omega](glm::vec2 tu, glm::vec3 v, glm::vec3 w) { return om(tu.x, tu.y)(v, w); }(tu, v, w); }

    Differential2FormPS operator*(float a) const { return Differential2FormPS([w=_omega, a](float t, float s) { return w(t, s)*a; }, surface); }
    Differential2FormPS operator/(float a) const { return (*this)*(1/a); }
    Differential2FormPS operator-() const { return (*this)*(-1); }
    Differential2FormPS operator+(const Differential2FormPS& eta) const { return Differential2FormPS([w1=_omega, w2=eta._omega](float t, float s) { return w1(t, s) + w2(t, s); }, surface); }
    Differential2FormPS operator-(const Differential2FormPS& eta) const {return *this + (-eta); }
    Differential2FormPS operator*(const RealFunctionPS& f) const { return Differential2FormPS([w=_omega, f](float t, float s) { return w(t, s)*f(t, s); }, surface); }
    Differential2FormPS operator/(const RealFunctionPS& f) const { return Differential2FormPS([w=_omega, f](float t, float s) { return w(t, s)/f(t, s); }, surface); }



};
