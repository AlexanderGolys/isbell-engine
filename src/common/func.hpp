#ifndef FUNC_HPP
#define FUNC_HPP

#include "mat.hpp"
#include <functional>
#include <memory>
#include <iostream>

#include "func.hpp"

#define endC std::function<Complex(Complex)>

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

std::function<float(float)> derivativeOperator(std::function<float(float)> f, float epsilon=0.01);
std::function<glm::vec2(float)> derivativeOperator(std::function<glm::vec2(float)> f, float epsilon=0.01);
std::function<glm::vec3(float)> derivativeOperator(std::function<glm::vec3(float)> f, float epsilon=0.01);


class VectorFieldR2 {
public:
	std::function<glm::vec2(glm::vec2)> field;
	VectorFieldR2();
	explicit VectorFieldR2(std::function<glm::vec2(glm::vec2)> field);
	glm::vec2 operator()(glm::vec2 v) const;
};

class VectorFieldR3;

class SmoothRealFunctionR3 {
	std::function<float(glm::vec3)> _f;
	std::function<glm::vec3(glm::vec3)> _df;
public:
	SmoothRealFunctionR3();
	SmoothRealFunctionR3(std::function<float(glm::vec3)> f, std::function<glm::vec3(glm::vec3)> df);
	SmoothRealFunctionR3(std::function<float(glm::vec3)> f, float epsilon=0.01);
	float operator()(glm::vec3 v) const;
	glm::vec3 df(glm::vec3 v) const;
	SmoothRealFunctionR3 operator+(SmoothRealFunctionR3 g) const;
	SmoothRealFunctionR3 operator-(SmoothRealFunctionR3 g) const;
	SmoothRealFunctionR3 operator*(SmoothRealFunctionR3 g) const;
	SmoothRealFunctionR3 operator/(SmoothRealFunctionR3 g) const;
	SmoothRealFunctionR3 operator*(float a) const;
    SmoothRealFunctionR3 operator+(float a) const;
	SmoothRealFunctionR3 operator-(float a) const;
	SmoothRealFunctionR3 operator/(float a) const;
	VectorFieldR3 gradient() const;

	static SmoothRealFunctionR3 linear(glm::vec3 v);
	static SmoothRealFunctionR3 constant(float a);
	static SmoothRealFunctionR3 projection(int i);
};


class SpaceEndomorphism {
protected:
	std::function<glm::vec3(glm::vec3)> _f;
	std::function<glm::mat3(glm::vec3)> _df;
public:
	SpaceEndomorphism(std::function<glm::vec3(glm::vec3)> f, std::function<glm::mat3(glm::vec3)> df);
	SpaceEndomorphism(std::function<glm::vec3(glm::vec3)> f, float epsilon=0.01);
	glm::vec3 directional_derivative(glm::vec3 x, glm::vec3 v) const;
	glm::mat3 df(glm::vec3 x) const { return _df(x); }
	glm::vec3 operator()(glm::vec3 v) const;
	SpaceEndomorphism compose(SpaceEndomorphism g) const;

	static SpaceEndomorphism linear(glm::mat3 A);
	static SpaceEndomorphism translation(glm::vec3 v);
	static SpaceEndomorphism scaling(float x, float y, float z);
	static SpaceEndomorphism affine(glm::mat3 A, glm::vec3 v);
};

class VectorFieldR3 {
	std::function<glm::vec3(glm::vec3)> _field;
public:
	VectorFieldR3();
	explicit VectorFieldR3(std::function<glm::vec3(glm::vec3)> field);
	explicit VectorFieldR3(SpaceEndomorphism f);
	VectorFieldR3(VectorFieldR2 f);
	glm::vec3 operator()(glm::vec3 v) const;
	VectorFieldR3 operator+(VectorFieldR3 g) const;
	VectorFieldR3 operator*(float a) const;
	VectorFieldR3 operator*(glm::mat3 A) const;
	VectorFieldR3 operator-(VectorFieldR3 g) const;
	VectorFieldR3 operator*(SmoothRealFunctionR3 f) const;

	static VectorFieldR3 constant(glm::vec3 v);
	static VectorFieldR3 linear(glm::mat3 A);
	static VectorFieldR3 radial(glm::vec3 scale);
	static VectorFieldR3 rotational();
	static VectorFieldR3 wirlpool();
	static VectorFieldR3 gradient(SmoothRealFunctionR3 f);

};

const VectorFieldR3 dabbaX = VectorFieldR3::constant(glm::vec3(1, 0, 0));
const VectorFieldR3 dabbaY = VectorFieldR3::constant(glm::vec3(0, 1, 0));
const VectorFieldR3 dabbaZ = VectorFieldR3::constant(glm::vec3(0, 0, 1));



class SpaceAutomorphism : public SpaceEndomorphism {
	std::function<glm::vec3(glm::vec3)> _f_inv;
public:
	SpaceAutomorphism(std::function<glm::vec3(glm::vec3)> f, std::function<glm::vec3(glm::vec3)> f_inv, std::function<glm::mat3(glm::vec3)> df);
	SpaceAutomorphism(std::function<glm::vec3(glm::vec3)> f, std::function<glm::vec3(glm::vec3)> f_inv, float epsilon=0.01);

	glm::vec3 inv(glm::vec3 v) const;
	SpaceAutomorphism operator~() const;
	SpaceAutomorphism inv() const { return ~(*this); }
	SpaceAutomorphism compose(SpaceAutomorphism g) const;
};

class ParametricCurve {
protected:
	std::function<glm::vec3(float)> _f;
public:
	Regularity regularity=UNKNOWN;
	ParametricCurve();
	explicit ParametricCurve(std::function<glm::vec3(float)> curve);
	glm::vec3 operator()(float t) const;
};

class AffinePlane;

class SmoothParametricCurve : public ParametricCurve {
protected:
	std::function<glm::vec3(float)> _df;
	std::function<glm::vec3(float)> _ddf;
	std::function<std::function<glm::vec3(float)>(int)> _der_higher =
		[this](int n){return n == 1 ? _df : n == 2 ? _ddf : derivativeOperator(_der_higher(n-1), this->eps);};
	float eps;
public:
	SmoothParametricCurve();
	SmoothParametricCurve(std::function<glm::vec3(float)> f, std::function<glm::vec3(float)> df, std::function<glm::vec3(float)> ddf, float epsilon=0.01);
	SmoothParametricCurve(std::function<glm::vec3(float)> f, std::function<glm::vec3(float)> df, float epsilon=0.01);
	SmoothParametricCurve(std::function<glm::vec3(float)> f, float epsilon=0.01);
	SmoothParametricCurve(std::function<glm::vec3(float)> f, std::function<std::function<glm::vec3(float)>(int)> derivativeOperator, float epsilon=0.01);
	SmoothParametricCurve(std::function<glm::vec3(float)> f, std::vector<std::function<glm::vec3(float)>> derivatives, float epsilon=0.01);
	glm::vec3 derivative(float t) const { return _df(t); }
	glm::vec3 df(float t) const { return derivative(t); }
	glm::vec3 second_derivative(float t) const { return _ddf(t); }
	glm::vec3 higher_derivative(float t, int n) const { return _der_higher(n)(t); }
	glm::vec3 ddf(float t) const { return second_derivative(t); }
	glm::vec3 tangent(float t) const { return normalise(_df(t)); }
	glm::vec3 normal(float t) const { return cross(tangent(t), binormal(t)); }
	glm::vec3 binormal(float t) const { return normalise(cross(_df(t), _ddf(t))); }
	float length(float t0, float t1, int n) const;
	SmoothParametricCurve operator*(SpaceEndomorphism g) const;
	void operator*=(SpaceEndomorphism g);
	SmoothParametricCurve compose(SpaceEndomorphism g) const {return (*this)*g;}
	void composeInPlace(SpaceEndomorphism g) {(*this)*=g;}
	glm::mat3 FrenetFrame(float t) const;
	float curvature(float t) const;
	float torsion(float t) const;
	glm::vec3 curvature_vector(float t) const;
	AffinePlane osculatingPlane(float t) const;

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

class SmoothParametricSurface {
	std::function <glm::vec3(float, float)> _f;
	std::function <glm::vec3(float, float)> _df_t;
	std::function <glm::vec3(float, float)> _df_u;
	float t0=NAN, t1=NAN, u0=NAN, u1=NAN;
	bool t_periodic, u_periodic;
public:
	SmoothParametricSurface(std::function <glm::vec3(float, float)> f, std::function <glm::vec3(float, float)> df_t,
		std::function <glm::vec3(float, float)> df_u, float t0=0, float t1=PI, float u0=0, float u1=PI,
		bool t_periodic=true, bool u_periodic=true);

	SmoothParametricSurface(std::function <glm::vec3(float, float)> f, float t0=0, float t1=PI, float u0=0, float u1=PI,
		bool t_periodic=true, bool u_periodic=true, float epsilon=.01);

	SmoothParametricSurface(std::function<SmoothParametricCurve(float)> pencil, float t0=0, float t1=PI, float curve_t0=0,
	                        float curve_t1=PI, bool periodic=false, bool curve_periodic=true, float eps=.01);


	glm::vec3 operator()(float t, float s) const;
};

class SmoothImplicitSurface {
	SmoothRealFunctionR3 _F;
public:
	SmoothImplicitSurface(SmoothRealFunctionR3 F);
	float operator()(glm::vec3 p) const;
};

class AffinePlane : public SmoothImplicitSurface {
	glm::vec3 n;
	float d; // (n, x) - d = 0
	glm::vec3 pivot, v1, v2;
public:
	AffinePlane(glm::vec3 n, float d);
	AffinePlane(glm::vec3 pivot, glm::vec3 v1, glm::vec3 v2);
	static AffinePlane spanOfPts(glm::vec3 p0, glm::vec3 p1, glm::vec3 p2);
	AffineLine intersection(AffinePlane &p) const;
	float distance (glm::vec3 p) const;
	glm::vec3 orthogonalProjection(glm::vec3 p) const;
	bool contains(glm::vec3 p, float eps=1e-6) const;
	glm::vec3 normal() const;
	float getD() const;
	std::pair<glm::vec3, float> equationCoefs() const;
	glm::vec3 intersection(AffineLine &l) const;
	SmoothRealFunctionR3 distanceField() const;
	operator SmoothImplicitSurface() const;
	glm::mat3 pivotAndBasis() const;
	glm::vec2 localCoordinates(glm::vec3 p) const;
};







class ParametricPlanarCurve {
public:
	std::unique_ptr<std::function<glm::vec2(float)>> f;
	float epsilon;
	std::unique_ptr<std::function<glm::vec2(float)>> df;
	std::unique_ptr<std::function<glm::vec2(float)>> ddf;
	std::unique_ptr<std::function<glm::vec2(float)>> N;
	bool cyclic;
	float period;
	float t0, t1;

	ParametricPlanarCurve(std::function<glm::vec2(float)> curve, float t0=0, float t1=1, float period=0, float epsilon=0.01);
	glm::vec2 operator()(float t);
	std::vector<glm::vec2> sample(float t0, float t1, int n);
	std::vector<glm::vec2> sample(int n);
	std::vector<glm::vec3> adjacency_lines_buffer(float t0, float t1, int n, float z) const;
	
};

class ComplexCurve
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
	ComplexCurve(ParametricPlanarCurve* curve);
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
	std::shared_ptr<std::function<glm::vec2(glm::vec2)>> _f;
	std::shared_ptr<std::function<glm::mat2(glm::vec2)>> _df;
public:
	PlaneSmoothEndomorphism();
	PlaneSmoothEndomorphism(std::shared_ptr<std::function<glm::vec2(glm::vec2)>> f, std::shared_ptr<std::function<glm::mat2(glm::vec2)>> _df);
	glm::vec2 operator()(glm::vec2 x) const;
	glm::vec2 df(glm::vec2 x, glm::vec2 v) const;
	glm::mat2 df(glm::vec2 x) const;
};

class PlaneAutomorphism : public PlaneSmoothEndomorphism {
protected:
	std::shared_ptr<std::function<glm::vec2(glm::vec2)>> _f_inv;
public:
	PlaneAutomorphism(std::shared_ptr<std::function<glm::vec2(glm::vec2)>> f, std::shared_ptr<std::function<glm::mat2(glm::vec2)>> _df, std::shared_ptr<std::function<glm::vec2(glm::vec2)>> f_inv);
	glm::vec2 f_inv(glm::vec2 x) const;
	PlaneAutomorphism operator~() const;
	PlaneAutomorphism inv() const { return ~(*this); }
};

class Meromorphism {
public:
	std::shared_ptr<endC> _f;
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

const auto EXP = Biholomorphism::_EXP();
const auto LOG = Biholomorphism::_LOG();
const auto Id = Biholomorphism::linear(ONE, ZERO);
const auto ADD1 = Biholomorphism::linear(ONE, ONE);
const auto SQUARE = Biholomorphism::power(2);
const auto SQRT = Biholomorphism::power(.5f);
const auto CAYLEY = Biholomorphism::mobius(Matrix<Complex, 2>(ONE, -I, ONE, I));

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






#endif