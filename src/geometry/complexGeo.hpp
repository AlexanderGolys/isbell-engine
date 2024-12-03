#pragma once

#include "planarGeometry.hpp"


const Mob CayleyTransform = Mob(ONE, -I, ONE, I);
const Mob CayleyTransformInv = ~Mob(ONE, -I, ONE, I);
const Mob Imob = Mob(1, ZERO, ZERO, 1);

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
    // explicit ComplexCurve(SmoothParametricPlaneCurve* curve);
    Complex operator()(float t) const;
    std::vector<Complex> sample(float t0, float t1, int n);
    std::vector<Complex> sample(int n);
    ComplexCurve disjointUnion(ComplexCurve &other);

    static ComplexCurve line(Complex z0, Complex z1);
    static ComplexCurve circle(Complex z0, float r);
    static ComplexCurve arc(Complex center, Complex z0, Complex z1);
};


class Meromorphism {
public:
	std::shared_ptr<endC> _f; // TODO: is this shared ptr even making sense
	std::shared_ptr<endC> _df;
	Meromorphism();
	Meromorphism(std::shared_ptr<endC> f, std::shared_ptr<endC> df);
	Meromorphism(std::shared_ptr<endC> f, float eps); // todo
	Complex operator()(Complex z) const;
	Complex operator()(vec2 z) const { return (*_f)(Complex(z)); }
	Complex df(Complex z) const;
	Complex df(vec2 z) const { return (*_df)(Complex(z)); }
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
