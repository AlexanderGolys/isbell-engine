#pragma once

#include <utility>

#include "planarGeometry.hpp"

using std::make_shared, std::shared_ptr, std::function, std::vector, std::pair, std::make_unique;


// const Mob CayleyTransform = Mob(1, -1.0i, 1, 1.0i);
// const Mob CayleyTransformInv = ~Mob(1, -1.0i, 1, -1.0i);
// const Mob Imob = Mob(1, 0, 0, 1);

class ComplexCurve // TODO: make this shit modern
{
public:
    HOM(float, Complex) _f;
    float epsilon = .01;
    HOM(float, Complex) _df;
    HOM(float, Complex) N;
    bool cyclic;
    float period;
    float t0, t1;

	explicit ComplexCurve(HOM(float, Complex) f, bool cyclic=true, float t0=0, float t1=TAU);
    explicit ComplexCurve(const SmoothParametricPlaneCurve& curve);
	Complex operator()(float t) const {return _f(t);}
    vector<Complex> sample(float x0, float x1, int n) {
		vector<Complex> res ={};
		for (int i=0; i<n; i++)
			res.push_back((*this)(lerp(x0, x1, 1.f*i/n)));
		return res;
	}
    vector<Complex> sample(int n) {
	    return sample(t0, t1, n);
    }
    ComplexCurve disjointUnion(ComplexCurve &other);

    static ComplexCurve line(Complex z0, Complex z1);
    static ComplexCurve circle(Complex z0, float r);
    static ComplexCurve arc(Complex center, Complex z0, Complex z1);
};


class Meromorphism {
public:
	endC _f; // TODO: is this shared ptr even making sense
	endC _df;

	Meromorphism(const Meromorphism &other);
	Meromorphism(Meromorphism &&other) noexcept;
	Meromorphism & operator=(const Meromorphism &other);
	Meromorphism & operator=(Meromorphism &&other) noexcept;
	Meromorphism() : _f([](Complex){return Complex(0.f);}), _df([](Complex){return Complex(0.f);}) {}
	Meromorphism(endC f, endC df) : _f(std::move(f)), _df(std::move(df)) {}
	explicit Meromorphism(endC f, float eps=.01f) : _f(std::move(f)), _df([F=f, e=eps](Complex z){return (F(z + Complex(e)) - F(z))/e;}) {}

	Complex operator()(Complex z) const {return _f(z);}
	Complex operator()(vec2 z) const { return _f(Complex(z)); }
	Complex df(Complex z) const { return _df(z); }
	Complex df(vec2 z) const { return _df(Complex(z)); }
	Meromorphism compose(const Meromorphism &g) const {return Meromorphism([f=_f, g=g._f](Complex z){return f(g(z));}, [f=_f, df=_df, g=g._f, dg=g._df](Complex z){return df(g(z)) * dg(z);});}
	Meromorphism operator&(const Meromorphism &g) const {return (this)->compose(g);}
	Meromorphism operator+(const Meromorphism& g) const {return Meromorphism([f=_f, g=g._f](Complex z){return f(z) + g(z);}, [df=_df, dg=g._df](Complex z){return df(z) + dg(z);});}
	Meromorphism operator*(const Meromorphism& g) const {return Meromorphism([f=_f, g=g._f](Complex z){return f(z) * g(z);}, [f=_f, df=_df, g=g._f, dg=g._df](Complex z){return f(z) * dg(z) + df(z) * g(z);});}
	Meromorphism operator-() const {return Meromorphism([f=_f](Complex z){return -f(z);}, [df=_df](Complex z){return -df(z);});}
	Meromorphism operator-(const Meromorphism& g) const {return *this + -g;}
};


class Biholomorphism : public Meromorphism {
public:
	Meromorphism f_inv;
	Biholomorphism(endC f, endC df, endC f_inv) : Meromorphism(std::move(f), std::move(df)), f_inv(std::move(f_inv)) {}
	Biholomorphism(endC f, endC f_inv, float eps=.01) : Meromorphism(std::move(f), eps), f_inv(std::move(f_inv), eps) {}
	Biholomorphism(Meromorphism f, Meromorphism f_inv) : Meromorphism(f), f_inv(f_inv) {}

	Complex inv(Complex z) const { return f_inv(z); }
	Biholomorphism operator~() const {return Biholomorphism(f_inv, f_inv._df, _f);}
	Biholomorphism operator*(Complex a) const;
	Biholomorphism operator+(Complex a) const;
	Biholomorphism operator-(Complex a) const;
	Biholomorphism operator/(Complex a) const;
	Biholomorphism inv() const { return ~(*this); }
	Complex operator()(Complex z) const { return _f(z); }
	Biholomorphism compose(Biholomorphism g) const;
	static Biholomorphism mobius(Matrix<Complex> m);
	static Biholomorphism linear(Complex a, Complex b);
	static Biholomorphism _LOG();
	static Biholomorphism _EXP();
	static Biholomorphism power(float a);
};







inline Biholomorphism Biholomorphism::operator*(Complex a) const {
	return linear(a, 0i).compose(*this);
}

inline Biholomorphism Biholomorphism::operator+(Complex a) const {
	return linear(1i, a).compose(*this);
}

inline Biholomorphism Biholomorphism::operator-(Complex a) const {
	return linear(1i, -a).compose(*this);
}

inline Biholomorphism Biholomorphism::operator/(Complex a) const {
	return linear(1/a, 0i).compose(*this);
}


inline Biholomorphism Biholomorphism::mobius(Matrix<Complex> m) {
	Complex a = m.at(0, 0);
	Complex b = m.at(0, 1);
	Complex c = m.at(1, 0);
	Complex d = m.at(1, 1);
	return Biholomorphism([a, b, c, d](Complex z) {return (a * z + b) / (c * z + d); },
		[a, b, c, d](Complex z) {return (a * d - b * c) / ((c * z + d) * (c * z + d)); },
		[a, b, c, d](Complex z) {return (d * z - b) / (-c * z + a); });
}
