#include "complexGeo.hpp"
#include <map>

using namespace glm;
using std::vector, std::string, std::shared_ptr, std::unique_ptr, std::pair, std::make_unique, std::make_shared;

ComplexCurve::ComplexCurve(std::function<Complex(float)> curve, float t0, float t1, float period, float epsilon)
{
	this->f = std::make_unique<std::function<Complex(float)>>(curve);
	this->cyclic = period > 0;
	this->epsilon = epsilon;
	this->t0 = t0;
	this->t1 = t1;

	this->df = std::make_unique<std::function<Complex(float)>>([this](float t) {
		return ((*this)(t + this->epsilon) - (*this)(t)) / this->epsilon;
	});

	this->ddf = std::make_unique<std::function<Complex(float)>>([this](float t) {
		return ((*this)(t + this->epsilon) - (*this)(t)*2.f + (*this)(t - this->epsilon)) / (this->epsilon * this->epsilon);
	});

	this->N = std::make_unique<std::function<Complex(float)>>([this](float t) {
		return Complex(normalise<vec2>(tvec2<float>((*ddf)(t)))); });

	this->period = period;
}


Complex ComplexCurve::operator()(float t) const {

    return (*f)(t);
}

std::vector<Complex> ComplexCurve::sample(float t0, float t1, int n)
{
	if (n == 1)
	{
		return { (*f)(t0) };
	}
	std::vector<Complex> result = std::vector<Complex>();
	result.reserve(n);
	float dt = (t1 - t0) / (n-1);
	for (int i = 0; i < n; i++)
	{
		result.push_back((*f)(t0 + i * dt));
	}
	printVector(result, "sampled complex curve from " + std::to_string(t0) + " to " + std::to_string(t1));
	return result;
}

std::vector<Complex> ComplexCurve::sample(int n)
{
	auto result = sample(t0, t1, n);
	return result;
}

ComplexCurve ComplexCurve::disjointUnion(ComplexCurve &other)
{
    return ComplexCurve([this, &other](float t) {
		if (t < this->t1)
		{
			return (*this->f)(t);
		}
		else
		{
			return (*other.f)(t - this->t1+other.t0);
		}
	}, this->t0, this->t1 + other.t1-other.t0, 0, this->epsilon);
}

ComplexCurve ComplexCurve::line(Complex z0, Complex z1)
{
	return ComplexCurve([&z0, &z1](float t) {
		return z0 + (z1 - z0) * t;
		}, 0, 1);
}

ComplexCurve ComplexCurve::circle(Complex z0, float r)
{
	return ComplexCurve([&z0, &r](float t) {
		return z0 + exp(1.0i * t) * r;
		}, 0, TAU, TAU, .001f);
}

ComplexCurve ComplexCurve::arc(Complex center, Complex z0, Complex z1)
{
	float r = abs(z0 - center);
	float phi0 = (z0 - center).arg();
	float phi1 = (z1 - center).arg();
	return ComplexCurve([&center, &r, &phi0, &phi1](float t) {
		return center + exp(1.0i * t) * r;
		}, phi0, phi1);
}

Meromorphism::Meromorphism() {
	_f = make_shared<endC>([](Complex z) {return z; });
	_df =  make_shared<endC>([](Complex z) {return 1; });
}

Meromorphism::Meromorphism(shared_ptr<endC> f, shared_ptr<endC> df) {
	_f = f;
	_df = df;
}


Meromorphism::operator PlaneSmoothEndomorphism() const {
	auto new_f = _f;
	auto new_df = _df;
	return PlaneSmoothEndomorphism(make_shared<std::function<vec2(vec2)>>([new_f](vec2 v) {return vec2((*new_f)(Complex(v))); }),
		make_shared<std::function<mat2(vec2) >>([new_df](vec2 x) {return mat2((*new_df)(Complex(x)).re(), (*new_df)(Complex(x)).im(), -(*new_df)(Complex(x)).im(), (*new_df)(Complex(x)).re()); }));
}

Meromorphism Meromorphism::compose(Meromorphism g) const {
	auto new_f = _f;
	auto new_df = _df;
	auto g_f = g._f;
	auto g_df = g._df;
	return Meromorphism(make_shared<endC>([new_f, g_f](Complex z) {return (*new_f)((*g_f)(z)); }),
		make_shared<endC>([new_f, new_df, g_f, g_df](Complex z) {return (*new_df)((*g_f)(z)) * (*g_df)(z); }));
}

Meromorphism Meromorphism::operator+(Meromorphism g) const {
	auto new_f = _f;
	auto new_df = _df;
	auto g_f = g._f;
	auto g_df = g._df;
	return Meromorphism(make_shared<endC>([new_f, g_f](Complex z) {return (*new_f)(z) + (*g_f)(z); }),
		make_shared<endC>([new_df, g_df](Complex z) {return (*new_df)(z) + (*g_df)(z); }));
}

Meromorphism Meromorphism::operator*(Meromorphism g) const {
	auto new_f = _f;
	auto new_df = _df;
	auto g_f = g._f;
	auto g_df = g._df;
	return Meromorphism(make_shared<endC>([new_f, g_f](Complex z) {return (*new_f)(z) * (*g_f)(z); }),
		make_shared<endC>([new_f, new_df, g_f, g_df](Complex z) {return (*new_df)(z) * (*g_f)(z) + (*new_f)(z) * (*g_df)(z); }));
}

Meromorphism Meromorphism::operator-() const {
	auto new_f = _f;
	auto new_df = _df;
	return Meromorphism(make_shared<endC>([new_f](Complex z) {return -(*new_f)(z); }),
		make_shared<endC>([new_df](Complex z) {return -(*new_df)(z); }));

}

Meromorphism Meromorphism::operator-(Meromorphism g) const {
	auto new_f = _f;
	auto new_df = _df;
	auto g_f = g._f;
	auto g_df = g._df;
	return Meromorphism(make_shared<endC>([new_f, g_f](Complex z) {return (*new_f)(z) - (*g_f)(z); }),
		make_shared<endC>([new_df, g_df](Complex z) {return (*new_df)(z) - (*g_df)(z); }));
}

Meromorphism Meromorphism::operator/(Meromorphism g) const {
	auto new_f = _f;
	auto new_df = _df;
	auto g_f = g._f;
	auto g_df = g._df;
	return Meromorphism(make_shared<endC>([new_f, g_f](Complex z) {return (*new_f)(z) / (*g_f)(z); }),
		make_shared<endC>([new_f, new_df, g_f, g_df](Complex z) {return ((*new_df)(z) * (*g_f)(z) - (*new_f)(z) * (*g_df)(z)) / ((*g_f)(z) * (*g_f)(z)); }));
}


Complex Meromorphism::operator()(Complex z) const {
	return (*_f)(z);
}

Complex Meromorphism::df(Complex z) const {
	return  (*_df)(z);
}


Biholomorphism::Biholomorphism() {
	_f = make_shared<endC>([](Complex z) {return z; });
	_df = make_shared<endC>([](Complex z) {return 1; });
	_f_inv = make_shared<endC>([](Complex z) {return z; });
}

Biholomorphism::Biholomorphism(shared_ptr<endC> f, shared_ptr<endC> df, shared_ptr<endC> f_inv) : Meromorphism(f, df) {
	_f_inv = f_inv;
}

Biholomorphism Biholomorphism::mobius(Matrix<Complex, 2> mobius) {
	auto _f = make_shared<endC>([mobius](Complex z) {return mobius.mobius(z); });
	auto _df = make_shared<endC>([mobius](Complex z) {return mobius.mobius_derivative(z); });
	auto _f_inv = make_shared<endC>([mobius](Complex z) {return (~mobius).mobius(z); });
	return Biholomorphism(_f, _df, _f_inv);
}

Biholomorphism Biholomorphism::linear(Complex a, Complex b) {
	return mobius(Matrix<Complex, 2>(a, b, 0, 1));
}

Biholomorphism Biholomorphism::_LOG() {
	auto _f = make_shared<endC>([](Complex z) {return log(z); });
	auto _df = make_shared<endC>([](Complex z) {return 1 / z; });
	auto _f_inv = make_shared<endC>([](Complex z) {return exp(z); });
	return Biholomorphism(_f, _df, _f_inv);
}

Biholomorphism Biholomorphism::_EXP() {
	auto _f = make_shared<endC>([](Complex z) {return exp(z); });
	auto _df = make_shared<endC>([](Complex z) {return exp(z); });
	auto _f_inv = make_shared<endC>([](Complex z) {return log(z); });
	return Biholomorphism(_f, _df, _f_inv);
}

Biholomorphism Biholomorphism::power(float a) {
	auto _f = make_shared<endC>([a](Complex z) {return z.pow(a); });
	auto _df = make_shared<endC>([a](Complex z) {return a!=0 ? z.pow(a - 1)*a : 0; });
	auto _f_inv = make_shared<endC>([a](Complex z) {return z.pow(1 / a); });
	return Biholomorphism(_f, _df, _f_inv);
}

Complex Biholomorphism::f_inv(Complex z) const {
    return (*_f_inv)(z);
}

Biholomorphism Biholomorphism::operator~() const {
    auto df_cpy = _df;
    auto f_inv_cpy = _f_inv;
    return Biholomorphism(_f_inv, make_shared<endC>([df_cpy, f_inv_cpy](Complex z) {return 1 / (*df_cpy)((*f_inv_cpy)(z)); }), _f);
}

Biholomorphism::operator PlaneAutomorphism() const {
    auto new_f = _f;
    auto new_df = _df;
    auto new_f_inv = _f_inv;
    return PlaneAutomorphism(make_shared<std::function<vec2(vec2)>>([new_f](vec2 v) {return vec2((*new_f)(Complex(v))); }),
        make_shared<std::function<mat2(vec2) >>([new_df](vec2 x) {return mat2((*new_df)(Complex(x)).re(), (*new_df)(Complex(x)).im(), -(*new_df)(Complex(x)).im(), (*new_df)(Complex(x)).re()); }),
        make_shared<std::function<vec2(vec2)>>([new_f_inv](vec2 v) {return vec2((*new_f_inv)(Complex(v))); }));
}

Biholomorphism Biholomorphism::compose(Biholomorphism g) const {
    auto new_f = _f;
    auto new_df = _df;
    auto new_inv = _f_inv;
    auto g_f = g._f;
    auto g_df = g._df;
    auto g_inv = g._f_inv;
    return Biholomorphism(make_shared<endC>([new_f, g_f](Complex z) {return (*new_f)((*g_f)(z)); }),
        make_shared<endC>([new_f, new_df, g_f, g_df](Complex z) {return (*new_df)((*g_f)(z)) * (*g_df)(z); }),
        make_shared<endC>([new_inv, g_inv](Complex z) {return (*g_inv)((*new_inv)(z)); }));
}

//ComplexCurve::operator SmoothParametricPlaneCurve()
//{
//	auto new_f = *_f;
//	return SmoothParametricPlaneCurve([new_f](float t) {
//		return vec2(new_f(t));
//	}, t0, t1, period, epsilon);
//}
