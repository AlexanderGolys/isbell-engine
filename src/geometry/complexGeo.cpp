#include "complexGeo.hpp"
#include <map>

using namespace glm;
using std::vector, std::string, std::shared_ptr, std::unique_ptr, std::pair, std::make_unique, std::make_shared;




ComplexCurve ComplexCurve::line(Complex z0, Complex z1)
{
	return ComplexCurve([&z0, &z1](float t) {
		return z0 + (z1 - z0) * t;
		}, 0, 1);
}

ComplexCurve ComplexCurve::circle(Complex z0, float r)
{
	return ComplexCurve([z0, r](float t) {
		return z0 + exp(1.0i * t) * r;
	});
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

Meromorphism::Meromorphism(const Meromorphism &other): _f(other._f),
													   _df(other._df) {}

Meromorphism::Meromorphism(Meromorphism &&other) noexcept: _f(std::move(other._f)),
														   _df(std::move(other._df)) {}

Meromorphism & Meromorphism::operator=(const Meromorphism &other) {
	if (this == &other)
		return *this;
	_f = other._f;
	_df = other._df;
	return *this;
}

Meromorphism & Meromorphism::operator=(Meromorphism &&other) noexcept {
	if (this == &other)
		return *this;
	_f = std::move(other._f);
	_df = std::move(other._df);
	return *this;
}

Biholomorphism Biholomorphism::compose(Biholomorphism g) const {
	return Biholomorphism([f=_f, g=g._f](Complex z) {
							  return f(g(z));
						  }, [f=_f, df=_df, g=g._f, dg=g._df](Complex z) {
							  return df(g(z)) * dg(z);
						  }, [f=f_inv, g=g.f_inv](Complex z) {
							  return f(g(z));
						  });
}

Biholomorphism Biholomorphism::linear(Complex a, Complex b) {
	return Biholomorphism([a, b](Complex z) {return a * z + b; }, [a, b](Complex z) {return a; }, [a, b](Complex z) {return (z - b) / a; });
}

Biholomorphism Biholomorphism::_LOG() {
	return Biholomorphism([](Complex z) {return log(z); }, [](Complex z) {return 1 / z; }, [](Complex z) {return exp(z); });
}

Biholomorphism Biholomorphism::_EXP() {
	return Biholomorphism([](Complex z) {return exp(z); }, [](Complex z) {return exp(z); }, [](Complex z) {return log(z); });
}

Biholomorphism Biholomorphism::power(float a) {
	return Biholomorphism([a](Complex z) {return pow(z, a); }, [a](Complex z) {return a * pow(z, a - 1); }, [a](Complex z) {return pow(z, 1 / a); });
}

ComplexCurve::ComplexCurve(std::function<Complex(float)> f, bool cyclic, float t0, float t1): _f(move(f)), cyclic(cyclic), t0(t0), t1(t1) {
	_df = [F=_f, eps=epsilon](float t){ return Complex(F(t+eps) - F(t))/eps;};
	N = [d=_df](float t){ return Complex(orthogonalComplement(d(t).z));};
	period = 0;
	if (cyclic) period = t1-t0;
}

ComplexCurve::ComplexCurve(const SmoothParametricPlaneCurve &curve)
: cyclic(false),
  t0(0),
  t1(0) {
	_f = [F=curve](float t) { return Complex(F(t)); };
	_df = [C=curve](float t) { return Complex(C.df(t)); };
	N = [d=_df](float t) { return Complex(orthogonalComplement(d(t).z)); };
	period = 0;
}
