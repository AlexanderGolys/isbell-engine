#include "accumulators.hpp"



float Distribution::cdf(float x) const {
	return cdf_func(x);
}

float Distribution::pdf(float x) const {
	return pdf_func(x);
}

Distribution Distribution::operator+(const Distribution& other) const {
	return Distribution(pdf_func + other.pdf_func, cdf_func + other.cdf_func, vec2(min(support.x, other.support.x), max(support.y, other.support.y)));
}

Distribution Distribution::operator*(float a) const {
	THROW_IF(a < 0, "Cannot scale distribution by negative number");
	return Distribution(pdf_func * a, cdf_func * a, support);
}

Distribution Distribution::operator/(float a) const {
	THROW_IF(a <= 0, "Cannot scale distribution by non-positive number");
	return Distribution(pdf_func / a, cdf_func / a, support);
}

Distribution Distribution::normalize() const {
	return *this / norm;
}



sptr<Distribution> Distribution::fromCDF(END(float) cdf, vec2 support, float eps) {
	RealFunction cdf_func = RealFunction(cdf, eps);
	RealFunction pdf_func = cdf_func.df();
	return make_shared<Distribution>(pdf_func, cdf_func, support);
}

sptr<Distribution> Distribution::fromPDF(END(float) pdf, vec2 support, float eps, uint prec) {
	RealFunction pdf_func = RealFunction(pdf, eps);
	RealFunction cdf_func = cdf_func.antiderivative(support.x, prec);
	return make_shared<Distribution>(pdf_func, cdf_func, support);
}

float ImpulseResponseAccumulator::step(float dt) {
	if (current_t >= distribution.get_support().y)
		return 0;
	current_t += dt;
	return distribution.cdf(current_t) - distribution.cdf(current_t-dt);
}

bool ImpulseResponseAccumulator::isFinished() const {
	return current_t >= distribution.get_support().y;
}

void DefferedResponseAccumulator::addImpulse(const Distribution& distribution) {
	impulses.emplace_back(distribution);
}

void DefferedResponseAccumulator::addImpulse(const Distribution& distribution, float start_t) {
	impulses.emplace_back(distribution, start_t);
}

float DefferedResponseAccumulator::step(float dt) {
	float total = 0.f;
	for (auto& impulse : impulses)
		total += impulse.step(dt);
	std::erase_if(impulses, [](const ImpulseResponseAccumulator& impulse) { return impulse.isFinished(); });
	return total;
}
