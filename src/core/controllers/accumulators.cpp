#include "accumulators.hpp"

#include "logging.hpp"


float Distribution::cdf(float x) const {
	return cdf_func(x);
}

float Distribution::pdf(float x) const {
	return pdf_func(x);
}

Distribution Distribution::operator+(const Distribution& other) const {
	return Distribution(pdf_func + other.pdf_func, cdf_func + other.cdf_func, vec2(min(support.x, other.support.x), max(support.y, other.support.y)));
}

Distribution Distribution::operator-(const Distribution& other) const {
	return *this + (other * -1.f);
}

Distribution Distribution::operator-() const {
	return *this * -1.f;
}

Distribution Distribution::operator*(float a) const {
	return Distribution(pdf_func * a, cdf_func * a, support);
}

Distribution Distribution::operator/(float a) const {
	return Distribution(pdf_func / a, cdf_func / a, support);
}

Distribution Distribution::normalize() const {
	return *this / norm;
}

Distribution Distribution::shift(float a) const {
	return Distribution(pdf_func & (X_R - a), cdf_func & (X_R - a), support + vec2(a, a));
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

sptr<Distribution> Distribution::normal(float mean, float sigma, float tail_cutoff_sigmas) {
	auto pdf = RealFunction([mean, sigma](float x) {
		return std::exp(-pow2(x-mean)/(2*pow2(sigma))) / (sigma * std::sqrt(TAU));
	});
	auto cdf = RealFunction([mean, sigma](float x) {
		return 0.5f * (1 + erf((x - mean) / (sigma * std::sqrt(2))));
	});
	return make_shared<Distribution>(pdf, cdf, vec2(mean - tail_cutoff_sigmas * sigma, mean + tail_cutoff_sigmas * sigma));
}

sptr<Distribution> Distribution::uniform(float a, float b) {
	auto pdf = RealFunction([a, b](float x) {
		return (x < a || x > b) ? 0.f : 1.f / (b - a);
	});
	auto cdf = RealFunction([a, b](float x) {
		if (x < a) return 0.f;
		if (x > b) return 1.f;
		return (x - a) / (b - a);
	});
	return make_shared<Distribution>(pdf, cdf, vec2(a, b));
}

sptr<Distribution> Distribution::exponential(float lambda, float cutoff) {
	auto pdf = RealFunction([lambda](float x) {
		return x < 0 ? 0.f : lambda * std::exp(-lambda * x);
	});
	auto cdf = RealFunction([lambda](float x) {
		return x < 0 ? 0.f : 1.f - std::exp(-lambda * x);
	});
	return make_shared<Distribution>(pdf, cdf, vec2(0.f, cutoff));
}

sptr<Distribution> Distribution::gamma(float shape, float scale, float cutoff, uint prec) {
	THROW_IF(shape <= 0, ValueError, "Shape parameter must be positive");
	THROW_IF(scale <= 0, ValueError, "Scale parameter must be positive");
	THROW_IF(cutoff <= 0, ValueError, "Cutoff must be positive");
	auto pdf = RealFunction([shape, scale](float x) {
		if (x <= 0) return 0.f;
		return powf(x, shape - 1) * std::exp(-x / scale) / (powf(scale, shape) * (float)tgamma(shape));
	});
	auto cdf = RealFunction([shape, scale, prec](float x) {
		if (x <= 0) return 0.f;
		return lowerIncompleteGamma(shape, x / scale, prec) / (float)tgamma(shape);
	});
	return make_shared<Distribution>(pdf, cdf, vec2(0.f, cutoff));
}

sptr<Distribution> Distribution::inverseGamma(float shape, float scale, float cutoff, uint prec) {
	THROW_IF(shape <= 0, ValueError, "Shape parameter must be positive");
	THROW_IF(scale <= 0, ValueError, "Scale parameter must be positive");
	THROW_IF(cutoff <= 0, ValueError, "Cutoff must be positive");
	auto pdf = RealFunction([shape, scale](float x) {
		if (x <= 0) return 0.f;
		return powf(scale, shape) * std::exp(-scale / x) / (powf(x, shape + 1) * (float)tgamma(shape));
	});
	auto cdf = RealFunction([shape, scale, prec](float x) {
		if (x <= 0) return 0.f;
		return upperIncompleteGamma(shape, scale / x, prec) / (float)tgamma(shape);
	});
	return make_shared<Distribution>(pdf, cdf, vec2(0.f, cutoff));
}

float ImpulseResponseAccumulator::step(float dt) {
	current_t += dt;
	if (current_t >= distribution.get_support().y)
		return 0;
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
	std::erase_if(impulses, [](const ImpulseResponseAccumulator& impulse) {
		return impulse.isFinished();
	});
	return total;
}

