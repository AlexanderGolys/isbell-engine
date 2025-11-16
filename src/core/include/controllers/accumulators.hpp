#pragma once
#include "func.hpp"
#include "mat.hpp"

class Distribution {
	CONST_PROPERTY(RealFunction, pdf_func);
	CONST_PROPERTY(RealFunction, cdf_func);
	CONST_PROPERTY(float, norm);
	CONST_PROPERTY(vec2, support);
public:
	Distribution(const RealFunction& pdf, const RealFunction& cdf, vec2 support) : pdf_func(pdf), cdf_func(cdf), norm(cdf(support.y)), support(support) {}

	float cdf(float x) const;
	float pdf(float x) const;

	Distribution operator+(const Distribution& other) const;
	Distribution operator-(const Distribution& other) const;
	Distribution operator-() const;
	Distribution operator*(float a) const;
	Distribution operator/(float a) const;
	Distribution normalize() const;
	Distribution shift(float a) const;

	static sptr<Distribution> fromCDF(END(float) cdf, vec2 support, float eps=.01f);
	static sptr<Distribution> fromPDF(END(float) pdf, vec2 support, float eps=.01f, uint prec=10000);

	static sptr<Distribution> normal(float mean, float sigma, float tail_cutoff_sigmas=3);
	static sptr<Distribution> uniform(float a, float b);
	static sptr<Distribution> exponential(float lambda, float cutoff);
	static sptr<Distribution> gamma(float shape, float scale, float cutoff, uint prec=1000);
	static sptr<Distribution> inverseGamma(float shape, float scale, float cutoff, uint prec=1000);


};

class ImpulseResponseAccumulator {
	Distribution distribution;
	float current_t;
public:
	explicit ImpulseResponseAccumulator(const Distribution& distribution) : distribution(distribution), current_t(distribution.get_support().x) {}
	ImpulseResponseAccumulator(const Distribution& distribution, float start_t) : distribution(distribution), current_t(start_t) {}
	float step(float dt);
	bool isFinished() const;
};

class DefferedResponseAccumulator {
	vector<ImpulseResponseAccumulator> impulses;
public:
	DefferedResponseAccumulator() = default;
	void addImpulse(const Distribution& distribution);
	void addImpulse(const Distribution& distribution, float start_t);
	float step(float dt);
};

