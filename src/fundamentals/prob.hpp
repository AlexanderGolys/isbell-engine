#pragma once
#include <iostream>
#include <iomanip>
#include <string>
#include <map>
#include <random>
#include "macros.hpp"

class RealDistribution {
	public:
	virtual ~RealDistribution() = default;
	virtual float operator()() { return 0; }
};

class UniformRealDistribution : public RealDistribution {
	float a, b;
	std::mt19937 gen;
	std::uniform_real_distribution<float> dist;

public:
	UniformRealDistribution(float a, float b) : a(a), b(b), gen(std::random_device{}()), dist(a, b) {}

	float operator()() override {
		return dist(gen);
	}
};


class R2Distribution {
	public:
	virtual ~R2Distribution() = default;
	virtual vec2 operator()() = 0;
};

class ProductDistribution : public R2Distribution {
	RealDistribution p_x, p_y;
public:
	ProductDistribution(const RealDistribution &p_x, const RealDistribution &p_y) : p_x(p_x), p_y(p_y) {}
	vec2 operator()() override {
		return vec2(p_x(), p_y());
	}
	vector<vec2> operator()(int n) {
		vector<vec2> res;
		for (int i = 0; i < n; i++)
			res.push_back((*this)());
		return res;
	}
};
