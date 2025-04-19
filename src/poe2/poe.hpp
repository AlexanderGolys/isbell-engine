#pragma once
#include "../fundamentals/func.hpp"

class Percent {
	float value;

public:
	explicit Percent(float value) : value(value) {}
	Percent() : value(0) {}
	Percent operator+(const Percent &other) const { return Percent(value + other.value); }
	Percent operator-(const Percent &other) const { return Percent(value - other.value); }
	Percent operator*(float other) const { return Percent(value * other); }
	Percent operator/(float other) const { return Percent(value / other); }
	Percent operator+(float other) const { return Percent(value + other); }
	Percent operator-(float other) const { return Percent(value - other); }

	void operator+=(const Percent &other) { value += other.value; }
	void operator-=(const Percent &other) { value -= other.value; }
	void operator+=(float other) { value += other; }
	void operator-=(float other) { value -= other; }
	void operator*=(float other) { value *= other; }
	void operator/=(float other) { value /= other; }
	std::ostream &operator<<(std::ostream &os) const {
		os << value*100 << "%";
		return os;
	}
	string str() const {
		return std::to_string(value * 100) + "%";
	}
};

class DamageRange {
	int min_v;
	int max_v;
public:
	DamageRange(int min, int max) : min_v(min), max_v(max) {}
	DamageRange() : min_v(0), max_v(0) {}
	DamageRange operator+(const DamageRange &other) const { return DamageRange(min_v + other.min_v, max_v + other.max_v); }
	DamageRange operator-(const DamageRange &other) const { return DamageRange(min_v - other.min_v, max_v - other.max_v); }
	DamageRange operator+(float other) const { return DamageRange(min_v + other, max_v + other); }
	DamageRange operator-(float other) const { return DamageRange(min_v - other, max_v - other); }
	DamageRange operator*(float other) const { return DamageRange(min_v * other, max_v * other); }
	DamageRange operator/(float other) const { return DamageRange(min_v / other, max_v / other); }
	void operator+=(const DamageRange &other) { min_v += other.min_v; max_v += other.max_v; }
	void operator-=(const DamageRange &other) { min_v -= other.min_v; max_v -= other.max_v; }
	void operator+=(float other) { min_v += other; max_v += other; }
	void operator-=(float other) { min_v -= other; max_v -= other; }
	void operator*=(float other) { min_v *= other; max_v *= other; }
	void operator/=(float other) { min_v /= other; max_v /= other; }
	std::ostream &operator<<(std::ostream &os) const {  os << min_v << "-" << max_v; return os; }
	string str() const { return std::to_string(min_v) + "-" + std::to_string(max_v); }
	int min() const { return min_v; }
	int max() const { return max_v; }
	int avg() const { return (min_v + max_v) / 2; }
};

enum ModifierType {
	FLAT, INCREASE, MORE
};

class Modifier {
public:
	vector<string> stats_affected;
	float value;
	ModifierType type;

	virtual ~Modifier() = default;
	Modifier(const vector<string> &stats_affected, float value, ModifierType type) : stats_affected(stats_affected), value(value), type(type) {}
	Modifier(const string &stat_affected, float value, ModifierType type) : stats_affected(vector{stat_affected}), value(value), type(type) {}
};

class FlatModifier : public Modifier {
	FlatModifier(const vector<string> &stats_affected, float value) : Modifier(stats_affected, value, FLAT) {}
	FlatModifier(const string &stat_affected, float value) : Modifier(stat_affected, value, FLAT) {}
};

class IncreaseModifier : public Modifier {
	IncreaseModifier(const vector<string> &stats_affected, float value) : Modifier(stats_affected, value, INCREASE) {}
	IncreaseModifier(const string &stat_affected, float value) : Modifier(stat_affected, value, INCREASE) {}
};

class MoreModifier : public Modifier {
	MoreModifier(const vector<string> &stats_affected, float value) : Modifier(stats_affected, value, MORE) {}
	MoreModifier(const string &stat_affected, float value) : Modifier(stat_affected, value, MORE) {}
};


template<typename U=float>
class Stat {
	U flat;
	float increase_coef = 1;
	float more_coef = 1;
	string name;

	public:
	virtual ~Stat() = default;
	Stat(U flat, const string &name) : flat(flat), name(name) {}
	Stat(U flat, float increase, float more, const string &name) : flat(flat), increase_coef(increase), more_coef(more), name(name) {}

	Stat add(U value) const { return Stat(flat + value, increase_coef, more_coef, name); }
	Stat increase(float increase) const { return Stat(flat, increase_coef + increase, more_coef, name); }
	Stat reduce(float reduce) const { return Stat(flat, increase_coef - reduce, more_coef, name); }
	Stat more(float more) const { return Stat(flat, increase_coef, more_coef*(more + 1), name); }
	Stat less(float less) const { return Stat(flat, increase_coef, more_coef*(1.f - less), name); }

	Stat modify(const Modifier &mod) const {
		if (mod.stats_affected.empty() || std::ranges::find(mod.stats_affected, name) != mod.stats_affected.end())
			return *this;

		if (mod.type == FLAT)
			return add(mod.value);

		if (mod.type == INCREASE)
			return increase(mod.value);

		return more(mod.value);
	}

	U value() const { return flat * increase_coef * more_coef; }
};


enum Tag {
	SPELL, ELEMENTAL, COLD, FIRE, LIGHTNING, PHYSICAL, CHAOS, AILMENT, DAMAGING_AILMENT, NON_DAMAGING_AILMENT, HIT, AOE
};


class Item {
};
