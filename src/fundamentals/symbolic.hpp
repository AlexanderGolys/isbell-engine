#pragma once
#include <map>
#include <utility>
#include <vector>
#include <string>
#include "func.hpp"


using std::vector, std::string;



class Polynomial {
public:
	vector<float> coefs;
	void reduce_leading_zeros();
	std::multiset<Complex> roots = std::multiset<Complex>();
	std::set<int> term_degrees() const;
	void try_to_factorize();

	Polynomial();
	explicit Polynomial(float term, int monomial_pow=0);
	explicit Polynomial(const vector<float> &c);
	Polynomial(const Polynomial &p);
	Polynomial(Polynomial &&other) noexcept;
	Polynomial & operator=(const Polynomial &other);
	Polynomial & operator=(Polynomial &&other) noexcept;

	int degree() const;
	bool monomial() const;
	bool constant() const;
	float operator[](int i) const;
	void set(int i, float val);

	Polynomial operator+(const Polynomial &p) const;
	Polynomial operator-(const Polynomial &p) const;
	Polynomial operator*(const Polynomial &p) const;
	Polynomial operator*(float a) const;
	Polynomial operator/(float f) const;
	Polynomial operator+(float p) const;
	Polynomial operator-(float p) const;
	Polynomial operator-() const;

	friend Polynomial operator+(float f, const Polynomial &p);
	friend Polynomial operator-(float f, const Polynomial &p);
	friend Polynomial operator*(float f, const Polynomial &p);

	template<typename T> T operator()(const T& x) const;
	Polynomial operator&(const Polynomial &p) const { return (*this)(p); }

	bool operator==(const Polynomial &p) const { return coefs == p.coefs; }
	bool operator==(float c) const { return constant() && coefs[0] == c; }
};

template<typename T>
T Polynomial::operator()(const T& x) const {
	T c = T((*this)[0]);
	for (int i = 1; i < coefs.size(); i++)
		c += pow(x, i) * coefs[i];
	return c;
}

class Rational {
	Polynomial num, den;
public:
	explicit Rational(Polynomial num, Polynomial den=Polynomial(1)) : num(std::move(num)), den(std::move(den)) {}
	explicit Rational(float c) : num(c), den(1) {}

	int degree() const { return num.degree() - den.degree(); }
	bool constant() const { return num.constant() && den.constant() || num == 0; }
	float operator[](int i) const;
	void set(int i, float val);

	Rational operator+(const Rational &p) const { return Rational(num*p.den + p.num*den, den*p.den); }
	Rational operator-(const Rational &p) const { return Rational(num*p.den - p.num*den, den*p.den); }
	Rational operator*(const Rational &p) const { return Rational(num*p.num, den*p.den); }
	Rational operator/(const Rational &p) const { return Rational(num*p.den, den*p.num); }
	Rational operator+(const Polynomial &p) const { return Rational(num + p*den, den); }
	Rational operator-(const Polynomial &p) const { return Rational(num - p*den, den); }
	Rational operator*(const Polynomial &p) const { return Rational(num*p, den); }
	Rational operator/(const Polynomial &p) const { return Rational(num, den*p); }
	Rational operator*(float a) const { return Rational(num*a, den); }
	Rational operator/(float f) const { return Rational(num, den*f); }
	Rational operator+(float p) const { return Rational(num + p*den, den); }
	Rational operator-(float p) const { return Rational(num - p*den, den); }
	Rational operator-() const { return Rational(-num, den); }

	friend Rational operator+(float f, const Rational &p) { return Rational(f) + p; }
	friend Rational operator-(float f, const Rational &p) { return Rational(f) - p; }
	friend Rational operator*(float f, const Rational &p) { return Rational(f) * p; }
	friend Rational operator+(const Polynomial &f, const Rational &p) { return Rational(f) + p; }
	friend Rational operator-(const Polynomial &f, const Rational &p) { return Rational(f) - p; }
	friend Rational operator*(const Polynomial &f, const Rational &p) { return Rational(f) * p; }

	template<DivisionRing T> T operator()(const T& x) const { return num(x) / den(x); }
	Rational operator&(const Rational &p) const { return (*this)(p); }

	bool operator==(const Rational &p) const { return num == p.num && den == p.den || num == 0 && p.num == 0; }
	bool operator==(float c) const { return constant() && num == c; }
	bool operator==(const Polynomial &p) const { return constant() && num == p; }
};


enum FunctorialConstrctorStrategy {
	EXACT,
	TERMINALLY_NUMERICAL,
	MIXED,
	COMPOSED_NUMERICAL,
};


enum SymbolicPrimitiveTransform {
	ADD, MULT, COMPOSE, DERIVATIVE, ELEMENTARY, CONSTANT, LINEAR
};


class ElementarySymbolicFunction : public RealFunction {
public:
	std::string elementary_codename;
	ElementarySymbolicFunction(const RealFunction &f, std::string elementary_codename) :
		RealFunction(f), elementary_codename(std::move(elementary_codename)) {}
};



class SymbolicTransformGraphNode {
	float node_parameter;
	std::unique_ptr<SymbolicTransformGraphNode> _left = nullptr;
	std::unique_ptr<SymbolicTransformGraphNode> _right = nullptr;
public:

	std::string elementary_terminal_codename = "";
	SymbolicPrimitiveTransform _type;
	bool binary;
	bool terminal = false;

	SymbolicTransformGraphNode(SymbolicPrimitiveTransform type, SymbolicTransformGraphNode left, SymbolicTransformGraphNode right);
	SymbolicTransformGraphNode(SymbolicPrimitiveTransform type, SymbolicTransformGraphNode left, float node_parameter=0);
};
