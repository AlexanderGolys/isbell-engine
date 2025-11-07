#pragma once
#include <ranges>

#include "exceptions.hpp"


template<typename R>
concept SymbolicRing = requires(R a, R b) {
	{ a + b } -> std::same_as<R>;
	{ a * b } -> std::same_as<R>;
	{ a - b } -> std::same_as<R>;
	{ R::one() } -> std::same_as<R>;
	{ R::zero()} -> std::same_as<R>;
	{ -a } -> std::same_as<R>;
};

template<typename R>
concept SymbolicField = SymbolicRing<R> && requires(R a, R b) {
	{ a / b } -> std::same_as<R>;
	{ a == b } -> std::convertible_to<bool>;
};

template<SymbolicRing R>
class Monomial {
	R coefficient;
	vector<pair<string, int>> varExpPairs;

public:
	explicit Monomial(const R& coefficient, const vector<pair<string, int>>& varExpPairs)
		: coefficient(coefficient), varExpPairs(varExpPairs) {
		std::sort(this->varExpPairs.begin(), this->varExpPairs.end(),
				  [](const pair<string, int>& a, const pair<string, int>& b) {
					  return a.first < b.first;
				  });
	}

	const R& getCoefficient() const {
		return coefficient;
	}

	const vector<pair<string, int>>& getVarExpPairs() const {
		return varExpPairs;
	}

	int degree() const {
		int deg = 0;
		for (const auto& exp : varExpPairs | std::views::values)
			deg += exp;
		return deg;
	}

	Monomial operator*(const Monomial& other) const {
		R newCoefficient = this->coefficient * other.coefficient;
		vector<pair<string, int>> newVarExpPairs = this->varExpPairs;
		THROW_IF(newVarExpPairs.size() != other.varExpPairs.size(), ValueError, "Multiplication of monomials with different variables is not supported in this implementation.");
		for (int i = 0; i < other.varExpPairs.size(); i++) {
			const auto& [x, a] = other.varExpPairs[i];
			const auto& [y, b] = newVarExpPairs[i];
			THROW_IF(x != y, ValueError, "Multiplication of monomials with different variables is not supported in this implementation.");
			newVarExpPairs[i].second = a + b;
		}
		return Monomial(newCoefficient, newVarExpPairs);
	}
	Monomial operator*(R other) const {
		return Monomial(this->coefficient * other, this->varExpPairs);
	}
	bool operator==(const Monomial& other) const {
		return this->coefficient == other.coefficient && this->varExpPairs == other.varExpPairs;
	}
	bool operator!=(const Monomial& other) const {
		return !(*this == other);
	}
	bool operator<(const Monomial& other) const {
		if (this->varExpPairs == other.varExpPairs)
			return false;
		if (degree() != other.degree())
			return degree() < other.degree();
		THROW_IF(this->varExpPairs.size() != other.varExpPairs.size(), ValueError, "Comparison of monomials with different variables is not supported in this implementation.");
		for (int i = 0; i < this->varExpPairs.size(), other.varExpPairs.size(); i++) {
			const auto& [x, a] = this->varExpPairs[i];
			const auto& [y, b] = other.varExpPairs[i];
			THROW_IF(x != y, ValueError, "Comparison of monomials with different variables is not supported in this implementation.");
			if (a < b)
				return true;
			if (a > b)
				return false;
		}
		return false;
	}
};

template<SymbolicRing R>
class Polynomial {
	vector<Monomial<R>> terms;
public:
	explicit Polynomial(const vector<Monomial<R>>& terms) : terms(terms) {
		std::sort(this->terms.begin(), this->terms.end(),
				  [](const Monomial<R>& a, const Monomial<R>& b) {
					  return a < b;
				  });
	}

	Polynomial operator+(const Polynomial& other) const {
		vector<Monomial<R>> newTerms = this->terms;
		newTerms.insert(newTerms.end(), other.terms.begin(), other.terms.end());
		return Polynomial(newTerms);
	}


};