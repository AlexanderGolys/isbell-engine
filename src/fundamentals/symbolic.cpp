
#include "symbolic.hpp"


using std::vector, std::string;



Polynomial::Polynomial(float term, int monomial_pow){
	if (nearlyEqual(term)) {
		coefs = vector<float>();
		return;
	}
	coefs = vector<float>(monomial_pow+1, 0);
	coefs[monomial_pow] = term;
}

void Polynomial::reduce_leading_zeros() {
	while (!coefs.empty() && nearlyEqual(coefs[coefs.size()-1])) coefs.pop_back();
}

std::set<int> Polynomial::term_degrees() const {
	auto res = std::set<int>();
	for (int i = 0; i < coefs.size(); i++)
		if (!nearlyEqual(coefs[i])) res.insert(i);
	return res;
}

void Polynomial::try_to_factorize() {
	if (degree() == 1) roots.insert(Complex(-coefs[0] / coefs[1]));
	if (degree() == 2) {
		float a = coefs[2], b = coefs[1], c = coefs[0];
		float d = b * b - 4 * a * c;
		roots.insert((-b + sqrt(Complex(d))) / (2 * a));
		roots.insert((-b - sqrt(Complex(d))) / (2 * a));
	}
}

Polynomial::Polynomial(const vector<float> &c) {
	coefs = c;
	reduce_leading_zeros();
}


Polynomial::Polynomial(const Polynomial &p){
	coefs = p.coefs;
	reduce_leading_zeros();
}

Polynomial::Polynomial(Polynomial &&other) noexcept: coefs(std::move(other.coefs)) {}

Polynomial & Polynomial::operator=(const Polynomial &other) {
	if (this == &other)
		return *this;
	coefs = other.coefs;
	roots = other.roots;
	return *this;
}

Polynomial & Polynomial::operator=(Polynomial &&other) noexcept {
	if (this == &other)
		return *this;
	coefs = std::move(other.coefs);
	roots = std::move(other.roots);
	return *this;
}

int Polynomial::degree() const { return coefs.size() - 1; }

bool Polynomial::monomial() const { return coefs.size() == 2; }

bool Polynomial::constant() const {
	return coefs.size() == 1 || coefs.empty();
}

float Polynomial::operator[](int i) const {
	if (i >= coefs.size()) return 0;
	return coefs[i];
}

void Polynomial::set(int i, float val) { coefs[i] = val; }

Polynomial Polynomial::operator+(const Polynomial &p) const {
	vector cc = vector<float>(std::max(degree(), p.degree()) + 1, 0);
	for (int i = 0; i <= degree(); i++)
		cc[i] += coefs[i];
	for (int i = 0; i <= p.degree(); i++)
		cc[i] += p[i];
	return Polynomial(cc);
}

Polynomial Polynomial::operator-(const Polynomial &p) const {
	return *this + -p;
}

Polynomial Polynomial::operator*(const Polynomial &p) const {
	vector cc = vector<float>(degree() + p.degree() + 1, 0);
	for (int i = 0; i <= degree(); i++)
		for (int j = 0; j <= p.degree(); j++)
			cc[i + j] += coefs[i] * p[j];

	auto res = Polynomial(cc);
	res.roots.insert(roots.begin(), roots.end());
	res.roots.insert(p.roots.begin(), p.roots.end());
	return res;
}

Polynomial Polynomial::operator*(float a) const {return *this * Polynomial(a);}

Polynomial Polynomial::operator/(float f) const { return *this * (1.f / f); }

Polynomial Polynomial::operator+(float p) const {
	if (coefs.empty()) return Polynomial(p);
	if (p == 0) return *this;
	auto c = coefs; c[0] += p; return Polynomial(c);
}

Polynomial Polynomial::operator-(float p) const { return *this + -p; }

Polynomial Polynomial::operator-() const { return *this * Polynomial(-1); }

SymbolicTransformGraphNode::SymbolicTransformGraphNode(SymbolicPrimitiveTransform type, SymbolicTransformGraphNode left, SymbolicTransformGraphNode right)
:	_type(type), node_parameter(0), binary(true){
	_left = std::make_unique<SymbolicTransformGraphNode>(std::move(left));
	_right = std::make_unique<SymbolicTransformGraphNode>(std::move(right));
}

SymbolicTransformGraphNode::SymbolicTransformGraphNode(SymbolicPrimitiveTransform type, SymbolicTransformGraphNode left, float node_parameter)
: node_parameter(node_parameter), _type(type), binary(false){
	_left = std::make_unique<SymbolicTransformGraphNode>(std::move(left));
}

Polynomial operator+(float f, const Polynomial &p) { return p + f; }

Polynomial operator-(float f, const Polynomial &p) { return -p + f; }

Polynomial operator*(float f, const Polynomial &p) { return p * f; }
