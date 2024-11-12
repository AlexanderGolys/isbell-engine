# pragma once


#include <string>
#include <vector>
#include <iostream>
#include <cmath>
#include <array>
#include <algorithm>
#include "macros.hpp"
#include "metaUtils.hpp"
#include "glmUX.hpp"




template <typename T>
concept AbelianSemigroup = requires (T a, T b) {
	{ a + b } -> std::same_as<T>;
};

template <typename T>
concept MultiplicativeSemigroup = requires (T a, T b) {
	{ a*b } -> std::same_as<T>;
};

template <typename T>
concept AbelianMonoid = AbelianSemigroup<T> && requires {
	T(0);
};

template <typename G>
concept Monoid = MultiplicativeSemigroup<G> &&  requires {
	G(1); 
};

template <typename G>
concept GroupConcept = Monoid<G> && requires (G g) 
{
	{ ~g } -> std::same_as<G>;
};

template <typename G>
concept AbelianGroupConcept = AbelianMonoid<G> && requires (G g)
{
	{ -g } -> std::same_as<G>;
};

template <typename T>
concept Rng = requires (T a, T b) {
	{ a + b } -> std::same_as<T>;
	-a;
	{ a* b } -> std::same_as<T>;
	T(0);
};

template <typename T>
concept RingConcept = Rng<T> && Monoid<T>;


template <typename T>
concept DivisionRing = RingConcept<T> && GroupConcept<T>;



template <typename V, typename K> 
concept VectorSpaceConcept = DivisionRing<K> && requires (V a, V b, K k) {
	{a + b} -> std::same_as<V>;
	-a;
	{a*k} -> std::same_as<V>;
};

template <typename T> 
concept Normed = requires (T a, float c) {
	norm2(a);
	a / c;
};

template <typename V, typename K=float>
concept EuclideanSpace = VectorSpaceConcept<V, K> && requires (V a, V b, K c) {
	{c} -> std::convertible_to<float>;
	{dot(a, b)} -> std::convertible_to<float>;
};


template<RingConcept T, int n, int m=n> 
class Matrix {
	std::array<std::array<T, n>, m> coefs;
public: 

	Matrix()
	{
		this->coefs = std::array<std::array<T, n>, m>();
		for (int i = 0; i < m; i++)
			for (int j = 0; j < n; j++)
				this->coefs[i][j] = 0;
	}
	explicit Matrix(std::array<std::array<T, n>, m> c)
	{
		this->coefs = c;
	}
	explicit Matrix(std::vector<std::vector<T>> c)
	{
		this->coefs = std::array<std::array<T, n>, m>();
		for (int i = 0; i < m; i++)
			for (int j = 0; j < n; j++)
				this->coefs[i][j] = c[i][j];
	}
	Matrix operator*(const T &f) const;
	Matrix operator/(const T &f) const requires DivisionRing<T>;
	Matrix operator+(const Matrix &M) const;
	Matrix operator-(const Matrix &M) const;
	Matrix<T, m, n> transpose() const;
	Matrix operator-() const;
	operator std::string() const;
	static Matrix zero() { return Matrix(); }
	T at(int i, int j) const { return this->coefs[i][j]; }
};


template <RingConcept T, int n>
class Matrix<T, n> {
	std::array<std::array<T, n>, n> coefs;
public:
	Matrix(T diag) : Matrix() {
        for (int i = 0; i < n; i++)
            this->coefs[i][i] = diag;
    }
	static Matrix identity() {return Matrix(T(1));}
	static Matrix Id() { return identity(); }
	T minor(int i, int j) const;
	T det() const;
	Matrix pow(int p) const;
	T trace() const;
	Matrix inv() const requires DivisionRing<T>;
	Matrix operator~() const requires DivisionRing<T> { return inv(); }
	Matrix GramSchmidtProcess() const requires EuclideanSpace<T>;
};


template <RingConcept T>
class Matrix<T, 2> {
	std::array<std::array<T, 2>, 2> coefs;
public:
	T a, b, c, d;
	Matrix(T a, T b, T c, T d) {
        this->a = a;
        this->b = b;
        this->c = c;
        this->d = d;
		this->coefs = std::array<std::array<T, 2>, 2>({ { {a, b}, {c, d} } });
    }
	explicit Matrix(std::array<std::array<T, 2>, 2> c) {
        this->a = c[0][0];
        this->b = c[0][1];
        this->c = c[1][0];
        this->d = c[1][1];
		this->coefs = c;
    }
	explicit Matrix(std::vector<std::vector<T>> c)
	{
		this->coefs = std::array<std::array<T, 2>, 2>();
		for (int i = 0; i < 2; i++)
			for (int j = 0; j < 2; j++)
				this->coefs[i][j] = c[i][j];
		this->a = c[0][0];
		this->b = c[0][1];
		this->c = c[1][0];
		this->d = c[1][1];
	}
	Matrix(T diag) : Matrix(diag, T(0), T(0), diag) {}
	T mobius_derivative(T z) const requires DivisionRing<T>;
	T det() const { return a * d - b * c; };
	T trace() const { return a + d; };
	Matrix inv() const requires DivisionRing<T> { return Matrix(d, -b, -c, a) / det(); }
	Matrix operator~() const requires DivisionRing<T> { return inv(); }
	Matrix pow(int p) const {
		if (p < 0)
            return ~(*this).pow(-p);
		if (p == 0)
            return Matrix(T(1));
		if (p == 1)
            return *this;
		if (p % 2 == 0)
            return this->pow(p / 2) * this->pow(p / 2);
		return (*this) * this->pow(p / 2) * this->pow(p / 2);
	}
	Matrix operator*(const T &f) const { return Matrix(a * f, b * f, c * f, d * f); }
	Matrix operator/(const T &f) const requires DivisionRing<T> { return *this * (T(1) / f); }
	Matrix operator+(const Matrix &M) const { return Matrix(a + M.a, b + M.b, c + M.c, d + M.d); }
	Matrix operator-(const Matrix &M) const { return *this + M * -1; }
	Matrix operator*(const Matrix &M) const { return Matrix(a * M.a + b * M.c, a * M.b + b * M.d, c * M.a + d * M.c, c * M.b + d * M.d); }
	Matrix transpose() const { return Matrix(a, c, b, d); }
	Matrix operator-() const { return *this * -1; }
	static Matrix zero() { return Matrix(); }
	T at(int i, int j) const { return this->coefs[i][j]; }
    T mobius(T z) const requires DivisionRing<T> { return (a * z + b) / (c * z + d); }
};




template <RingConcept T, int n, int m, int k>
Matrix<T, n, k> operator*(Matrix<T, n, m> A, Matrix<T, m, k> B);




float norm2(glm::vec2 v);
float norm2(glm::vec3 v);
float norm2(glm::vec4 v);
inline float det(glm::mat3 m) {return glm::determinant(m);}
inline float det(glm::mat2 m) {return glm::determinant(m);}
inline float det(glm::mat4 m) {return glm::determinant(m);}



template <typename T>
float norm(T v) { return sqrt(norm2(v)); };

template <typename T>
T normalise(T v) { return v / float(sqrt(norm2(v))); };


glm::vec2 intersectLines(glm::vec2 p1, glm::vec2 p2, glm::vec2 q1, glm::vec2 q2);
glm::vec2 orthogonalComplement(glm::vec2 v);
std::pair<glm::vec3, glm::vec3> orthogonalComplementBasis(glm::vec3 v);
glm::mat3 GramSchmidtProcess(glm::mat3 m);


class Complex {
public: 
	glm::vec2 z;
	float x;
	float y;
	Complex();
	Complex(glm::vec2 z); // NOLINT(*-explicit-constructor)
	Complex(float x, float y);
	Complex(float x); // NOLINT(*-explicit-constructor)
	auto operator+(Complex c) const -> Complex;
	auto operator-(Complex c) const -> Complex;
	auto operator*(Complex c) const -> Complex;
	auto operator/(Complex c) const -> Complex;
	auto operator~() const -> Complex;
	auto operator-() const -> Complex;
	void operator+=(Complex c);
	void operator-=(Complex c);
	void operator*=(Complex c);
	void operator/=(Complex c);

	Complex inv() const { return ~(*this); }
	auto operator==(Complex c) const -> bool;
	auto operator==(float f) const -> bool;
	static auto one() -> Complex;
	static auto zero() -> Complex;
	auto expForm() const -> glm::vec2;
	auto arg() const -> float;
	auto conj() const -> Complex;
	auto pow(float exponent) const -> Complex;
	auto re() const -> float;
	auto im() const -> float;

	auto real() const -> float { return re(); }
	auto imag() const -> float { return im(); }
	auto square() const -> Complex;
	auto sqrt() const -> Complex;
	friend auto operator<<(std::ostream &_stream, Complex const &z) -> std::ostream &;
	operator glm::vec2() const; // NOLINT(*-explicit-constructor)
	explicit operator std::string() const;


	bool nearlyEqual(Complex c) const;
};

auto norm2(Complex c) -> float;
auto abs(Complex c) -> float;

inline Complex operator*(float f, Complex c) { return c * f; }
inline Complex operator/(float f, Complex c) { return Complex(f) / c; }
inline Complex operator+(float f, Complex c) { return c + f; }
inline Complex operator-(float f, Complex c) { return Complex(f) - c; }


template<int n, int m>
bool nearlyEqual(Matrix<Complex, n, m> a, Matrix<Complex, n, m> b) {
    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
            if (!a.at(i, j).nearlyEqual(b.at(i, j)))
                return false;
    return true;
}

inline bool nearlyEqual(glm::mat3 a, glm::mat3 b) {
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            if (abs(a[i][j] - b[i][j]) > 1e-6)
                return false;
    return true;
}

inline bool nearlyEqual(float a, float b) {
    return abs(a - b) < 1e-6;
}


Complex fromExpForm(glm::vec2 e);
Complex rootOfUnity(int n, int k);
std::vector<Complex> allRootsOf1(int n);
Complex exp(Complex c);
Complex log(Complex c);
Complex sin(Complex c);
Complex cos(Complex c);
Complex tan(Complex c);
Complex sinh(Complex c);
Complex cosh(Complex c);
Complex tanh(Complex c);
Complex asin(Complex c);
Complex acos(Complex c);
Complex atan(Complex c);
Complex asinh(Complex c);
Complex acosh(Complex c);
Complex atanh(Complex c);




class CP1 {
public:
	Complex z;
	bool inf;
	CP1();
	CP1(Complex z);
	CP1(Complex z, bool inf);
	CP1(Complex z1, Complex z2);
	CP1(float x);
	CP1(float x, bool inf);

	CP1 operator+(CP1 c) const;
	CP1 operator-(CP1 c) const;
	CP1 operator*(CP1 c) const;
	CP1 operator/(CP1 c) const;
	CP1 operator+(Complex c) const;
	CP1 operator-(Complex c) const;
	CP1 operator*(Complex c);
	CP1 operator/(Complex c);
	CP1 operator+(float f);
	CP1 operator-(float f);
	CP1 operator*(float f);
	CP1 operator/(float f);
	CP1 operator-();
	CP1 inv();
	CP1 square();
	CP1 sqrt();
	CP1 pow(float f);
	CP1 conj();
	operator Complex();
};

CP1 oneCP1();
CP1 zeroCP1();
CP1 infCP1();
CP1 iCP1();

Complex intersectLines(Complex p1, Complex p2, Complex q1, Complex q2);


template<typename V>
V lerp(V a, V b, float t)
{
	return b * t + a * (1.f - t);
}

template<typename V>
V mean(std::vector<V> points)
{
	V m = V(0);
	for (V p : points) {
		m = m + p;
	}
	return m / (1.f*points.size());
}

glm::mat3 scaleMatrix3(glm::vec3 s);
glm::mat3 scaleMatrix3(float s);
glm::mat3 changeOfBasis(glm::vec3 target1, glm::vec3 target2, glm::vec3 target3);
glm::mat3 changeOfBasis(glm::vec3 source1, glm::vec3 source2, glm::vec3 source3, glm::vec3 target1, glm::vec3 target2, glm::vec3 target3);
glm::mat3 rotationMatrix3(float angle);
glm::mat3 rotationMatrix3(glm::vec3 axis, float angle);

template<typename T>
int binSearch(std::vector<T> v, T x)
{
    int l = 0;
    int r = v.size() - 1;
    while (l <= r) {
        int m = l + (r - l) / 2;
        if (v[m] == x)
            return m;
        if (v[m] < x)
            l = m + 1;
        else
            r = m - 1;
    }
    return 0;
}

float frac(float x);

float pseudorandomizer(float x, float seed=0.f);
























template <RingConcept T, int n, int m>
Matrix<T, n, m> Matrix<T, n, m>::operator*(const T& f) const
{
	Matrix result = Matrix();
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			result.coefs[i][j] = this->coefs[i][j] * f;
		}
	}
	return result;
};


template <RingConcept T, int n, int m>
Matrix<T, n, m> Matrix<T, n, m>::operator+(const Matrix& M) const
{
	Matrix result = Matrix();
	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++)
			result.coefs[i][j] = this->coefs[i][j] + M.coefs[i][j];
	return result;
}

template <RingConcept T, int n, int m>
Matrix<T, n, m> Matrix<T, n, m>::operator-(const Matrix<T, n, m>& M) const
{
	Matrix<T, n, m> result = Matrix<T, n, m>();
	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++)
			result.coefs[i][j] = this->coefs[i][j] - M.coefs[i][j];
	return result;
}

template <RingConcept T, int n, int m>
Matrix<T, m, n> Matrix<T, n, m>::transpose() const
{
	Matrix<T, m, n> result;
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			result.coefs[i][j] = this->coefs[j][i];
		}
	}
	return result;
}

template <RingConcept T, int n, int m>
Matrix<T, n, m> Matrix<T, n, m>::operator-() const
{
	return *this * -1;
}



template<RingConcept T, int n>
T Matrix<T, n>::minor(int i, int j) const {
	if (n <= 1) {
		throw std::invalid_argument("Matrix must be at least 2x2");
	}
	auto cfs = std::vector<std::vector<T>>();
	for (int k = 0; k < n; k++) {
		if (k == i) continue;
		auto row = std::vector<T>();
		for (int l = 0; l < n; l++) {
			if (l == j) continue;
			row.push_back(this->coefs[k][l]);
		}
		cfs.push_back(row);
	}
	return Matrix<T, std::max(1,n - 1)>(cfs).det();
}


template<RingConcept T, int n, int m, int k>
Matrix<T, n, k> operator*(Matrix<T, n, m> A, Matrix<T, m, k> B)
{
	Matrix<T, n, k> result;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < k; j++) {
			T sum = 0;
			for (int l = 0; l < m; l++) {
				sum += A.coefs[i][l] * B.coefs[l][j];
			}
			result.coefs[i][j] = sum;
		}
	}
	if (n == 2 && k == 2) {
		result.a = result.coefs[0][0];
		result.b = result.coefs[0][1];
		result.c = result.coefs[1][0];
		result.d = result.coefs[1][1];
	}
	return result;
}

template <RingConcept T, int n, int m>
Matrix<T, n, m> Matrix<T, n, m>::operator/(const T& f) const requires DivisionRing<T> {
	Matrix<T, n, m> result = Matrix<T, n, m>();
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			result.coefs[i][j] = this->coefs[i][j] / f;
		}
	}
	return result;
}




template <RingConcept T, int n>
T Matrix<T, n>::det() const
{
	if (n == 1)
		return this->coefs[0][0];

	if (n <= 0)
		return 0;

	if (n == 2)
		return this->coefs[0][0] * this->coefs[1][1] - this->coefs[0][1] * this->coefs[1][0];

	T sum = 0;
	for (int i = 0; i < n; i++) {
		int sign = (i % 2 == 0) ? 1 : -1;
		sum = sum +  this->coefs[0][i] * this->minor(0, i) * sign;
	}
	return sum;
}



template<RingConcept T, int n>
Matrix<T, n> Matrix<T, n>::pow(int p) const
{
	if (p<0)
		return ~(*this).pow(-p);
	if (p == 0)
		return Matrix::identity();
	if (p == 1)
		return *this;
	if (p % 2 == 0)
		return this->pow(p / 2) * this->pow(p / 2);
	return (*this) * this->pow(p / 2) * this->pow(p / 2);
}

template <RingConcept T, int n>
T Matrix<T, n>::trace() const
{
	T tr = 0;
	for (int i = 0; i < n; i++) {
		tr = tr +  this->coefs[i][i];
	}
	return tr;
}

template <RingConcept T, int n, int m>
Matrix<T, n, m>::operator std::string() const
{
	std::string s = "";
	for (int i = 0; i < m; i++) {
		s = s + "|";
		for (int j = 0; j < n; j++) {
			s = s + std::string(this->coefs[i][j]);
		}
		s = s + "|\n";
	}
	return s;
}

template <RingConcept T, int n>
Matrix<T, n> Matrix<T, n>::inv() const requires DivisionRing<T>
{
	if (n == 1)
		return Matrix<T, n>(T(1) / this->coefs[0][0]);
	Matrix<T, n> I;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			int sign = ((i + j) % 2 == 0) ? 1 : -1;
			I.coefs[j][i] = this->minor(i, j) * sign;
		}
	}
	return I / this->det();
}

template<RingConcept T, int n>
Matrix<T, n> Matrix<T, n>::GramSchmidtProcess() const requires EuclideanSpace<T> {
	Matrix<T, n> result;
    for (int i = 0; i < n; i++) {
        result.coefs[i][i] = this->coefs[i][i];
        for (int j = 0; j < i; j++) {
            T proj = dot(result.coefs[j], this->coefs[i]) / norm2(result.coefs[j]);
            result.coefs[i] = result.coefs[i] - proj * result.coefs[j];
        }
    }
    return result;
}


template<RingConcept T>
T Matrix<T, 2>::mobius_derivative(T z) const requires DivisionRing<T> {
	return (a * d - b * c) / ((c * z + d) * (c * z + d));
}


const Complex ONE = Complex(1, 0);
const Complex ZERO = Complex(0, 0);
const Complex I = Complex(0, 1);

template<typename vec>
std::vector<float> vecToVecHeHe(vec v) {
	std::vector<float> res;
	res.reserve(v.length());
	for (int i = 0; i < v.length(); i++)
		res.push_back(v[i]);
	return res;
};

template <typename vec>
vec barycenter(vec a, vec b, vec c) {
    return (a + b + c) / 3.f;
}





