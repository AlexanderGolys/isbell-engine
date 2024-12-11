# pragma once

#include <cmath>
#include <vector>
#include <iostream>
#include <cmath>
#include <array>
#include <algorithm>
#include <format>
#include <memory>


#include "metaUtils.hpp"



template <typename T>
concept AbelianSemigroup = requires (T a, T b) {
	{ a + b } -> std::convertible_to<T>;
};

template <typename T>
concept Semigroup = requires (T a, T b) {
	{ a*b } -> std::convertible_to<T>;
};

template <typename T>
concept TotalOrder =  requires (T a, T b) {
	{ a < b } -> std::convertible_to<bool>;
};

template <typename T>
concept TotallyOrderedAbelianSemigroup = AbelianSemigroup<T> && TotalOrder<T>;

template <typename T>
concept AbelianMonoid = AbelianSemigroup<T> && requires {
	{T(0)} -> std::convertible_to<T>;
};

template <typename G>
concept Monoid = Semigroup<G> &&  requires {
	{G(1)} -> std::convertible_to<G>;
};

template <typename G>
concept GroupConcept = Monoid<G> && (requires (G g) {
	{ ~g } -> std::convertible_to<G>; }
	|| requires (G g) {
	{ G(1)/g } -> std::convertible_to<G>; }
	|| requires (G g) {
	{ 1.f/g } -> std::convertible_to<G>; });

template <typename G>
concept AbelianGroupConcept = AbelianMonoid<G> && requires (G g) {
	{ -g } -> std::convertible_to<G>; };

template <typename T>
concept Rng = AbelianGroupConcept<T> && Semigroup<T>;


template <typename T>
concept RingConcept = Rng<T> && Monoid<T>;


template <typename T>
concept DivisionRing = RingConcept<T> && GroupConcept<T>;

template <typename A, typename R>
concept ModuleConcept = Rng<R> && AbelianGroupConcept<A> && requires (A a, R r) {
    {a*r} -> std::convertible_to<A>; };

template <typename A, typename R>
concept Algebra = ModuleConcept<A, R> && Rng<A>;


template <typename V, typename K> 
concept VectorSpaceConcept = ModuleConcept<V, K> && DivisionRing<K>;

template <typename T> 
concept Normed = requires (T a, float c) {
	norm(a);
	a/c; };

template <typename V>
concept EuclideanSpaceConcept = VectorSpaceConcept<V, float> && requires (V a, V b) {
	{dot(a, b)} -> std::convertible_to<float>; };


template<typename vec>
std::vector<float> vecToVecHeHe(vec v) {
	std::vector<float> res;
	for (int i = 0; i < v.length(); i++)
		res.push_back(v[i]);
	return res;
};




template<typename domain, typename codomain=float>
class Morphism {
protected:
    std::function<codomain(domain)> _f;
public:
    Morphism(std::function<codomain(domain)> f) : _f(f) {} // NOLINT(*-explicit-constructor)
    codomain operator()(domain x) const { return _f(x); }
    Morphism operator+(const Morphism &g) const requires AbelianSemigroup<codomain> { return Morphism([f_=_f, g_=g._f](domain x) { return f_(x) + g_(x); }); }
    // friend Morphism operator+(codomain a, const Morphism &f) requires AbelianSemigroup<codomain> { return Morphism([f_ = _f, a](domain x) { return f_(x) + a; }); }

     Morphism operator-(const Morphism &g) const requires AbelianGroupConcept<codomain> { return Morphism([this, g](domain x) { return (*this)(x) - g(x); }); }
     Morphism operator*(const Morphism &g) const requires Semigroup<codomain> { return Morphism([this, g](domain x) { return (*this)(x) * g(x); }); }
     Morphism operator/(const Morphism &g) const requires DivisionRing<codomain> { return Morphism([this, g](domain x) { return (*this)(x) / g(x); }); }
     Morphism operator*(codomain a) const requires Semigroup<codomain> { return Morphism([this, a](domain x) { return (*this)(x) * a; }); }
    // friend Morphism operator*(codomain a, const Morphism &f) requires GroupConcept<codomain> { return Morphism([f, a](domain x) { return a*f(x); }); }

    template<DivisionRing K>  Morphism operator/(K a) const requires VectorSpaceConcept<codomain, K>      { return Morphism([this, a](domain x) { return (*this)(x) / a; }); }
    template<Rng R>           Morphism operator*(R a) const requires ModuleConcept<codomain, R>           { return Morphism([this, a](domain x) { return (*this)(x) * a; }); }
};

template<typename domain, typename codomain>
class Isomorphism : public Morphism<domain, codomain> {
    Morphism<codomain, domain> _inverse;
public:
    Isomorphism(std::function<codomain(domain)> f, std::function<domain(codomain)> g) : Morphism<domain, codomain>(f), _inverse(g) {}
    Isomorphism inverseMorphism() const { return Isomorphism(_inverse, this->_f); }
    Isomorphism operator~() const { return inverseMorphism(); }
    domain inv(codomain y) const { return _inverse(y); }
};

template<typename domain>
class Endomorphism : public Morphism<domain, domain> {
public:
    using Morphism<domain, domain>::Morphism;
    static Endomorphism id() { return Endomorphism([](domain x) { return x; }); }
    Endomorphism pow(int p) const {
        if (p < 0) return ~(*this).pow(-p);
        if (p == 0) return Endomorphism::id();
        if (p == 1) return *this;
        if (p % 2 == 0) return this->pow(p / 2) * this->pow(p / 2);
        return (*this) * this->pow(p / 2) * this->pow(p / 2); }
    Endomorphism operator^(int p) const { return pow(p); }
};



// composition, used with operator f&g
template <typename X, typename Y, typename Z>
Morphism<X, Z> compose(const Morphism<Y, Z> &f, const Morphism<X, Y> &g) {
    return Morphism<X, Z>([f_=f, g_=g](X x) { return f_(g_(x)); });
}
template <typename X, typename Y, typename Z>
Morphism<X, Z> operator&(const Morphism<Y, Z> &f, const Morphism<X, Y> &g) {
    return Morphism<X, Z>([f_=f, g_=g](X x) { return f_(g_(x)); });
}









// -----------------------------  MATRICES  ----------------------------------


template<RingConcept T, int n, int m=n> 
class Matrix {
	std::array<std::array<T, n>, m> coefs;

public:
	Matrix();
	explicit Matrix(std::array<std::array<T, n>, m> c);
	explicit Matrix(std::vector<std::vector<T>> c);

	Matrix operator*(const T &f) const;
	Matrix operator/(const T &f) const requires DivisionRing<T>;
	Matrix operator+(const Matrix &M) const;
	Matrix operator-(const Matrix &M) const;
	Matrix<T, m, n> transpose() const;
	Matrix operator-() const;
	explicit operator std::string() const;
	static Matrix zero() { return Matrix(); }
	T at(int i, int j) const { return this->coefs[i][j]; }
};


template <RingConcept T, int n>
class Matrix<T, n> {
	std::array<std::array<T, n>, n> coefs;

public:
	explicit Matrix(T diag);

	T minor(int i, int j) const;
	T det() const;
	Matrix pow(int p) const;
	T trace() const;
	Matrix inv() const requires DivisionRing<T>;
	Matrix operator~() const requires DivisionRing<T> { return inv(); }
	Matrix GramSchmidtProcess() const requires EuclideanSpaceConcept<T>;

	static Matrix identity() {return Matrix(T(1));}
	static Matrix Id() { return identity(); }
};


template <>
class Matrix<float, 3> {
public:
	explicit Matrix(mat3 m);
	Matrix(vec3 v1, vec3 v2, vec3 v3) : Matrix(mat3(v1, v2, v3)) {}
};

template <>
class Matrix<float, 4> {
public:
	explicit Matrix(mat4 m);
	Matrix(vec4 v1, vec4 v2, vec4 v3, vec4 v4) : Matrix(mat4(v1, v2, v3, v4)) {}
};

template <RingConcept T>
class Matrix<T, 2> {
	std::array<std::array<T, 2>, 2> coefs;

public:
	T a, b, c, d;
	Matrix(T a, T b, T c, T d);
	explicit Matrix(std::array<std::array<T, 2>, 2> c);
	explicit Matrix(std::vector<std::vector<T>> c);
	explicit Matrix(T diag) : Matrix(diag, T(0), T(0), diag) {}

	T mobius_derivative(T z) const requires DivisionRing<T>;
	T det() const { return a * d - b * c; };
	T trace() const { return a + d; };
	Matrix inv() const requires DivisionRing<T> { return Matrix(d, -b, -c, a) / det(); }
	Matrix pow(int p) const;
	Matrix transpose() const { return Matrix(a, c, b, d); }
	T at(int i, int j) const { return this->coefs[i][j]; }
	T mobius(T z) const requires DivisionRing<T> { return (a * z + b) / (c * z + d); }

	Matrix operator~() const requires DivisionRing<T> { return inv(); }
	Matrix operator*(const T &f) const { return Matrix(a * f, b * f, c * f, d * f); }
	Matrix operator/(const T &f) const requires DivisionRing<T> { return *this * (T(1) / f); }
	Matrix operator+(const Matrix &M) const { return Matrix(a + M.a, b + M.b, c + M.c, d + M.d); }
	Matrix operator-(const Matrix &M) const { return *this + M * -1; }
	Matrix operator*(const Matrix &M) const { return Matrix(a * M.a + b * M.c, a * M.b + b * M.d, c * M.a + d * M.c, c * M.b + d * M.d); }
	Matrix operator-() const { return *this * -1; }

	static Matrix zero() { return Matrix(T(0)); }
};



template <RingConcept T, int n, int m, int k>
Matrix<T, n, k> operator*(Matrix<T, n, m> A, Matrix<T, m, k> B);




class SparseMatrix {
  std::vector<std::vector<std::pair<int, float>>> data;
    int n, m;

public:
    SparseMatrix(int n, int m);

    void set(int i, int j, float val);
    float get(int i, int j);
	glm::ivec2 size() const { return glm::ivec2(n, m); }

    float operator()(int i, int j) { return get(i, j); }
    SparseMatrix operator* (float f);
    SparseMatrix operator+ (SparseMatrix M);
    SparseMatrix operator-(SparseMatrix M);
};




class BigVector {
	std::vector<float> data;
public:
	explicit BigVector(const std::vector<float> &data) { this->data = data; }
	explicit BigVector(vec2 data) { this->data = vecToVecHeHe(data); }
	explicit BigVector(vec3 data) { this->data = vecToVecHeHe(data); }
	explicit BigVector(vec4 data) { this->data = vecToVecHeHe(data); }
	explicit BigVector(const std::vector<std::vector<float>> &data);
	BigVector(int n, float val);
	explicit BigVector(int n) : BigVector(n, 0) {}

	BigVector(const BigVector &other) = default;
	BigVector(BigVector &&other) noexcept : data(std::move(other.data)) {}
	BigVector & operator=(const BigVector &other);
	BigVector & operator=(BigVector &&other) noexcept;


	float operator[](int i) const { return this->data[i]; }
	BigVector operator*(float f) const;
	BigVector operator/(float f) const { return *this * (1.f / f); }
	BigVector operator+(const BigVector& v) const;
	BigVector operator-(const BigVector& v) const;
	BigVector operator-() const { return *this * -1; }
	void operator+=(const BigVector &v) { for (int i = 0; i < this->data.size(); i++) this->data[i] += v[i]; }
	void operator-=(const BigVector &v) { for (int i = 0; i < this->data.size(); i++) this->data[i] -= v[i]; }
	void operator*=(float f) { for (float & i : this->data) i *= f; }
	void operator/=(float f) { *this *= 1.f / f; }

	void append (float f) { this->data.push_back(f); }
	void append (const BigVector &v) { for (float f : v.data) this->data.push_back(f); }
	void append (const vec69 &v) { for (float f : v) this->data.push_back(f); }

	std::vector<float> getVec() const { return this->data; }
	int size() { return this->data.size(); }

	friend float dot(BigVector a, BigVector b);
	friend BigVector concat(const BigVector& a, const BigVector& b) { BigVector res = a; res.append(b); return res; }
};


class BigMatrix {
  MATR$X data;
public:
  BigMatrix(int n, int m);
  explicit BigMatrix(const MATR$X &data);
  explicit BigMatrix(const vector<float> &data);
  explicit BigMatrix(const vector<vec2> &data);
  explicit BigMatrix(const vector<vec3> &data);
  explicit BigMatrix(const vector<vec4> &data);
  explicit BigMatrix(const vector<mat3> &data);

  explicit BigMatrix(MATR$X &&data);
  explicit BigMatrix(vec69 &&data) : BigMatrix({data}) {}
//  explicit BigMatrix(const mat3 &data) : BigMatrix(vecToVecHeHe(data)) {}

  void set(int i, int j, float val) { this->data[i][j] = val; }
  float get(int i, int j) const { return this->data[i][j]; }
  bool isSquare() const { return this->n() == this->m(); }
  float det();
  BigMatrix transpose() const;
  BigMatrix operator*(float f) const;
  BigMatrix operator+(const BigMatrix &M) const;
  BigMatrix operator-(const BigMatrix &M) const;
  BigMatrix operator*(const BigMatrix &M) const;
  BigMatrix operator*(const MATR$X &M) const;

  vec69 operator*(const vec69 &v) const { return (*this * vector<vec69>{v})[0]; }
  BigVector operator*(const BigVector &v) const { return BigVector(*this * v.getVec()); }

  BigMatrix operator-() const { return *this * (-1.f); }
  BigMatrix operator/(float x) const { return *this * (1.f / x); }
  BigMatrix inv();
  BigMatrix pow(int p);
  BigMatrix operator~() { return inv(); }
  BigMatrix GramSchmidtProcess();
  BigMatrix submatrix(int i, int j);
	BigMatrix diagonalComponent() const;
  BigMatrix invertedDiagonal() const;
	BigMatrix subtractedDiagonal() const { return *this - diagonalComponent(); }

  vec69 operator[] (int i) { return this->data[i]; }
	friend BigMatrix operator*(const MATR$X &M, const BigMatrix &B);

  glm::ivec2 size() const {return  glm::ivec2(n(), m());}
  int n() const { return this->data.size(); }
  int m() const { return this->data[0].size(); }

    explicit operator float () { if (n() != m() || n() != 1) throw std::format_error("wrong dimension of matrix (not 1x1)"); return this->data[0][0]; }
    explicit operator vec2 ();
    explicit operator vec3 ();
    explicit operator vec4 ();
};


template <Rng R=float, ModuleConcept<R> M=R>
class GenericTensor {
	std::vector<M> data;
	int dim_;

public:
	GenericTensor(int length, M fill, int dim=1) : dim_(dim) {data = std::vector<M>(); data.reserve(length); for (int i = 0; i < length; i++) data.push_back(fill);}
	explicit GenericTensor(vector<M> data, int dim=1) : data(data), dim_(dim) {}


	M operator[](int i) const { return this->data[i]; }

	template <ModuleConcept<R> E0=R> E0 at(int i, int j) const { return this->data[i][j]; }
	template <ModuleConcept<R> E0=R> E0 at(int i, int j, int k) const { return this->data[i][j][k]; }
	template <ModuleConcept<R> E0=R> E0 at(vector<int> indices) const;
	template <ModuleConcept<R> E0=R> E0 at(int i) const { return this->data[i]; }
	template <ModuleConcept<R> E0=R> void set(vector<int> ind, E0 val);
	template <ModuleConcept<R> E0=R> void set(int i, E0 val) { this->data[i] = val; }
	template <ModuleConcept<R> E0=R> void set(int i, int j, E0 val) { this->data[i][j] = val; }
	template <ModuleConcept<R> E0=R> void set(int i, int j, int k, E0 val) { this->data[i][j][k] = val; }


	GenericTensor<R, GenericTensor> pretendFat() const { return GenericTensor<R, GenericTensor>({*this}, dim+1); }
	int dim() const { return this->dim_; }
	int size() const { return this->data.size(); }

	GenericTensor operator*(R c) const { return GenericTensor(map(this->data, [c](M x) { return x * c; }), dim); }
	GenericTensor operator/(R c) const requires DivisionRing<R> { return *this * (1.f / c); }
	GenericTensor operator+(const GenericTensor &T) const { return GenericTensor(map(this->data, T.data, [](M x, M y) { return x + y; }), dim); }
	GenericTensor operator-(const GenericTensor &T) const { return *this + T * -1; }
	GenericTensor operator*(const GenericTensor &T) const { return GenericTensor(map(this->data, T.data, [](M x, M y) { return x * y; }), dim); }
	GenericTensor operator-() const { return *this * -1; }

	template <typename dom>
	Morphism<dom, GenericTensor> switchEvaluationOrder(GenericTensor<R, Morphism<dom, M>> M_f) const {
		return Morphism<dom, GenericTensor>([this, M_f](dom x) { return GenericTensor(map(this->data, [x, M_f](M f) { return M_f(x)(f); }), dim); });
	}

	template <typename dom>
	GenericTensor<R, Morphism<dom, M>> switchEvaluationOrder( Morphism<dom, GenericTensor> f_M) const {
		return Morphism<dom, GenericTensor>([this, f_M](dom x) { return GenericTensor(map(this->data, [x, f_M](M f) { return f_M(x)(f); }), dim); });
	}
};

template <Rng R>
glm::ivec2 mat_size(const GEN_MAT(R) &M) { return glm::ivec2(M.size(), M[0].size()); }



template <Rng R>
GEN_MAT(R) operator*(const GEN_MAT(R) &M1, const GEN_MAT(R) &M2) {
	if (mat_size(M1).y != mat_size(M2).x) throw std::format_error("wrong dimensions of matrices");
	GEN_MAT(R) res = GEN_MAT(R)(mat_size(M1).x, GEN_VEC(R)(mat_size(M2).y, R(0)));
	for (int i = 0; i < mat_size(M1).x; i++)
		for (int j = 0; j < mat_size(M2).y; j++)
			for (int k = 0; k < mat_size(M1).y; k++)
				res.set(i, j, res.at(i, j) +  M1[i][k] * M2[k][j]);
	return res;
}





class Complex {
public: 
	vec2 z;
	float x;
	float y;
	Complex();
	// ReSharper disable once CppNonExplicitConvertingConstructor
	explicit Complex(vec2 z); // NOLINT(*-explicit-constructor)
	Complex(float x, float y);
	// ReSharper disable once CppNonExplicitConvertingConstructor
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

	friend Complex operator*(float f, Complex c) { return c * f; }
	friend Complex operator/(float f, Complex c) { return Complex(f) / c; }
	friend Complex operator+(float f, Complex c) { return c + f; }
	friend Complex operator-(float f, Complex c) { return Complex(f) - c; }

	Complex inv() const { return ~(*this); }
	auto operator==(Complex c) const -> bool;
	auto operator==(float f) const -> bool;
	static auto one() -> Complex;
	static auto zero() -> Complex;
	auto expForm() const -> vec2;
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
	// ReSharper disable once CppNonExplicitConversionOperator
	explicit operator glm::vec2() const;
	explicit operator std::string() const;

	bool nearlyEqual(Complex c) const;
};

auto norm2(Complex c) -> float;
auto abs(Complex c) -> float;


Complex fromExpForm(vec2 e);
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
Complex cot(Complex c);


const Complex ONE = Complex(1, 0);
const Complex ZERO = Complex(0, 0);
const Complex I = Complex(0, 1);

Complex intersectLines(Complex p1, Complex p2, Complex q1, Complex q2);

class Quaternion {
	vec4 q;
public:
	explicit Quaternion(vec4 q) : q(q) {}
	Quaternion(float x, float y, float z, float w) : q(vec4(x, y, z, w)) {}
	explicit Quaternion(float x) : q(vec4(x, 0, 0, 0)) {}
	explicit Quaternion(vec3 im) : q(vec4(0, im)) {}

	Quaternion operator*(Quaternion r) const;
	Quaternion operator*(float f) const { return Quaternion(q*f); }
	Quaternion operator/(float f) const { return *this * (1.f / f); }
	Quaternion operator+(Quaternion r) const { return Quaternion(q + r.q); }
	Quaternion operator-(Quaternion r) const { return *this + r * -1; }
	Quaternion operator-() const { return *this * -1; }
	float norm2() const { return dot(q, q); }
	float norm() const { return std::sqrt(norm2()); }
	Quaternion inv() const { return conj() / norm2(); }
	Quaternion operator~() const { return inv(); }
	Quaternion pow(int p) const;
	Quaternion operator+(float f) const { return *this + Quaternion(f); }
	Quaternion operator-(float f) const { return *this - Quaternion(f); }
	Quaternion operator/(Quaternion r) const { return *this * r.inv(); }
	friend Quaternion operator*(float f, Quaternion q) { return q * f; }
	friend Quaternion operator/(float f, Quaternion q) { return Quaternion(f) * q.inv(); }
	friend Quaternion operator+(float f, Quaternion q) { return q + f; }
	friend Quaternion operator-(float f, Quaternion q) { return Quaternion(f) - q; }
	explicit operator vec4() const { return q; }
	explicit operator std::string() const { return std::format("({0}, {1}, {2}, {3})", q.x, q.y, q.z, q.w); }
	float re() const { return q.w; }
	vec3 im() const { return vec3(q.x, q.y, q.z); }

	Quaternion conj() const { return Quaternion(q.w, -q.x, -q.y, -q.z); }
	Quaternion normalise() const { return *this / norm(); }
	Quaternion rotate(Quaternion r) const { return r * *this * ~r; }
	vec4 rotateS3(vec4 r) const { return vec4(rotate(Quaternion(r))); }

	static Quaternion one() { return Quaternion(1, 0, 0, 0); }
	static Quaternion zero() { return Quaternion(0, 0, 0, 0); }
	static Quaternion i() { return Quaternion(0, 1, 0, 0); }
	static Quaternion j() { return Quaternion(0, 0, 1, 0); }
	static Quaternion k() { return Quaternion(0, 0, 0, 1); }
};

template<DivisionRing F>
class P1 {
public:
	F z;
	bool inf;
	P1() : F(0), inf(false) {}
	explicit P1(F z) : F(z), inf(false) {}
	P1(F z, bool inf) : F(z), inf(inf) {}
	P1(F z1, F z2);

	explicit P1(float x) : F(x), inf(false) {}
	P1(float x, bool inf) : F(x), inf(inf) {}

	bool operator==(P1 c) const {return inf && c.inf || z == c.z && !inf && !c.inf;}
	bool operator==(F c) const {return !inf && z == c;}

	P1 operator+(P1 c) const {return P1(z + c.z, inf || c.inf);}
	P1 operator-() {return P1(-z, inf);}
	P1 operator-(P1 c) const {return (*this) + -c;}
	P1 operator*(P1 c) const;
	P1 operator~() const;
	P1 inv() const {return ~(*this);}
	P1 operator/(P1 c) const {return (*this) * ~c;}
	P1 operator+(F c) const {return (*this) + P1(c);}
	P1 operator-(F c) const {return (*this) - P1(c);}
	P1 operator*(F c) const {return (*this) * P1(c);}
	P1 operator/(F c) const {return (*this) / P1(c);}
	P1 square() const {return *this * *this;}
	P1 pow(int n) const;
	explicit operator F() const;
};

inline float frac(float x) { return x - std::floor(x); }
inline int sgn(float x) { return x < 0 ? -1 : 1; }
inline int sign(float x) { return sgn(x); }

float pseudorandomizer(float x, float seed=0.f);

int binomial(int n, int k);


template <typename M>
float det(const M &m) {return determinant(m);}

template <Normed T>
T normalise(T v) { if (norm(v) > 1e-6) return v / norm(v); return v*0.f; };

template<Normed T>
bool nearlyEqual(T a, T b) { return norm(a - b) < 1e-6; }

inline bool nearlyEqual(float a, int b) { return norm(a - b) < 1e-6; }
inline bool nearlyEqual(int a, float b) { return norm(a - b) < 1e-6; }

template<Normed T>
bool nearlyEqual(T a) { return norm(a) < 1e-6; }


template <VectorSpaceConcept<float> vec>
vec barycenter(vec a, vec b, vec c) { return (a + b + c) / 3.f; };



vec2 intersectLines(vec2 p1, vec2 p2, vec2 q1, vec2 q2);
inline vec2 projectVectorToVector(vec2 v, vec2 n) { return v - dot(v, n) * n; }
inline mat2 scaleMatrix2(vec2 s) { return mat2(s.x, 0, 0, s.y); }
inline mat2 scaleMatrix2(float s) { return scaleMatrix2(vec2(s)); }
inline mat2 changeOfBasis(vec2 t1, vec2 t2) { return inverse(mat2(t1, t2)); }
inline mat2 changeOfBasis(vec2 s1, vec2 s2, vec2 t1, vec2 t2) { return changeOfBasis(t1, t2) * mat2(s1, s2); }
mat2 rotationMatrix2(float angle);
mat2 rotationBetween(vec2 v0, vec2 v1);

inline vec2 orthogonalComplement(vec2 v) { return vec2(-v.y, v.x); }
inline mat2 GramSchmidtProcess(mat2 m) { return mat2( normalise(m[0]), normalise(m[1] - normalise(m[0]) * dot(normalise(m[0]), m[1]))); }

template <VectorSpaceConcept<float> V, VectorSpaceConcept<float> M>
float bilinearForm(V a, V b, M m) {
	float res=0;
	for (int i = 0; i < a.length(); i++)
		for (int j = 0; j < b.length(); j++)
			res += a[i] * m[i][j] * b[j];
	return res; }


std::pair<Complex, Complex> eigenvalues(mat2 m);
inline vec2 eigenvaluesReal(mat2 m) { return vec2(eigenvalues(m).first.re(), eigenvalues(m).second.re()); }
inline bool eigenvectorExists(mat2 m) { return !isClose(m, mat2(0)) && isClose(eigenvalues(m).first.im(), 0); }
bool eigenbasisExists(mat2 m);
vec2 eigenvector(mat2 m, float eigenvalue);
std::pair<vec2, mat2> eigendecomposition(mat2 m);






template <VectorSpaceConcept<float> V, VectorSpaceConcept<float> M>
class EuclideanSpace {
	V value;
	std::shared_ptr<M> metric;
public:
	int dim = value.length();

	EuclideanSpace(V value, std::shared_ptr<M> metric) : value(value), metric(metric) {}
	EuclideanSpace(V value, M&& metric) : value(value), metric(std::make_shared<M>(metric)) {}
	M metricTensor() { return *metric; }

	explicit operator V() { return value; }
	EuclideanSpace operator+(EuclideanSpace other) { return EuclideanSpace(value + other.value, metric); }
	EuclideanSpace operator-(EuclideanSpace other) { return EuclideanSpace(value - other.value, metric); }
	EuclideanSpace operator*(float f) { return EuclideanSpace(value * f, metric); }
	EuclideanSpace operator/(float f) { return EuclideanSpace(value / f, metric); }
	EuclideanSpace operator-() { return EuclideanSpace(-value, metric); }

	friend EuclideanSpace operator*(float f, EuclideanSpace v) { return v * f; }
	friend EuclideanSpace operator/(float f, EuclideanSpace v) { return EuclideanSpace(f / v.value, v.metric); }
	friend float dot(EuclideanSpace a, EuclideanSpace b) { return bilinearForm(a.value, b.value, *a.metric); }

	float norm2() { return dot(*this, *this); }
	float norm() { return sqrt(norm2()); }
	EuclideanSpace normalise() { return *this / norm(); }
	EuclideanSpace orthogonalProjection(EuclideanSpace v) { return v - normalise(*this) * dot(v, *this); }
	M GSProcess(const M &basis);

	friend float norm(EuclideanSpace v) { return v.norm(); }
	friend float norm2(EuclideanSpace v) { return v.norm2(); }
	friend EuclideanSpace normalise(EuclideanSpace v) { return v / norm(v); }

	vec2 toVec2() { if (dim != 2) throw std::format_error("wrong dimension of vector (not 2)"); return vec2(value); }
	vec3 toVec3() { if (dim != 3) throw std::format_error("wrong dimension of vector (not 2)"); return vec3(value); }
	vec4 toVec4() { if (dim != 4) throw std::format_error("wrong dimension of vector (not 2)"); return vec4(value); }

	std::vector<Complex> eigenvalues(M m) {
		if (dim != 2) throw std::format_error("eigendecomposition currently implemented only in dimension 2");
		return eigenvalues(static_cast<mat2>(m)); }

	M orthogonalEigenbasis() {
        if (dim != 2) throw std::format_error("eigendecomposition currently implemented only in dimension 2");
        return GSProcess(static_cast<M>(eigendecomposition(static_cast<mat2>(metric)).second)); }

	bool orthogonalEigenbasisExists() {
        if (dim != 2) throw std::format_error("eigendecomposition currently implemented only in dimension 2");
        return eigenbasisExists(static_cast<mat2>(metric)); }
};

vec3 projectVectorToPlane(vec3 v, vec3 n);
mat3 scaleMatrix3(vec3 s);
mat3 scaleMatrix3(float s);
inline mat3 changeOfBasis(vec3 t1, vec3 t2, vec3 t3) {return inverse(mat3(t1, t2, t3));}
inline mat3 changeOfBasis(vec3 s1, vec3 s2, vec3 s3, vec3 t1, vec3 t2, vec3 t3) {return changeOfBasis(t1, t2, t3) * mat3(s1, s2, s3);}
mat3 rotationMatrix3(float angle);
mat3 rotationMatrix3(vec3 axis, float angle);
mat3 rotationBetween(vec3 v0, vec3 v1);

std::pair<vec3, vec3> orthogonalComplementBasis(vec3 v);
mat3 GramSchmidtProcess(mat3 m);


inline float min(float a, float b) { return a < b ? a : b; }
inline float max(float a, float b) { return a > b ? a : b; }

template<typename V>
V lerp(V a, V b, float t) { return b * t + a * (1.f - t); }
template<typename V>
V clamp(V a, V b, float t) { return min(b, max(a, t)); }
template<typename V>
V lerpClamp(V a, V b, float t, V max_v = V(1), V min_v=V(0)) { return lerp(a, b, clamp(min_v, max_v, t)); }

template<typename V>
V mean(std::vector<V> points) {
	V m = V(0);
	for (V p : points)
		m = m + p;
	return m / (1.f*points.size()); }


template<VectorSpaceConcept<float> V>
vector<V> linspace(V a, V b, int n, bool includeEnd=true) {
	vector<V> res;
	res.reserve(n);
	if (includeEnd)
		for (int i = 0; i < n; i++)
			res.push_back(lerp(a, b, i*1.f / (n - 1.f)));
	else
		for (int i = 0; i < n; i++)
			res.push_back(lerp(a, b, i*1.f / n));
	return res; }

template<AbelianSemigroup V>
vector<vector<V>> linspace2D(V p0, V v0, V v1, int n, int m) {
	vector<vector<V>> res;
	res.reserve(n);
	for (int i = 0; i < n; i++) {
		res.push_back(vector<V>());
		res[i].reserve(m);
		for (int j = 0; j < m; j++)
			res[i].push_back(p0 + v0 * i + v1 * j);
	}
	return res;
}



vector<int> range(int a, int b, int step=1);
inline vector<int> range(int n) { return range(0, n); }

template<TotallyOrderedAbelianSemigroup A>
vector<A> arange(A a, A b, A step) {
	vector<A> res;
	for (A i = a; i < b; i += step) {
		res.push_back(i);
		if (res.size() > 1000000) throw std::format_error("arange didn't terminate in milion iterations");
	}
	return res; }

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




template<RingConcept T, int n, int m>
Matrix<T, n, m>::Matrix() {
	this->coefs = std::array<std::array<T, n>, m>();
	for (int i = 0; i < m; i++)
		for (int j            = 0; j < n; j++)
			this->coefs[i][j] = 0;
}

template<RingConcept T, int n, int m>
Matrix<T, n, m>::Matrix(std::array<std::array<T, n>, m> c) {
	this->coefs = c;
}

template<RingConcept T, int n, int m>
Matrix<T, n, m>::Matrix(vector<vector<T>> c) {
	this->coefs = std::array<std::array<T, n>, m>();
	for (int i = 0; i < m; i++)
		for (int j            = 0; j < n; j++)
			this->coefs[i][j] = c[i][j];
}

template <RingConcept T, int n, int m>
Matrix<T, n, m> Matrix<T, n, m>::operator*(const T& f) const
{
	Matrix result = Matrix();
	for (int i = 0; i < m; i++)
		for (int j             = 0; j < n; j++)
			result.coefs[i][j] = this->coefs[i][j] * f;
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
	for (int i = 0; i < m; i++)
		for (int j             = 0; j < n; j++)
			result.coefs[i][j] = this->coefs[j][i];
	return result;
}

template <RingConcept T, int n, int m>
Matrix<T, n, m> Matrix<T, n, m>::operator-() const
{
	return *this * -1;
}



template<RingConcept T, int n>
Matrix<T, n, n>::Matrix(T diag): Matrix() {
	for (int i            = 0; i < n; i++)
		this->coefs[i][i] = diag;
}

template<RingConcept T, int n>
T Matrix<T, n>::minor(int i, int j) const {
	if (n <= 1)
		throw std::invalid_argument("Matrix must be at least 2x2");
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
	for (int i = 0; i < n; i++)
		for (int j = 0; j < k; j++) {
			T sum = 0;
			for (int l = 0; l < m; l++)
				sum += A.coefs[i][l] * B.coefs[l][j];
			result.coefs[i][j] = sum;
		}
	if (n == 2 && k == 2) {
		result.a = result.coefs[0][0];
		result.b = result.coefs[0][1];
		result.c = result.coefs[1][0];
		result.d = result.coefs[1][1];
	}
	return result;
}

template<Rng R, ModuleConcept<R> M>
template<ModuleConcept<R> E0>
E0 GenericTensor<R, M>::at(vector<int> indices) const {
	if (indices.size() == 1) return this->data[indices[0]];
	vector<int> reduced_indices = indices;
	reduced_indices.erase(reduced_indices.begin());
	return data[indices[0]].multiindex(reduced_indices);
}

template<Rng R, ModuleConcept<R> M>
template<ModuleConcept<R> E0>
void GenericTensor<R, M>::set(vector<int> ind, E0 val) {
	if (ind.size() == 1) this->data[ind[0]] = val;
	else if (ind.size() == 2) this->data[ind[0]][ind[1]] = val;
	else if (ind.size() == 3) this->data[ind[0]][ind[1]][ind[2]] = val;
	else {
		vector<int> ind2 = vector(ind.begin() + 1, ind.end());
		this->data[ind[0]].set(ind2, val);
	}
}

template <RingConcept T, int n, int m>
Matrix<T, n, m> Matrix<T, n, m>::operator/(const T& f) const requires DivisionRing<T> {
	Matrix<T, n, m> result = Matrix<T, n, m>();
	for (int i = 0; i < m; i++)
		for (int j             = 0; j < n; j++)
			result.coefs[i][j] = this->coefs[i][j] / f;
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
template<DivisionRing F>
P1<F>::operator F() const {
	if (inf) throw ValueError("cannot convert infinity to field element");
	return z;
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
	for (int i = 0; i < n; i++)
		tr     = tr +  this->coefs[i][i];
	return tr;
}

template <RingConcept T, int n, int m>
Matrix<T, n, m>::operator std::string() const
{
	std::string s = "";
	for (int i = 0; i < m; i++) {
		s = s + "|";
		for (int j = 0; j < n; j++)
			s      = s + std::string(this->coefs[i][j]);
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
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++) {
			int sign      = ((i + j) % 2 == 0) ? 1 : -1;
			I.coefs[j][i] = this->minor(i, j) * sign;
		}
	return I / this->det();
}

template<RingConcept T, int n>
Matrix<T, n> Matrix<T, n, n>::GramSchmidtProcess() const requires EuclideanSpaceConcept<T> {
    Matrix result;
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
Matrix<T, 2, 2>::Matrix(T a, T b, T c, T d) {
	this->a     = a;
	this->b     = b;
	this->c     = c;
	this->d     = d;
	this->coefs = std::array<std::array<T, 2>, 2>({ { {a, b}, {c, d} } });
}

template<RingConcept T>
Matrix<T, 2, 2>::Matrix(std::array<std::array<T, 2>, 2> c) {
	this->a     = c[0][0];
	this->b     = c[0][1];
	this->c     = c[1][0];
	this->d     = c[1][1];
	this->coefs = c;
}

template<RingConcept T>
Matrix<T, 2, 2>::Matrix(vector<vector<T>> c) {
	this->coefs = std::array<std::array<T, 2>, 2>();
	for (int i = 0; i < 2; i++)
		for (int j            = 0; j < 2; j++)
			this->coefs[i][j] = c[i][j];
	this->a = c[0][0];
	this->b = c[0][1];
	this->c = c[1][0];
	this->d = c[1][1];
}


template<RingConcept T>
T Matrix<T, 2>::mobius_derivative(T z) const requires DivisionRing<T> {
	return (a * d - b * c) / ((c * z + d) * (c * z + d));
}

template<RingConcept T>
Matrix<T, 2> Matrix<T, 2, 2>::pow(int p) const {
	if (p < 0)
		return ~(*this).pow(-p);
	if (p == 0)
		return Matrix(T(1));
	if (p == 1)
		return *this;
	if (p % 2 == 0)
		return this->pow(p / 2) * this->pow(p / 2);
	return *this * this->pow(p / 2) * this->pow(p / 2);
}









template<int n, int m>
bool nearlyEqual(Matrix<Complex, n, m> a, Matrix<Complex, n, m> b) {
	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++)
			if (!a.at(i, j).nearlyEqual(b.at(i, j)))
				return false;
	return true;
}

inline bool nearlyEqual(mat3 a, mat3 b) {
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			if (abs(a[i][j] - b[i][j]) > 1e-6)
				return false;
	return true;
}


template<DivisionRing F>
P1<F>::P1(F z1, F z2) {
	if (z2 == F(0)) {
		if (z1 == F(0)) throw ValueError("0/0 not well defined");
		z   = z1;
		inf = true;
	}
	else {
		z   = z1/z2;
		inf = false;
	}
}

template<DivisionRing F>
P1<F> P1<F>::operator*(P1 c) const {
	if (inf && c==F(0)) throw ValueError("0*inf is undefined");
	if (z == F(0) && c.inf) throw ValueError("0*inf is undefined");
	return P1(z * c.z, inf || c.inf);
}

template<DivisionRing F>
P1<F> P1<F>::operator~() const {
	if (inf) return P1(0);
	if (z == 0) return P1(0, true);
	return P1(1/z);}

template<DivisionRing F>
P1<F> P1<F>::pow(int n) const {return n == 0 ? P1(1) : n == 1 ? *this : n == -1 ? inv() : square().pow(n/2);}

template<VectorSpaceConcept<float> V, VectorSpaceConcept<float> M>
M EuclideanSpace<V, M>::GSProcess(const M &basis) {
	M result = basis;
    for (int i = 0; i < dim; i++) {
        EuclideanSpace v = EuclideanSpace(basis[i], metric);
    	for (int j = 0; j < i; j++) {
    		EuclideanSpace u = EuclideanSpace(result[j], metric);
    		v = v - normalise(u) * dot(v, u);
    	}
        result[i] = static_cast<V>(normalise(v));
    }
    return result;
}



template <typename T>
vector<T> reverse(vector<T> v) {
	std::reverse(v.begin(), v.end());
	return v;
}

template <AbelianMonoid V>
vector<V> arange(V a, V b, V step) {
	if (step == V(0)) throw std::invalid_argument("step cannot be 0");

	vector<V> res = vector<V>();
	res.reserve(abs(b - a) / abs(step));

	if (step > 0) {
		for (V i = a; i < b; i += step)
			res.push_back(i);
		return res; }

	for (V i = b; i > a; i += step)
		res.push_back(i);
	return res;
}

template <typename T>
vector<T> smartRange(vector<T> v, int a, int b, int step=1) {
	if (step == 0) throw std::invalid_argument("step cannot be 0");

	vector<T> res = vector<T>();
	res.reserve(abs(b - a) / abs(step));

	if (step > 0) {
		for (int i = a; i < b; i += step)
			res.push_back(v[i]);
		return res; }

	for (int i = b; i > a; i += step)
		res.push_back(v[i]);
	return res;
}

template <typename T>
vector<T> range(vector<T> v, int a, int b, int step=1) {
	return smartRange(v, a, b, step);
}

template <typename T>
vector<T> rangeFrom(vector<T> v, int a) {
	return smartRange(v, a, v.size());
}

template <typename T>
vector<T> rangeTo(vector<T> v, int a) {
	return smartRange(v, 0, a);
}

template <typename T>
vector<T> rangeStep(vector<T> v, int s) {
	return smartRange(v, 0, v.size(), s);
}

float randomUniform(float a, float b);
vec2 randomUniform(vec2 a, vec2 b);
vec3 randomUniform(vec3 a, vec3 b);

inline float angle(vec3 v1, vec3 v2) {
	return atan2(length(cross(v1, v2)), dot(v1, v2));
}

float rotationAngle(const mat3 &M);
float rotationCosAngle(mat3 M);

vec3 rotationAxis(mat3 M);
mat3 spinTensor(vec3 omega);
mat3 rotationMatrix(vec3 axis, float angle);
mat3 rotationMatrix(vec3 omega);

float polarAngle(vec3 v, vec3 t1, vec3 t2);
float polarAngle(vec3 v, vec3 n);

float cot(float x);
inline float square(float x) { return x * x; }
inline float cube(float x) { return x * x * x; }
inline float pow4(float x) { return x * x * x * x; }
inline float pow5(float x) { return x * x * x * x * x; }
inline float pow6(float x) { return x * x * x * x * x * x; }
inline float pow7(float x) { return x * x * x * x * x * x * x; }
inline float pow8(float x) { return x * x * x * x * x * x * x * x; }
inline float pow3(float x) { return x * x * x; }
inline float pow2(float x) { return x * x; }
inline float sq(float x) { return x * x; }
