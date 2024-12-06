#include "mat.hpp"

#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <set>
#include <sstream>
#include <vector>
#include <GL/glew.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

using namespace glm;
using std::vector, std::string, std::set, std::array, std::pair, std::exp, std::log, std::cos, std::sin, std::cosh, std::sinh, std::sqrt, std::pow, std::atan2, std::abs;


BigMatrix::operator tvec2<float>() { if (std::min(n(), m()) != 1 || std::max(m(), n()) != 2) throw std::format_error("wrong dimension of matrix (" + std::to_string(n()) + ", " + std::to_string(m()) + ")" );
	return (*this)[0].size() > 1 ? vec2((*this)[0][0], (*this)[0][1]) : vec2((*this)[0][0], (*this)[1][0]); }

BigMatrix::operator tvec3<float>() { if (std::min(n(), m()) != 1 || std::max(m(), n()) != 3) throw std::format_error("wrong dimension of matrix (" + std::to_string(n()) + ", " + std::to_string(m()) + ")" );
	return (*this)[0].size() > 1 ? vec3((*this)[0][0], (*this)[0][1], (*this)[0][2]) : vec3((*this)[0][0], (*this)[1][0], (*this)[2][0]); }

BigMatrix::operator tvec4<float>() { if (std::min(n(), m()) != 1 || std::max(m(), n()) != 4) throw std::format_error("wrong dimension of matrix (" + std::to_string(n()) + ", " + std::to_string(m()) + ")" );
	return (*this)[0].size() > 1 ? vec4((*this)[0][0], (*this)[0][1], (*this)[0][2], (*this)[0][3]) : vec4((*this)[0][0], (*this)[1][0], (*this)[2][0], (*this)[3][0]); }

Complex::Complex() {
	z = vec2(0, 0);
	x = 0;
	y = 0;
}

Complex::Complex(vec2 z)
{
	this->z = z;
	this->x = z.x;
	this->y = z.y;
}

Complex::Complex(float x, float y)
{
	z = vec2(x, y);
	this->x = x;
	this->y = y;
}

Complex::Complex(float x)
{
	z = vec2(x, 0);
	this->x = x;
	this->y = 0;
}

Complex Complex::operator+(Complex c) const {
	return Complex(z + c.z);
}

Complex Complex::operator-(const Complex c) const {
	return Complex(z - c.z);
}



Complex Complex::operator*(Complex c) const {
	return Complex(vec2(z.x * c.z.x - z.y * c.z.y, z.x * c.z.y + z.y * c.z.x));
}


Complex Complex::conj() const {
	return Complex(vec2(z.x, -z.y));
}

Complex Complex::operator~() const {
	return Complex(z.x / norm2(z), -z.y / norm2(z));
}

void Complex::operator+=(Complex c)
{
	z += c.z;
}

void Complex::operator-=(Complex c)
{
	z -= c.z;
}

void Complex::operator*=(Complex c)
{
	z = vec2(z.x * c.z.x - z.y * c.z.y, z.x * c.z.y + z.y * c.z.x);
}	

void Complex::operator/=(Complex c)
{
	*this *= c.inv();
}




Complex Complex::operator-() const {
	return Complex(-z);
}


float Complex::arg() const {
	return atan2(z.y, z.x);
}

float Complex::re() const {
	return z.x;
}

float Complex::im() const {
	return z.y;
}




float norm(Complex z)
{
	return norm(z.z);
}

std::ostream& operator<<(std::ostream& _stream, Complex const& z)
{
	_stream << z.z.x << " + " << z.z.y << "i";
	return _stream;
}

Complex toComplex(vec2 v)
{
	return Complex(v);
}

float abs(Complex z)
{
	return norm(z.z);
}

float norm2(Complex z)
{
	return norm2(z.z);
}




Complex Complex::operator/(Complex c) const {
	return c.inv() * (*this);
}


Complex Complex::square() const {
	return (*this) * (*this);
}

vec2 Complex::expForm() const {
	return vec2(norm(z), arg());
}

Complex fromExpForm(vec2 e)
{
	return Complex(e.x * cos(e.y), e.x * sin(e.y));
}

bool Complex::operator==(Complex c) const {
	return z == c.z;
}

bool Complex::operator==(float f) const {
	return z == vec2(f, 0);
}

Complex Complex::one()
{
    return Complex(1, 0);
}

Complex Complex::zero()
{
	return Complex(0, 0);
}

Complex Complex::pow(float exponent) const
{
	if (exponent == 0)
	{
		return Complex(1, 0);
	}
	if (exponent == 1)
	{
		return *this;
	}
	if (exponent == 2)
	{
		return square();
	}
	if (exponent == 0.5)
	{
		return sqrt();
	}
	if (exponent == -1)
	{
		return inv();
	}
	if (exponent == -2)
	{
		return inv().square();
	}
	if (norm(z) == 0)
	{
		return Complex(0, 0);
	}
	return fromExpForm(vec2(std::pow(norm(z), exponent), arg() * exponent));
}

Complex Complex::sqrt() const {
	if (z.x == 0 && z.y == 0)
	{
		return Complex(0, 0);
	}
	return fromExpForm(vec2(std::sqrt(norm(z)), arg() / 2));
}

Complex::operator glm::vec2() const {
	return z;
}

Complex::operator std::string() const {
	std::ostringstream out;
    out.precision(2);
    out << std::fixed << z.x << " + " << z.y << "i";
    return std::move(out).str();

}

bool Complex::nearlyEqual(Complex c) const {
	return abs(norm(z - c.z)) < 1e-6;
}


Complex rootOfUnity(int n, int k)
{
	return Complex(std::cos(k*TAU / n), std::sin(k*TAU / n));
}

vector<Complex> allRootsOf1(int n)
{
	std::vector<Complex> roots;
	for (int k = 0; k < n; k++)
	{
		roots.push_back(rootOfUnity(n, k));
	}
	return roots;
}

Complex exp(Complex c)
{
	return Complex(cos(c.y), sin(c.y))*std::exp(c.x);
}

Complex log(Complex c)
{
	return Complex(std::log(norm(c.z)), c.arg());
}

Complex sin(Complex c)
{
	return Complex(sin(c.real()) * cosh(c.imag()), cos(c.real()) * sinh(c.imag()));
}

Complex cos(Complex c)
{
	return Complex(cos(c.real()) * cosh(c.imag()), -sin(c.real()) * sinh(c.imag()));
}

Complex tan(Complex c)
{
	return sin(c) / cos(c);
}

Complex sinh(Complex c)
{
	return Complex(sinh(c.real()) * cos(c.imag()), cosh(c.real()) * sin(c.imag()));
}

Complex cosh(Complex c)
{
	return Complex(cosh(c.real()) * cos(c.imag()), sinh(c.real()) * sin(c.imag()));
}

Complex tanh(Complex c)
{
	return sinh(c) / cosh(c);
}

Complex asin(Complex c)
{
	return I * log(I * c + (ONE - c.square()).sqrt())*(-1);
}

Complex acos(Complex c)
{
	return I * log(c + I * (ONE - c.square()).sqrt()) * (-1);
}

Complex atan(Complex c)
{
	return I / 2 * (log(ONE - I * c) - log(ONE + I * c));
}

Complex asinh(Complex c)
{
	return log(c + (c.square() + ONE).sqrt());
}

Complex acosh(Complex c)
{
	return log(c + (c + ONE).sqrt() * (c - ONE).sqrt());
}

Complex atanh(Complex c)
{
	return (log(ONE + c) - log(ONE - c))/2;
}

CP1::CP1()
{
	z = Complex(0, 0);
	inf = false;
}

CP1::CP1(Complex z)
{
	this->z = z;
	inf = false;
}

CP1::CP1(Complex z, bool inf)
{
	this->z = z;
	this->inf = inf;
}

CP1::CP1(Complex z1, Complex z2)
{
	if (z2 == 0)
	{
		this->z = Complex();
		this->inf = true;
	}
	else
	{
		this->z = z1 / z2;
		this->inf = false;
	}
}

CP1::CP1(float x)
{
	z = Complex(x, 0);
	inf = false;
}

CP1::CP1(float x, bool inf)
{
	z = Complex(x, 0);
	this->inf = inf;
}

CP1 CP1::operator+(CP1 c) const {
	if (inf || c.inf)
	{
		return CP1(Complex(1, 0), true);
	}
	return CP1(z + c.z);
}

CP1 CP1::operator-(CP1 c) const {
	if (inf || c.inf)
	{
		return CP1(Complex(1, 0), true);
	}
	return CP1(z - c.z);
}

CP1 CP1::operator*(CP1 c) const {
	if (inf || c.inf)
	{
		return CP1(Complex(1, 0), true);
	}
	return CP1(z * c.z);
}

CP1 CP1::operator/(CP1 c) const {
	if (abs(c.z) == 0.f && !c.inf)
	{
		return CP1(Complex(1, 0), true);
	}
	if (c.inf)
	{
		return CP1(Complex(0, 0), false);
	}
	if (inf)
	{
		return CP1(Complex(1, 0), true);
	}
	return CP1(z / c.z, false);
}

CP1 CP1::operator+(Complex c) const {
	if (inf)
	{
		return CP1(Complex(1, 0), true);
	}
	return CP1(z + c);
}

CP1 CP1::operator-(Complex c) const {
	if (inf)
	{
		return CP1(Complex(1, 0), true);
	}
	return CP1(z - c);
}

CP1 CP1::operator*(Complex c)
{
	if (inf)
	{
		return CP1(Complex(1, 0), true);
	}
	return CP1(z * c);
}

CP1 CP1::operator/(Complex c)
{
	if (abs(c) == 0.f)
	{
		return CP1(Complex(1, 0), true);
	}
	if (inf)
	{
		return CP1(Complex(0, 0), false);
	}
	return CP1(z / c, false);
}

CP1 CP1::operator+(float f)
{
	if (inf)
		return CP1(Complex(1, 0), true);
	return CP1(z + f);
}

CP1 CP1::operator-(float f)
{
	if (inf)
	{
		return CP1(Complex(1, 0), true);
	}
	return CP1(z -f);
}

CP1 CP1::operator*(float f)
{
	if (inf)
	{
		return CP1(Complex(1, 0), true);
	}
	return CP1(z * f);
}

CP1 CP1::operator/(float f)
{
	if (abs(f) == 0.f)
	{
		return CP1(Complex(1, 0), true);
	}
	if (inf)
	{
		return CP1(Complex(0, 0), false);
	}
	return CP1(z / f);
}

CP1 CP1::operator-()
{
	return CP1(-z, inf);
}

CP1 CP1::inv()
{
	if (inf)
	{
		return CP1(Complex(0, 0), false);
	}
	if (abs(z) == 0.f)
	{
		return CP1(Complex(1, 0), true);
	}
	return CP1(z.inv(), false);
}

CP1 CP1::square()
{
	if (inf)
	{
		return CP1(Complex(1, 0), true);
	}
	return CP1(z.square());
}

CP1 CP1::sqrt()
{
	if (inf)
	{
		return CP1(Complex(1, 0), true);
	}
	return CP1(z.sqrt());
}

CP1 CP1::pow(float f)
{
	if (inf)
	{
		return CP1(Complex(1, 0), true);
	}
	return CP1(z.pow(f));
}

CP1 CP1::conj()
{
	if (inf)
	{
		return CP1(Complex(1, 0), true);
	}
	return CP1(z.conj());
}

CP1::operator Complex()
{
	if (inf)
	{
		std::cout << "Trying to cast infinity to complex number" << std::endl;
		return Complex(NAN, NAN);
	}
	return z;
}

CP1 oneCP1()
{
	return CP1(ONE);
}

CP1 zeroCP1()
{
	return CP1(ZERO);
}

CP1 infCP1()
{
	return CP1(Complex(0, 0), true);
}

CP1 iCP1()
{
	return CP1(I);
}





vec2 intersectLines(vec2 p1, vec2 p2, vec2 q1, vec2 q2)
{
	return (p1 * (q1.y - q2.y) - p2 * (q1.y - q2.y) - q1 * (p1.y - p2.y) + q2 * (p1.y - p2.y)) / ((p1.x - p2.x) * (q1.y - q2.y) - (p1.y - p2.y) * (q1.x - q2.x));
}


Complex intersectLines(Complex p1, Complex p2, Complex q1, Complex q2)
{
	return Complex(intersectLines(p1.z, p2.z, q1.z, q2.z));
}

mat3 scaleMatrix3(vec3 s) {
    return mat3(s.x, 0, 0, 0, s.y, 0, 0, 0, s.z);
}

mat3 scaleMatrix3(float s) {
	return scaleMatrix3(vec3(s, s, s));
}

mat3 rotationMatrix3(float angle) {
    return mat3(cos(angle), -sin(angle), 0, sin(angle), cos(angle), 0, 0, 0, 1);
}

mat2 rotationMatrix2(float angle) {
	return mat2(cos(angle), -sin(angle), sin(angle), cos(angle));
}


mat3 rotationMatrix3(vec3 axis, float angle)
{
    auto planeBasis = orthogonalComplementBasis(axis);
	mat3 change = changeOfBasis(planeBasis.first, planeBasis.second, axis);
	return inverse(change) * rotationMatrix3(angle) * change;
}

mat3 rotationBetween(vec3 v0, vec3 v1) {
    vec3 axis = normalise(cross(v0, v1));
    float angle = acos(dot(v0, v1) / (norm(v0) * norm(v1)));
    return rotationMatrix3(axis, angle);
}

mat2 rotationBetween(vec2 v0, vec2 v1) {
	float angle = acos(dot(v0, v1) / (norm(v0) * norm(v1)));
	return rotationMatrix2(angle);
}


mat3 GramSchmidtProcess(mat3 m) {
    vec3 u1 = normalise(m[0]);
    vec3 u2 = normalise(m[1] - u1 * dot(u1, m[1]));
    vec3 u3 = normalise(m[2] - u1 * dot(u1, m[2]) - u2 * dot(u2, m[2]));
	mat3 ma = mat3(u1, u2, u3);
    return ma;
}

vector<int> range(int a, int b, int step) {
	if (step == 0) throw std::invalid_argument("step cannot be 0");

	vector<int> res = vector<int>();
	res.reserve(abs(b - a) / abs(step));

	if (step > 0) {
		for (int i = a; i < b; i += step)
			res.push_back(i);
		return res; }

	for (int i = b; i > a; i += step)
		res.push_back(i);
	return res;
}

float randomUniform(float a, float b) {
	return a + (b - a) * (rand() / (RAND_MAX + 1.0));
}

vec2 randomUniform(vec2 a, vec2 b) {
	return vec2(randomUniform(a.x, b.x), randomUniform(a.y, b.y));
}

vec3 randomUniform(vec3 a, vec3 b) {
	return vec3(randomUniform(a.x, b.x), randomUniform(a.y, b.y), randomUniform(a.z, b.z));
}

float rotationAngle(const mat3 &M) {
	return acos(rotationCosAngle(M));
}

float rotationCosAngle(mat3 M) {
	return (M[0][0] + M[1][1] + M[2][2] - 1) / 2;
}

vec3 rotationAxis(mat3 M) {
	float sinAngle = sqrt(1-pow(rotationCosAngle(M), 2));
	return vec3(M[2][1] - M[1][2], M[0][2] - M[2][0], M[1][0] - M[0][1]) / (2 * sinAngle);
}

mat3 spinTensor(vec3 omega) {
	return mat3(0, -omega.z, omega.y, omega.z, 0, -omega.x, -omega.y, omega.x, 0);
}

mat3 rotationMatrix(vec3 axis, float angle) {
	vec3 n   = normalise(axis);
	float nx = n.x;
	float ny = n.y;
	float nz = n.z;

	return mat3(
			cos(angle) + nx * nx * (1 - cos(angle)), nx * ny * (1 - cos(angle)) - nz * sin(angle), nx * nz * (1 - cos(angle)) + ny * sin(angle),
			ny * nx * (1 - cos(angle)) + nz * sin(angle), cos(angle) + ny * ny * (1 - cos(angle)), ny * nz * (1 - cos(angle)) - nx * sin(angle),
			nz * nx * (1 - cos(angle)) - ny * sin(angle), nz * ny * (1 - cos(angle)) + nx * sin(angle), cos(angle) + nz * nz * (1 - cos(angle)));
}

mat3 rotationMatrix(vec3 omega) { return rotationMatrix(normalise(omega), length(omega)); }


vec3 projectVectorToPlane(vec3 v, vec3 n) {
	return v - dot(v, n) * n;
}

std::pair<vec3, vec3> orthogonalComplementBasis(vec3 v)
{
    vec3 b1 = norm(cross(v, normalise(e1-e2+e3))) > norm(cross(v, normalise(e1+e2))) ?  normalise(e1-e2+e3) : normalise(e1+e2);
    vec3 b2 = cross(normalise(v), b1);
    mat3 frame = GramSchmidtProcess(mat3(normalise(v), b1, b2 ));
    return std::make_pair(frame[1], frame[2]);

}

float pseudorandomizer(float x, float seed) { return frac(sin(x + seed) * 43758.5453f + seed); }


SparseMatrix::SparseMatrix(int n, int m) {
	this->n    = n;
	this->m    = m;
	this->data = std::vector<std::vector<std::pair<int, float>>>(n);
}

void SparseMatrix::set(int i, int j, float val) {
	this->data[i].emplace_back(j, val);
}

float SparseMatrix::get(int i, int j) {
	for (auto p : this->data[i])
		if (p.first == j)
			return p.second;
	return 0;
}

SparseMatrix SparseMatrix::operator*(float f) {
	SparseMatrix result = SparseMatrix(this->n, this->m);
	for (int i = 0; i < this->n; i++)
		for (auto p : this->data[i])
			result.set(i, p.first, p.second * f);
	return result;
}

SparseMatrix SparseMatrix::operator+(SparseMatrix M) {
	SparseMatrix result = SparseMatrix(this->n, this->m);
	for (int i = 0; i < this->n; i++)
		for (auto p : this->data[i])
			result.set(i, p.first, p.second);
	for (int i = 0; i < M.n; i++)
		for (auto p : M.data[i])
			result.set(i, p.first, result.get(i, p.first) + p.second);
	return result;
}

SparseMatrix SparseMatrix::operator-(SparseMatrix M) {
	SparseMatrix result = SparseMatrix(this->n, this->m);
	for (int i = 0; i < this->n; i++)
		for (auto p : this->data[i])
			result.set(i, p.first, p.second);
	for (int i = 0; i < M.n; i++)
		for (auto p : M.data[i])
			result.set(i, p.first, result.get(i, p.first) - p.second);
	return result;
}

BigVector::BigVector(const std::vector<std::vector<float>> &data) {
	this->data = data[0];
	for (int i = 1; i < data.size(); i++)
		for (float j : data[i])
			this->data.push_back(j);
}

BigVector::BigVector(int n, float val) {
	data = std::vector<float>();
	data.reserve(n);
	for (int i = 0; i < n; i++)
		this->data.push_back(val);
}

BigVector & BigVector::operator=(const BigVector &other) {
	if (this == &other)
		return *this;
	data = other.data;
	return *this;
}

BigVector & BigVector::operator=(BigVector &&other) noexcept {
	if (this == &other)
		return *this;
	data = std::move(other.data);
	return *this;
}

BigVector BigVector::operator*(float f) const {
	std::vector<float> result = std::vector<float>(this->data.size());
	for (int i    = 0; i < this->data.size(); i++)
		result[i] = this->data[i] * f;
	return BigVector(result);
}

BigVector BigVector::operator+(const BigVector& v) const{
	std::vector<float> res;
	for (int i = 0; i < this->data.size(); i++)
		res.push_back(this->data[i] + v[i]);
	return BigVector(res);
}

BigVector BigVector::operator-(const BigVector& v)const {
	std::vector<float> res;
	for (int i = 0; i < this->data.size(); i++)
		res.push_back(this->data[i] - v[i]);
	return BigVector(res);
}



BigMatrix::BigMatrix(int n, int m) { this->data = MATR$X(n, vec2137(m, 0)); }

BigMatrix BigMatrix::submatrix(int i, int j) {
    MATR$X sub = MATR$X();
    for (int k = 0; k < n(); k++) {
        if (k == i) continue;
        std::vector<float> row = std::vector<float>();
        for (int l = 0; l < m(); l++) {
            if (l == j) continue;
            row.push_back(this->data[k][l]);
        }
        sub.push_back(row);
    }
    return BigMatrix(sub);
}

BigMatrix BigMatrix::diagonalComponent() const { BigMatrix res = BigMatrix(this->n(), this->m()); for (int i = 0; i < this->n(); i++) res.set(i, i, this->get(i, i)); return res; }
BigMatrix BigMatrix::invertedDiagonal() const { BigMatrix res = BigMatrix(this->n(), this->m()); for (int i = 0; i < this->n(); i++) res.set(i, i, 1/this->get(i, i)); return res; }

int binomial(int n, int k) {
	if (k > n) return 0;
	if (k == 0 || k == n) return 1;
	return binomial(n - 1, k - 1) + binomial(n - 1, k);
}

BigMatrix::BigMatrix(const MATR$X &data) { this->data = data; };

float BigMatrix::det() {
    if (n() != m())  throw std::invalid_argument("Matrix must be square");
    if (n() == 1)  return this->data[0][0];
    float result = 0;
    for (int j = 0; j < m(); j++)  if (this->data[0][j] != 0)
        result += this->data[0][j] * submatrix(0, j).det() * (j % 2 == 0 ? 1 : -1);
    return result;
}


    BigMatrix BigMatrix::transpose() const {
    MATR$X data = MATR$X();
    data.reserve(n());
    for (int i = 0; i < n(); i++) data[i].reserve(m());
    for (int i = 0; i < m(); i++) {
        std::vector<float> row = std::vector<float>();
        for (int j = 0; j < n(); j++)
            row.push_back(this->data[j][i]);
        data.push_back(row);
    }
	return BigMatrix(data);
}


    BigMatrix BigMatrix::operator*(float f) const {
            BigMatrix result = BigMatrix(this->data);
            for (int j = 0; j < n(); j++)
                for (int i = 0; i < m(); i++)
                    result.data[j][i] *= f;
        return result;
    }
    BigMatrix BigMatrix::operator+(const BigMatrix &M)  const{
    BigMatrix result = BigMatrix(this->data);
    for (int j = 0; j < n(); j++)
        for (int i = 0; i < m(); i++)
            result.data[j][i] += M.data[j][i];
    return result;
}

BigMatrix BigMatrix::operator-(const BigMatrix &M) const {
    BigMatrix result = BigMatrix(this->data);
    for (int j = 0; j < n(); j++)
        for (int i = 0; i < m(); i++)
            result.data[j][i] -= M.data[j][i];
    return result;
}

BigMatrix BigMatrix::operator*(const BigMatrix &M) const {
    if (m() != M.n())  throw std::invalid_argument("Matrix dimensions must agree");
    BigMatrix result = BigMatrix(n(), M.m());
    for (int i = 0; i < n(); i++)
        for (int j = 0; j < M.m(); j++)
            for (int k = 0; k < m(); k++)
                result.data[i][j] += this->data[i][k] * M.data[k][j];
    return result;
}

BigMatrix BigMatrix::operator*(const std::vector<std::vector<float>> &M) const {
    BigMatrix result = BigMatrix(n(), M[0].size());
    for (int i = 0; i < n(); i++)
        for (int j = 0; j < M[0].size(); j++)
            for (int k = 0; k < m(); k++)
                result.data[i][j] += this->data[i][k] * M[k][j];
    return result;
}


float dot(BigVector a, BigVector b) {
	float res = 0;
	for (int i = 0; i < a.data.size(); i++)
		res += a[i] * b[i];
	return res;
}

BigMatrix operator*(const std::vector<std::vector<float>> &M, const BigMatrix &B) {
    BigMatrix result = BigMatrix(M.size(), B.m());
    for (int i = 0; i < M.size(); i++)
        for (int j = 0; j < B.m(); j++)
            for (int k = 0; k < M[0].size(); k++)
                result.data[i][j] += M[i][k] * B.data[k][j];
    return result;
}

BigMatrix BigMatrix::inv() {
    if (n() != m())  throw std::invalid_argument("Matrix must be square");
    float d = det();
    if (d == 0)  throw std::invalid_argument("Matrix must be invertible");
    BigMatrix result = BigMatrix(n(), m());
    for (int i = 0; i < n(); i++)
        for (int j = 0; j < m(); j++)
            result.data[i][j] = submatrix(j, i).det() * ((i + j) % 2 == 0 ? 1 : -1);
    return result/d;
}
    BigMatrix BigMatrix::pow(int p) {
    if (n() != m())  throw std::invalid_argument("Matrix must be square");
    if (p < 0)  return ~(*this).pow(-p);
    if (p == 0)  return BigMatrix(n(), m());
    if (p == 1)  return *this;
    if (p == 2)  return *this * *this;
    if (p % 2 == 0)  return this->pow(p / 2) * this->pow(p / 2);
    return (*this) * this->pow(p / 2) * this->pow(p / 2);
}
    BigMatrix BigMatrix::GramSchmidtProcess() {
    BigMatrix result = BigMatrix(this->data);
    for (int i = 0; i < n(); i++) {
        for (int j = 0; j < m(); j++) {
            if (i == 0)  result.data[i][j] = this->data[i][j];
            else {
                float sum = 0;
                for (int k = 0; k < i; k++)  sum += result.data[k][j] * this->data[i][k];
                result.data[i][j] = this->data[i][j] - sum;
            }
        }
    }
    return result;
}
//    BigMatrix BigMatrix::operator*(const std::vector<float> &vec) {
//    return (*this)*BigMatrix({vec}).transpose();
//}
BigMatrix::BigMatrix(MATR$X && data) { data = std::move(data); }

std::pair<Complex, Complex> eigenvalues(mat2 m) {
	float tr = m[0][0] + m[1][1];
	float det = m[0][0] * m[1][1] - m[0][1] * m[1][0];
	float D = tr * tr - 4 * det;
	if (D < 0)  return {(tr + I*sqrt(-D)) / 2, (tr - I*sqrt(-D)) / 2};
	return {(tr + sqrt(D)) / 2, (tr - sqrt(D)) / 2};
}


std::pair<vec2, mat2> eigendecomposition(mat2 m) {
	vec2 lambda = eigenvaluesReal(m);
	if (isClose(m[0][1]*m[1][0], 0))
		lambda = vec2(m[0][0], m[1][1]);
	vec2 v1 = vec2(lambda.x - m[1][1], m[0][1]);
	vec2 v2 = vec2(m[1][0], lambda.y - m[0][0]);
	return std::make_pair(lambda, mat2(v1, v2));
}
