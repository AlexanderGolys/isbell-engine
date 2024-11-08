#include "mat.hpp"

#include <iostream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/transform.hpp>
#include <cmath>
#include <vector>
#include <set>
#include <array>

using namespace glm;
using std::vector, std::string, std::set, std::array, std::pair, std::exp, std::log, std::cos, std::sin, std::cosh, std::sinh, std::sqrt, std::pow, std::atan2, std::abs;


Complex::Complex()
{
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
	{
		return CP1(Complex(1, 0), true);
	}
	return CP1(z + vec2(f, 0));
}

CP1 CP1::operator-(float f)
{
	if (inf)
	{
		return CP1(Complex(1, 0), true);
	}
	return CP1(z - vec2(f, 0));
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




float norm2(vec2 v)
{
	return dot(v, v);
}

float norm2(vec3 v)
{
	return dot(v, v);
}

float norm2(vec4 v)
{
	return dot(v, v);
}


vec2 intersectLines(vec2 p1, vec2 p2, vec2 q1, vec2 q2)
{
	return (p1 * (q1.y - q2.y) - p2 * (q1.y - q2.y) - q1 * (p1.y - p2.y) + q2 * (p1.y - p2.y)) / ((p1.x - p2.x) * (q1.y - q2.y) - (p1.y - p2.y) * (q1.x - q2.x));
}


Complex intersectLines(Complex p1, Complex p2, Complex q1, Complex q2)
{
	return Complex(intersectLines(p1.z, p2.z, q1.z, q2.z));
}

mat3 scaleMatrix3(vec3 s)
{
    return mat3(s.x, 0, 0, 0, s.y, 0, 0, 0, s.z);
}

mat3 scaleMatrix3(float s)
{
	return scaleMatrix3(vec3(s, s, s));
}

mat3 changeOfBasis(vec3 target1, vec3 target2, vec3 target3)
{
    return inverse(mat3(target1, target2, target3));
}

mat3 changeOfBasis(vec3 source1, vec3 source2, vec3 source3, vec3 target1, vec3 target2, vec3 target3)
{
	return inverse(mat3(target1, target2, target3)) * mat3(source1, source2, source3);
}

mat3 rotationMatrix3(float angle)
{
    return mat3(cos(angle), -sin(angle), 0, sin(angle), cos(angle), 0, 0, 0, 1);
}

mat3 rotationMatrix3(vec3 axis, float angle)
{
    auto planeBasis = orthogonalComplementBasis(axis);
	mat3 change = changeOfBasis(planeBasis.first, planeBasis.second, axis);
	return inverse(change) * rotationMatrix3(angle) * change;
}

float frac(float x)
{
    return x - floor(x);
}

vec2 orthogonalComplement(vec2 v)
{
	return (vec2(v.y, -v.x))/norm(v);
}

mat3 GramSchmidtProcess(mat3 m) {
    vec3 u1 = normalise(m[0]);
    vec3 u2 = normalise(m[1] - u1 * dot(u1, m[1]));
    vec3 u3 = cross(u1, u2);
    return mat3(u1, u2, u3);
}

std::pair<vec3, vec3> orthogonalComplementBasis(vec3 v)
{
    vec3 b1 = norm(cross(v, e1)) > norm(cross(v, e3)) ? e1 : e3;
    vec3 b2 = cross(v, b1);
    mat3 frame = GramSchmidtProcess(mat3(b1, b2, v));
    return std::make_pair(frame[0], frame[1]);

}

float pseudorandomizer(float x, float seed)
{
	return frac(sin(x + seed) * 43758.5453f + seed);
}



//template <RingConcept T, int n, int m>
//T Matrix<T, n, m>::mobius(T z) requires DivisionRing<T>
//{
//	if (n != 2 || m != 2) {
//		throw std::invalid_argument("Matrix must be 2x2");
//	}
//	return (this->coefs[0][0] * z + this->coefs[0][1]) / (this->coefs[1][0] * z + this->coefs[1][1]);
//}

//template <RingConcept T, int n, int m>
//std::function<T(T)> Matrix<T, n, m>::mobius() requires DivisionRing<T>
//{
//	if (n != 2 || m != 2) {
//		throw std::invalid_argument("Matrix must be 2x2");
//	}
//	return [this](T z) {return  (this->coefs[0][0] * z + this->coefs[0][1]) / (this->coefs[1][0] * z + this->coefs[1][1]); };
//}

