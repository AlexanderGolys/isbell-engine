#pragma once

#include <exception>
#include <map>
#include <set>
#include <string>
#include "macros.hpp"
#include <glm/glm.hpp>


inline std::string vecToString(glm::vec2 v) {
	return "(" + std::to_string(v.x) + ", " + std::to_string(v.y) + ")"; }

inline std::string vecToString(glm::vec3 v) {
	return "(" + std::to_string(v.x) + ", " + std::to_string(v.y) + ", " + std::to_string(v.z) + ")"; }

inline std::string vecToString(glm::vec4 v) {
	return  "(" + std::to_string(v.x) + ", " + std::to_string(v.y) + ", " + std::to_string(v.z) + ", " + std::to_string(v.w) + ")"; }

inline std::string vecToString(glm::ivec2 v) {
	return "(" + std::to_string(v.x) + ", " + std::to_string(v.y) + ")"; }

inline std::string vecToString(glm::ivec3 v) {
	return "(" + std::to_string(v.x) + ", " + std::to_string(v.y) + ", " + std::to_string(v.z) + ")"; }

inline std::string vecToString(glm::ivec4 v) {
	return  "(" + std::to_string(v.x) + ", " + std::to_string(v.y) + ", " + std::to_string(v.z) + ", " + std::to_string(v.w) + ")"; }


inline std::string plural (std::string word) {
    if (word.size() == 0)
      return word;
    bool caps = word.back() >= 'A' && word.back() <= 'Z' && word.size() > 1;
    bool es = word.back() == 's'||
              word.back() == 'x' ||
              word.back() == 'z' ||
              (word.size() > 1 && word[word.size()-2] == 'c' && word.back() == 'h') ||
              (word.size() > 1 && word[word.size()-2] == 's' && word.back() == 'h') ||
              word.back() == 'S'||
              word.back() == 'X' ||
              word.back() == 'Z' ||
              (word.size() > 1 && word[word.size()-2] == 'C' && word.back() == 'H') ||
              (word.size() > 1 && word[word.size()-2] == 'S' && word.back() == 'H');
    if (es && !caps)
      return word + "es";
    if (es)
      return word + "ES";
    if (caps)
      return word + "S";
    return word + "s";
}

class NotImplementedError : public std::exception {
    std::string msg_;
public:
    explicit NotImplementedError(const std::string& notImplementedMethodName, const std::string& lackingType="Method")
        : msg_(lackingType + " " + notImplementedMethodName + " is not implemented yet.") {}

    const char* what() const noexcept override {
        return msg_.c_str();
    }
};

class IndexOutOfBounds : public std::exception {
	std::string msg_;
public:
	IndexOutOfBounds(int index, int size, const std::string &indexName="i")
		: msg_("Index " + indexName + " is out of bounds [" + std::to_string(index) + "/" + std::to_string(size) + "].") {}
	IndexOutOfBounds(const std::string &index, const std::string &size, const std::string &indexName="i")
		: msg_("Index " + indexName + " is out of bounds [" + index + "/" + size + "].") {}
	IndexOutOfBounds(const ivec2 &index, const ivec2 &size, const std::string &indexName="i")
		: msg_("Index " + indexName + " is out of bounds [" + vecToString(index) + "/" + vecToString(size) + "].") {}
	IndexOutOfBounds(const ivec3 &index, const ivec3 &size, const std::string &indexName="i")
	: msg_("Index " + indexName + " is out of bounds [" + vecToString(index) + "/" + vecToString(size) + "].") {}
	IndexOutOfBounds(const ivec4 &index, const ivec4 &size, const std::string &indexName="i")
	: msg_("Index " + indexName + " is out of bounds [" + vecToString(index) + "/" + vecToString(size) + "].") {}

	const char* what() const noexcept override {
		return msg_.c_str();
	}
};

class NotImplementedMethodError : public NotImplementedError {
public:
  explicit NotImplementedMethodError(const std::string& methodName)
      :NotImplementedError(methodName, "Method") {}
};

class NotImplementedFunctionError : public NotImplementedError {
public:
  explicit NotImplementedFunctionError(const std::string& name)
      :NotImplementedError(name, "Function") {}
};

class NotImplementedVariantError : public NotImplementedError {
public:
  NotImplementedVariantError(const std::string& variant, const std::string& ofWhat)
      :NotImplementedError(variant + " of " + ofWhat, "Variant") {}
};

class UnknownVariantError : public std::exception {
  std::string msg_;
public:
  UnknownVariantError(const std::string& variant, const std::string& ofWhat)
      : msg_("Variant " + variant + " of  " + ofWhat + " is an unknown type in this context.") {}
  explicit UnknownVariantError(const std::string& msg)
      : msg_(msg) {}


  const char* what() const noexcept override {
    return msg_.c_str();
  }
};

class IllegalVariantError : public std::exception {
  std::string msg_;
public:
  IllegalVariantError(const std::string& variant, const std::string& ofWhat, const std::string& rejectingMethod)
      : msg_(plural(ofWhat) + " in variant " + variant + " are considered invalid by method " + rejectingMethod + ".") {}
  explicit IllegalVariantError(const std::string& msg)
      : msg_(msg) {}

  const char* what() const noexcept override {
    return msg_.c_str();
  }
};

class IllegalArgumentError : public std::exception {
	std::string msg_;
public:
	explicit IllegalArgumentError(const std::string& msg) : msg_(msg) {}
	const char* what() const noexcept override {
		return msg_.c_str();
	}
};

class ValueError : public std::exception {
	std::string msg_;
public:
	explicit ValueError(const std::string& msg) : msg_(msg) {}
	const char* what() const noexcept override {
		return msg_.c_str();
	}
};


class COLOR_PALETTE {
public:
    glm::vec4 mainColor;
    glm::vec4 second;
    glm::vec4 third;
    glm::vec4 accent;
    glm::vec4 accent2;

    COLOR_PALETTE(glm::vec4 mainColor, glm::vec4 second, glm::vec4 third, glm::vec4 accent, glm::vec4 accent2);
    COLOR_PALETTE(glm::vec3 mainColor, glm::vec3 second, glm::vec3 third, glm::vec3 accent, glm::vec3 accent2);
    COLOR_PALETTE(glm::ivec3 mainColor, glm::ivec3 second, glm::ivec3 third, glm::ivec3 accent, glm::ivec3 accent2);
    std::vector<glm::vec4> colors();
    glm::vec4 operator[] (int i);

};

class COLOR_PALETTE10 {
public:
    std::array<glm::vec4, 10> cls;

    explicit COLOR_PALETTE10(std::array<glm::vec4, 10> colors) : cls(colors) {}
    COLOR_PALETTE10(COLOR_PALETTE p1, COLOR_PALETTE p2) : cls({p1[0], p1[1], p1[2], p1[3], p1[4], p2[0], p2[1], p2[2], p2[3], p2[4]}) {}
    COLOR_PALETTE10(glm::ivec3 c1, glm::ivec3 c2, glm::ivec3 c3, glm::ivec3 c4, glm::ivec3 c5, glm::ivec3 c6, glm::ivec3 c7, glm::ivec3 c8, glm::ivec3 c9, glm::ivec3 c10);
    std::vector<glm::vec4> colors() const { return std::vector<glm::vec4>(cls.begin(), cls.end()); }
    glm::vec4 operator[] (int i) const { return cls[i]; }
};


namespace glm {
    inline vec3 xyz(const vec4& v) {
        return vec3(v.x, v.y, v.z);
    }
    inline vec3 yzw(const vec4& v) {
        return vec3(v.y, v.z, v.w);
    }
    inline vec2 xy(const vec4& v) {
        return vec2(v.x, v.y);
    }
    inline vec2 yz(const vec4& v) {
        return vec2(v.y, v.z);
    }
    inline vec2 zw(const vec4& v) {
        return vec2(v.z, v.w);
    }
    inline vec2 xy(const vec3& v) {
        return vec2(v.x, v.y);
    }
    inline vec2 yz(const vec3& v) {
        return vec2(v.y, v.z);
    }

    inline vec4 xyz1(const vec3& v) {
        return vec4(v, 1);
    }

	inline vec2 operator*(const vec2& v, int i) { return v*(1.f*i); }
	inline vec3 operator*(const vec3& v, int i) { return v*(1.f*i); }
	inline vec4 operator*(const vec4& v, int i) { return v*(1.f*i); }
	inline mat2 operator*(const mat2& v, int i) { return v*(1.f*i); }
	inline mat3 operator*(const mat3& v, int i) { return v*(1.f*i); }
	inline mat4 operator*(const mat4& v, int i) { return v*(1.f*i); }
	inline vec2 operator/(const vec2& v, int i) { return v/(1.f*i); }
	inline vec3 operator/(const vec3& v, int i) { return v/(1.f*i); }
	inline vec4 operator/(const vec4& v, int i) { return v/(1.f*i); }
	inline mat2 operator/(const mat2& v, int i) { return v/(1.f*i); }
	inline mat3 operator/(const mat3& v, int i) { return v/(1.f*i); }
	inline mat4 operator/(const mat4& v, int i) { return v/(1.f*i); }

}

    template<typename T>
    void printVector(std::vector<T> v, std::string title="vector")
    {
        std::cout << title << " [" << v.size() << "]: ";
        for (int i = 0; i < v.size(); i++)
        {
            std::cout << v[i] << ", ";
        }
        std::cout << std::endl;
    }


const PolyGroupID DEFAULT_POLY_GROUP_ID = PolyGroupID(0);
constexpr RP1 inf = std::nullopt;
constexpr RP1 unbounded = std::nullopt;
const glm::vec3 e1 = glm::vec3(1, 0, 0);
const glm::vec3 e2 = glm::vec3(0, 1, 0);
const glm::vec3 e3 = glm::vec3(0, 0, 1);
const glm::vec3 ORIGIN = glm::vec3(0, 0, 0);
const glm::vec2 PLANE_ORIGIN = glm::vec2(0, 0);
const PolyGroupID DFLT_CURV = PolyGroupID(420);


class ivec8 {
	glm::ivec4 a, b;
public:
	ivec8(glm::ivec4 a, glm::ivec4 b) : a(a), b(b) {}
	ivec8(int a, int b, int c, int d, int e, int f, int g, int h) : a(a, b, c, d), b(e, f, g, h) {}
	explicit ivec8(int i) : a(i), b(i) {}
	int operator[](int i) const { return i < 4 ? a[i] : b[i - 4]; }
	ivec8 operator+(int i) const { return ivec8(a + i, b + i); }
	ivec8 operator-(int i) const { return ivec8(a - i, b - i); }
	ivec8 operator*(int i) const { return ivec8(a * i, b * i); }
	ivec8 operator/(int i) const { return ivec8(a / i, b / i); }
	ivec8 operator+(ivec8 v) const { return ivec8(a + v.a, b + v.b); }
	ivec8 operator-(ivec8 v) const { return ivec8(a - v.a, b - v.b); }
	ivec8 operator*(ivec8 v) const { return ivec8(a * v.a, b * v.b); }
	ivec8 operator/(ivec8 v) const { return ivec8(a / v.a, b / v.b); }
	ivec8 operator-() const { return ivec8(-a, -b); }
	ivec8 operator%(ivec8 v) const { return ivec8(a % v.a, b % v.b); }
	static int size() { return 8; }
//	friend int dot(ivec8 a, ivec8 b) { return dot(1.F*a.a, b.a) + dot(a.b, b.b); }
	glm::ivec4 xyzw() const { return a; }
	glm::ivec4 stuv() const { return b; }
	glm::ivec3 xyz() const { return glm::ivec3(a); }
	glm::ivec3 yzw() const { return glm::ivec3(a.y, a.z, a.w); }
	glm::ivec3 stv() const { return glm::ivec3(b); }
	glm::ivec3 tvu() const { return glm::ivec3(b.y, b.z, b.w); }
	glm::ivec2 xy() const { return glm::ivec2(a); }
	glm::ivec2 yz() const { return glm::ivec2(a.y, a.z); }
	glm::ivec2 zw() const { return glm::ivec2(a.z, a.w); }
	glm::ivec2 st() const { return glm::ivec2(b); }
	glm::ivec2 tu() const { return glm::ivec2(b.y, b.z); }
	glm::ivec2 uv() const { return glm::ivec2(b.z, b.w); }
	int x() const { return a.x; }
	int y() const { return a.y; }
	int z() const { return a.z; }
	int w() const { return a.w; }
};

class ivec6 {
	glm::ivec3 a, b;
public:
	ivec6(glm::ivec3 a, glm::ivec3 b) : a(a), b(b) {}
	ivec6(int a, int b, int c, int d, int e, int f) : a(a, b, c), b(d, e, f) {}
	explicit ivec6(int i) : a(i), b(i) {}
	int operator[](int i) const { return i < 3 ? a[i] : b[i - 3]; }
	ivec6 operator+(int i) const { return ivec6(a + i, b + i); }
	ivec6 operator-(int i) const { return ivec6(a - i, b - i); }
	ivec6 operator*(int i) const { return ivec6(a * i, b * i); }
	ivec6 operator/(int i) const { return ivec6(a / i, b / i); }
	ivec6 operator+(ivec6 v) const { return ivec6(a + v.a, b + v.b); }
	ivec6 operator-(ivec6 v) const { return ivec6(a - v.a, b - v.b); }
	ivec6 operator*(ivec6 v) const { return ivec6(a * v.a, b * v.b); }
	ivec6 operator/(ivec6 v) const { return ivec6(a / v.a, b / v.b); }
	ivec6 operator-() const { return ivec6(-a, -b); }
	ivec6 operator%(ivec6 v) const { return ivec6(a % v.a, b % v.b); }
	static int size() { return 6; }
	bool operator==(ivec6 v) const { return a == v.a && b == v.b; }
//	friend int dot(ivec6 a, ivec6 b) { return dot(a.a, b.a) + dot(a.b, b.b); }
	glm::ivec3 xyz() const { return glm::ivec3(a); }
	glm::ivec3 stv() const { return glm::ivec3(b); }
	glm::ivec2 xy() const { return glm::ivec2(a); }
	glm::ivec2 yz() const { return glm::ivec2(a.y, a.z); }
	glm::ivec2 st() const { return glm::ivec2(b); }
	glm::ivec2 tu() const { return glm::ivec2(b.y, b.z); }
	int x() const { return a.x; }
	int y() const { return a.y; }
	int z() const { return a.z; }
	int s() const { return b.x; }
	int t() const { return b.y; }
	int u() const { return b.z; }

	std::vector<int> toVec() const { return {a.x, a.y, a.z, b.x, b.y, b.z}; }
};


template<typename T>
std::vector<T> flattened2DVector(std::vector<std::vector<T>> v)
{
	std::vector<T> res;
	for (auto& row : v)
		res.insert(res.end(), row.begin(), row.end());
	return res;
}


inline int flattened2DVectorIndex(int i, int j, glm::ivec2 size) {
	if (i < 0) return flattened2DVectorIndex(size.x + i, j, size);
	if (j < 0) return flattened2DVectorIndex(i, size.y + j, size);
	if (i >= size.x) throw IndexOutOfBounds(ivec2(i, j), size, "i");
	if (j >= size.y) throw IndexOutOfBounds(ivec2(i, j), size, "j");
	return i * size.y + j;
}


template<typename T>
T flattened2DVectorSample(std::vector<T> flattened, int i, int j, glm::ivec2 size) {
	return flattened[flattened2DVectorIndex(i, j, size)];
}



template<typename T>
std::vector<T> flattened3DVector(std::vector<std::vector<std::vector<T>>> v)
{
	std::vector<T> res;
	for (auto& row : v)
		for (auto& col : row)
			res.insert(res.end(), col.begin(), col.end());
	return res;
}


inline int flattened3DVectorIndex(int i, int j, int k,  glm::ivec3 size) {
	if (i < 0) return flattened3DVectorIndex(size.x + i, j, k, size);
	if (j < 0) return flattened3DVectorIndex(i, size.y + j, k, size);
	if (k < 0) return flattened3DVectorIndex(i, j, size.z + k, size);
	if (i >= size.x) throw IndexOutOfBounds(ivec3(i, j, k), size, "i");
	if (j >= size.y) throw IndexOutOfBounds(ivec3(i, j, k), size, "j");
	if (k >= size.z) throw IndexOutOfBounds(ivec3(i, j, k), size, "k");
	return i * size.y * size.z + j * size.z + k;
}


template<typename T>
T flattened3DVectorSample(std::vector<T> flattened, int i, int j, int k, glm::ivec3 size) {
	return flattened[flattened3DVectorIndex(i, j, k, size)];
}

inline int flattened4DVectorIndex(int i, int j, int k, int m,  glm::ivec4 size) {
	if (i < 0) return flattened4DVectorIndex(size.x + i, j, k,m, size);
	if (j < 0) return flattened4DVectorIndex(i, size.y + j, k,m, size);
	if (k < 0) return flattened4DVectorIndex(i, j, size.z + k,m, size);
	if (m < 0) return flattened4DVectorIndex(i, j, k, size.w + m, size);
	if (i >= size.x) throw IndexOutOfBounds(ivec4(i, j, k, m), size, "i");
	if (j >= size.y) throw IndexOutOfBounds(ivec4(i, j, k, m), size, "j");
	if (k >= size.z) throw IndexOutOfBounds(ivec4(i, j, k, m), size, "k");
	if (m >= size.w) throw IndexOutOfBounds(ivec4(i, j, k, m), size, "m");
	return i * size.y * size.z * size.w + j * size.z * size.w + k * size.w + m;
}


template<typename T>
T flattened4DVectorSample(std::vector<T> flattened, int i, int j, int k, int m, glm::ivec4 size) {
	return flattened[flattened4DVectorIndex(i, j, k, m, size)];
}

template<typename T>
bool contains(const std::vector<T> &v, T x) {
	return std::find(v.begin(), v.end(), x) != v.end();
}

template<typename T>
bool containsAll(const std::vector<T> &big, const std::vector<T> &subset) {
	for (auto& x : subset)
		if (!contains(big, x))
			return false;
	return true;
}

template<typename T, typename U>
std::vector<T> keys(const std::map<T, U> &m) {
	std::vector<T> res;
	for (auto& [k, v] : m)
		res.push_back(k);
	return res;
}

template<typename T, typename U>
std::set<U> values(const std::map<T, U> &m) {
	std::set<U> res;
	for (auto& [k, v] : m)
		res.insert(v);
	return res;
}

template<typename T, typename U>
std::vector<U> valuesVec(const std::map<T, U> &m) {
	std::vector<U> res;
	for (auto& [k, v] : m)
		res.push_back(v);
	return res;
}

template<typename T>
std::vector<T> setToVector(const std::set<T> &s) {
	std::vector<T> res;
	for (auto& x : s)
		res.push_back(x);
	return res;
}

template<typename T, int n>
std::array<T, n> vectorToArray(const std::vector<T> &s) {
	std::array<T, n> res;
	for (int i = 0; i < n; i++)
		res[i] = s[i];
	return res;
}

template<typename T, int n, int m=n>
std::array<std::array<T, n>, m> vecVecToArray(const std::vector<T> &s) {
	std::array<std::array<T, n>, m> res;
	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++)
			res[i][j] = s[i * n + j];
	return res;
}

template<typename T, int n, int m=n>
std::vector<std::vector<T>> arrayToVecVec(const std::array<std::array<T, n>, m> &s) {
	std::vector<std::vector<T>> res;
	res.reserve(m);
	for (int i = 0; i < m; i++)
		res.emplace_back(s[i].begin(), s[i].end());
	return res;
}

template<typename T, int n>
std::array<T, n> setToArray(const std::set<T> &s) {
	return vectorToArray<T, n>(setToVector(s));
}

inline float dot(float x, float y) { return x * y; }

template<typename  V>
float norm2(V v) { return dot(v, v); }

inline float norm2(mat2 m) { return m[0][0]*m[0][0] + m[0][1]*m[0][1] + m[1][0]*m[1][0] + m[1][1]*m[1][1]; }

template <typename  T>
float norm(T v) { return sqrt(norm2(v)); }

inline std::array<std::array<float, 2>, 2> mat2ToArray(mat2 m) {
	return vectorToArray<std::array<float, 2>, 2>({{m[0][0], m[0][1]}, {m[1][0], m[1][1]} }); }

inline std::array<std::array<float, 3>, 3> mat3ToArray(mat3 m) {
	return vectorToArray<std::array<float, 3>, 3>({{m[0][0], m[0][1], m[0][2]}, {m[1][0], m[1][1], m[1][2]}, {m[2][0], m[2][1], m[2][2]} }); }

//inline std::array<std::array<float, 4>, 4> mat4ToArray(mat4 m) {
//	return {{{m[0][0], m[0][1], m[0][2], m[0][3]}, {m[1][0], m[1][1], m[1][2], m[1][3]}, {m[2][0], m[2][1], m[2][2], m[2][3]}, {m[3][0], m[3][1], m[3][2], m[3][3]} }};}


inline MATR$X matToVecVec(const mat2 &m) { return { {m[0][0], m[0][1]}, {m[1][0], m[1][1]} }; }
inline MATR$X matToVecVec(const mat3 &m) { return { {m[0][0], m[0][1], m[0][2]}, {m[1][0], m[1][1], m[1][2]}, {m[2][0], m[2][1], m[2][2]} }; }
//inline MATR$X matToVecVec(const mat4 &m) { return arrayToVecVec(mat4ToArray(m)); }

template<typename T>
vector<T> flatten  (vector<vector<T>> v) {
	vector<T> res;
	for (auto& row : v)
		res.insert(res.end(), row.begin(), row.end());
	return res;
}

template<typename T>
vector<vector<T>> grid (T x0, T d1, T d2, int n1, int n2) {
	vector<vector<T>> res;
	for (int i = 0; i < n1; i++) {
		vector<T> row;
		for (int j = 0; j < n2; j++)
			row.push_back(x0 + d1 * i + d2 * j);
		res.push_back(row);
	}
	return res;
}

template<typename T>
vector<vector<vector<T>>> grid (T x0, T d1, T d2, T d3, int n1, int n2, int n3) {
	vector<vector<vector<T>>> res;
	for (int i = 0; i < n1; i++) {
		vector<vector<T>> row;
		for (int j = 0; j < n2; j++) {
			vector<T> col;
			for (int k = 0; k < n3; k++)
				col.push_back(x0 + d1 * i + d2 * j + d3 * k);
			row.push_back(col);
		}
		res.push_back(row);
	}
	return res;
}

inline bool contains(int i, ivec3 x) { return i == x.x || i == x.y || i == x.z; }
inline bool contains(int i, ivec2 x) { return i == x.x || i == x.y; }
inline bool contains(int i, ivec4 x) { return i == x.x || i == x.y || i == x.z || i == x.w; }
inline bool contains(ivec3 x, int i) { return contains(i, x);}
inline bool contains(ivec2 x, int i) { return contains(i, x);}
inline bool contains(ivec4 x, int i) { return contains(i, x);}

inline int setMinus(ivec3 x, ivec2 y) {
	if (x.x == y.x) return x.y == y.y ? x.z : x.y;
	if (x.x == y.y) return x.z;
	return x.x;
}

template<typename T>
std::vector<T> setMinus(std::vector<T> x, T y) {
	std::vector<T> res;
	res.reserve(x.size()-1);
	for (auto& z : x)
		if (z != y)
			res.push_back(z);
	return res;
}

template<typename T>
std::vector<T> setMinus(std::vector<T> x, std::vector<T> y) {
	std::vector<T> res;
	for (auto& z : x)
		if (!contains<T>(y, z))
			res.push_back(z);
	return res;
}

inline int setMinus(ivec4 x, ivec3 y) {
	if (x.x == y.x) return x.y == y.y ? x.z == y.z ? x.w : x.z : x.y;
	if (x.x == y.y) return x.z == y.z ? x.w : x.z;
	if (x.x == y.z) return x.w;
	return x.x;
}

inline int setMinus(ivec2 x, int y) {
	return x.x == y ? x.y : x.x;
}

inline ivec2 setMinus(ivec3 x, int y) {
	if (x.x == y) return ivec2(x.y, x.z);
	if (x.y == y) return ivec2(x.x, x.z);
	return ivec2(x.x, x.y);
}

inline ivec3 setMinus(ivec4 x, int y) {
	if (x.x == y) return ivec3(x.y, x.z, x.w);
	if (x.y == y) return ivec3(x.x, x.z, x.w);
	if (x.z == y) return ivec3(x.x, x.y, x.w);
	return ivec3(x.x, x.y, x.z);
}

inline ivec2 setMinus(ivec4 x, ivec2 y) {
	if (x.x == y.x) return x.y == y.y ? ivec2(x.z, x.w) : ivec2(x.y, x.z);
	if (x.x == y.y) return ivec2(x.y, x.w);
	return ivec2(x.x, x.y);
}

inline std::string hash_ivec3(ivec3 v) {
	return std::to_string(v.x) + "--" + std::to_string(v.y) + "--" + std::to_string(v.z);
}
