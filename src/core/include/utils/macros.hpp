#pragma once

#include <filesystem>
#include <variant>
#include <string>
#include <optional>
#include <functional>
#include <iostream>
#include <format>
#include <map>
#include <memory>
#include <initializer_list>
#include <string_view>
#include <exception>
#include <expected>

#include "glm/glm.hpp"


using glm::vec2, glm::vec3, glm::vec4, glm::dvec2, glm::dvec3, glm::dvec4;
using glm::mat2, glm::mat3, glm::mat4, glm::mat2x3, glm::mat2x4, glm::mat3x2, glm::mat3x4, glm::mat4x2, glm::mat4x3;
using glm::ivec2, glm::ivec3, glm::ivec4;
using std::vector, std::variant, std::optional, std::array, std::unordered_map, std::pair, std::initializer_list;
using std::shared_ptr, std::unique_ptr, std::make_unique, std::make_shared, std::weak_ptr;
using std::string, std::format, std::to_string, std::printf, std::size_t;
using std::expected, std::unexpected, std::bad_expected_access;
using std::exp, std::log, std::cos, std::sin, std::cosh, std::sinh, std::sqrt, std::pow, std::atan2, std::abs;
namespace filesystem = std::filesystem;

constexpr float PI = 3.14159265359f;
constexpr float TAU = 6.28318530718f;

using array_len = size_t;
using array_index = size_t;
using byte_size = size_t;
using vs_dim = unsigned short;
using raw_data_ptr = const void*;
using uint = unsigned int;
using uchar = unsigned char;
using string_cr = const string&;

template <typename T>
using sptr = shared_ptr<T>;

template <typename T>
using uptr = unique_ptr<T>;

template <typename T>
using data_ptr = const T*;

template <typename T>
using data_ptr_mut = T*;

// template <typename A, typename B>
// using dict = std::map<A, B>;

#define INT(A) static_cast<int>(A)
#define UNDEFINED std::nullopt

using PolyGroupID = variant<int, string>;
#define Mat2C Matrix<Complex, 2>
#define isClose nearlyEqual

#define MAYBE(A) optional<A>
#define maybe(A) optional<A>

#define HOM(B, A) std::function<A(B)>
#define BIHOM(A, B, C) std::function<C(A,B)>
#define TRIHOM(A, B, C, D) std::function<D(A,B,C)>
#define QUADHOM(A, B, C, D, E) std::function<E(A,B,C,D)>
#define END(A) HOM(A, A)

#define Fooo END(float)
#define Foo12 HOM(float, vec2)
#define Foo13 HOM(float, vec3)
#define Foo33 HOM(vec3, vec3)
#define Foo31 HOM(vec3, float)
#define Foo21 HOM(vec2, float)
#define Foo22 HOM(vec2, vec2)
#define Foo113 BIHOM(float, float, vec3)
#define Foo112 BIHOM(float, float, vec2)
#define Foo111 BIHOM(float, float, float)

#define Foo1111 TRIHOM(float, float, float, float)
#define Foo1113 TRIHOM(float, float, float, vec3)
#define Foo313 BIHOM(vec3, float, vec3)
#define Foo311 BIHOM(vec3, float, float)


#define makeUptr std::make_unique
#define makeSptr std::make_shared

#define dict(K, V) std::map<K, V>

#define INPLACE_HOM(A) std::function<void(A)>
#define INPLACE_BIHOM(A, B) std::function<void(A, B)>
#define INPLACE_TRIHOM(A, B, C) std::function<void(A, B, C)>
#define IP_HOM(A) std::function<void(A)>
#define IP_BIHOM(A, B) std::function<void(A, B)>
#define IP_TRIHOM(A, B, C) std::function<void(A, B, C)>
#define betterInFamily(A) HOM(float, A)
#define procrastinateIn(A) HOM(A, void)
#define alboLeniwyAlboCwaniak HOM(void, void)
#define pencilCurv betterInFamily(SmoothParametricCurve)
#define pencilSurf betterInFamily(SmoothParametricSurface)
#define Foo3Foo33  HOM(vec3, mat3)
#define Foo2Foo22  HOM(glm::vec2, glm::mat2)
#define COPROD(X, Y) std::variant<X, Y>
#define PROD(X, Y) std::pair<X, Y>

#define BUFF2  vector<vec2>
#define BUFF3  vector<vec3>
#define BUFF4  vector<vec4>
#define IBUFF3 vector<ivec3>

#define pack(cap, F, vec, type) [cap](type a, type b){return F(vec(a, b));}
#define unpack(cup, F, vec) [cup](vec v){return F(v[0], v[1]);}
#define unpackPair(F, vec, ...) [F](vec v, ...){return F(v.first, v.second, ...);}
#define unpackPairWeird(_F, F, vec, ...) [F](vec v, ...){return F(v.first, v.second, ...);}
#define unpackTriple(F, vec, ...) [F](vec v, ...){return F(v[0], v[1], v[2], ...);}
#define unpackTripleWeird(_F, F, vec, ...) [_F](vec v, ...){return F(v[0], v[1], v[2], ...);}
#define unpackQuad(F, vec, ...) [F](vec v, ...){return F(v[0], v[1], v[2], v[3], ...);}

#define THROW(ErrType, ...) throw ErrType(__VA_ARGS__ __VA_OPT__(,) __FILE__, __LINE__)
#define THROW_IF(cond, ErrType, ...) if (cond) THROW(ErrType, __VA_ARGS__)

#define CHECK_OUT_OF_BOUNDS(idx, size) THROW_IF(idx < 0 or idx >= size, IndexOutOfBounds, idx, size)
#define CHECK_OUT_OF_BOUNDS_NAME(idx, size, dsname) THROW_IF(idx < 0 or idx >= size, IndexOutOfBounds, idx, size, dsname)

#define DIRTY_FLAG \
	private: \
	bool dirty = true; \
	public: \
	void markDirty() { dirty = true; } \
	void markClean() { dirty = false; } \
	bool isDirty() const { return dirty; }

#define GETTER(type, name) type get_##name() const { return name; }
#define SETTER(type, name) void set_##name(const type& value) { this->name = value; }

#define PROPERTY(type, name) \
private: \
	type name; \
public: \
	GETTER(type, name) \
	SETTER(type, name)

#define CONST_PROPERTY(type, name) \
private: \
	type name; \
public: \
	GETTER(type, name) \

#define SETTER_DIRTY(type, name) void set_##name(const type& value) { name = value; markDirty(); }

#define PROPERTY_DIRTY(type, name)\
private: \
	type name; \
public: \
	GETTER(type, name) \
	SETTER_DIRTY(type, name)

#ifdef _WIN32
    #define IS_WINDOWS true
#else
    #define IS_WINDOWS false
#endif
