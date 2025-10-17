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
#include <filesystem>
#include <string>
#include <string_view>
#include <exception>
#include <expected>

#include "colors.hpp"


/**
* @note: the only commonly used std keywo♥rds that require namespace prefix are std::map and std::function
*/

using glm::vec2, glm::vec3, glm::vec4;
using glm::mat2, glm::mat3, glm::mat4, glm::mat2x3;
using glm::ivec2, glm::ivec3, glm::ivec4;
using std::vector, std::variant, std::optional, std::array, std::unordered_map, std::pair;
using std::shared_ptr, std::unique_ptr, std::make_unique, std::make_shared, std::weak_ptr;
using std::string, std::endl, std::format, std::to_string, std::cout, std::printf;
using std::expected, std::unexpected, std::bad_expected_access;
using std::exp, std::log, std::cos, std::sin, std::cosh, std::sinh, std::sqrt, std::pow, std::atan2, std::abs;
namespace filesystem = std::filesystem;

constexpr float PI = 3.14159265359f;
constexpr float TAU = 6.28318530718f;

// template <typename A, typename B>
// using dict = std::map<A, B>;

#define RP1 optional<float>
#define INT(A) static_cast<int>(A)
#define UNDEFINED std::nullopt

#define PolyGroupID variant<int, string>
#define Mat2C Matrix<Complex, 2>
#define maybeMaterial optional<MaterialPhongConstColor>
#define Mob Matrix<Complex>
#define isClose nearlyEqual

#define GEN_VEC(R) GenericTensor<R, R>
#define GEN_MAT(R) GenericTensor<R, GEN_VEC(R)>

// namespace std {

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
#define Foo32 HOM(vec3, vec2)
#define Foo23 HOM(vec2, vec3)
#define Foo113 BIHOM(float, float, vec3)
#define Foo112 BIHOM(float, float, vec2)
#define Foo111 BIHOM(float, float, float)

#define Foo1111 TRIHOM(float, float, float, float)
#define Foo1113 TRIHOM(float, float, float, vec3)
#define Foo313 BIHOM(vec3, float, vec3)
#define Foo311 BIHOM(vec3, float, float)

#define upt std::unique_ptr
#define spt std::shared_ptr
#define wpt std::weak_ptr
#define mupt std::make_unique
#define mspt std::make_shared

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



#define BUFF2  std::vector<vec2>
#define BUFF3  std::vector<vec3>
#define BUFF4  std::vector<vec4>
#define IBUFF3 std::vector<glm::ivec3>







#define pack(cap, F, vec, type) [cap](type a, type b){return F(vec(a, b));}
#define unpack(cup, F, vec) [cup](vec v){return F(v[0], v[1]);}
#define unpackPair(F, vec, ...) [F](vec v, ...){return F(v.first, v.second, ...);}
#define unpackPairWeird(_F, F, vec, ...) [F](vec v, ...){return F(v.first, v.second, ...);}
#define unpackTriple(F, vec, ...) [F](vec v, ...){return F(v[0], v[1], v[2], ...);}
#define unpackTripleWeird(_F, F, vec, ...) [_F](vec v, ...){return F(v[0], v[1], v[2], ...);}

#define unpackQuad(F, vec, ...) [F](vec v, ...){return F(v[0], v[1], v[2], v[3], ...);}
// #define curry(F, arg, ...) [F](arg a, ...){return F(a, ...);}
// #define uncarry(F, arg, arg2, ...) [F](HOM(arg, HOM(arg2, ...)) _f, arg a, ...){return _f(a)(...);}
// }

// #define Foo(A)_Foo33 Foo33([](float){return hom(A, glm::mat3)(v);};)
#define vec420 std::vector<float>
#define vec69 std::vector<float>
#define vec2137 std::vector<float>

#define V3CTOR$(A) std::vector<std::vector<A>>
#define MATR$X V3CTOR$(float)

template <typename A, typename B, typename C>
std::function<A(std::function<B(C)>)> curr(std::function<A(B, C)> f) {
	return [f](std::function<B(C)> g) { return [f, g](C c) { return f(g(c), c); }; };
}





#define float_hm optional<float>




// namespace color
// {


// }

#define THROW(ErrType, ...) throw ErrType(__VA_ARGS__ __VA_OPT__(,) __FILE__, __LINE__)
#define THROW_IF(cond, ErrType, ...) if (cond) throw ErrType(__VA_ARGS__ __VA_OPT__(,) __FILE__, __LINE__)


#ifdef _WIN32
    #define IS_WINDOWS true
#else
    #define IS_WINDOWS false
#endif
