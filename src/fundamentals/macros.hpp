#pragma once
#include <glm/glm.hpp>

#include <variant>
#include <string>
#include <optional>
#include <functional>
#include <ranges>


using std::vector, glm::vec2, glm::vec3, glm::vec4, glm::mat3, glm::mat4, glm::mat2, glm::mat2x3, glm::ivec2, glm::ivec3;
using std::variant, std::string, std::optional;

#define PI 3.14159265359f
#define TAU 6.28318530718f
#define RP1 optional<float>
#define maybeFloat std::optional<float>
#define PolyGroupID std::variant<int, std::string>
#define Mat2C Matrix<Complex, 2>
#define INT(x) static_cast<int>(x)
#define End1P std::function<SpaceEndomorphism(float)>
#define End2P std::function<SpaceEndomorphism(float, float)>
#define maybeMaterial std::optional<MaterialPhongConstColor>
#define Mob Matrix<Complex, 2>
#define isClose nearlyEqual

#define GEN_VEC(R) GenericTensor<R, R>
#define GEN_MAT(R) GenericTensor<R, GEN_VEC(R)>

namespace std {
#define HOM(B,A) std::function<A(B)>
#define BIHOM(A, B, C) std::function<C(A,B)>

#define FUNC(B, ...) std::function<B(...)>
#define END(A) HOM(A, A)
#define HOM$$(B,A) std::vector<HOM(B, A)>
#define FUNC$$(B, A) std::vector<FUNC(B, A)>
#define END$$(A) std::vector<END(A)>

#define endC END(Complex)
#define Fooo END(float)
#define Foo12 HOM(float, vec2)
#define Foo13  HOM(float, vec3)
#define Foo33  HOM(vec3, vec3)
#define Foo31  HOM(vec3, float)
#define Foo21  HOM(vec2, float)
#define Foo22  HOM(vec2, vec2)
#define Foo32  HOM(vec3, vec2)
#define Foo23  HOM(vec2, vec3)
#define Foo113 std::function<vec3(float, float)>
#define Foo112 std::function<vec2(float, float)>
#define Foo111 std::function<float(float, float)>
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
#define IBUFF3 vector<glm::ivec3>





#define timePassesDoShit2(A, B) std::function<void(A, B)>
#define C0(A)         HOM(A, float)
#define Cinf(A)       RealFunctionR3
#define O_(A)         HOM(A, float)
#define burdl(A)      HOM(float, A)
#define curvy(A)      HOM(float, A)
#define path(A)       HOM(float, A)
#define TRG           const IndexedTriangle&
#define triple_O(A)   HOM(A, vec3)
#define O_x3(A)       HOM(A, vec3)
#define H0(A, F)      HOM(A, F)


#define pack(cap, F, vec, type) [cap](type a, type b){return F(vec(a, b));}
#define unpack(cup, F, vec) [cup](vec v){return F(v[0], v[1]);}
#define unpackPair(F, vec, ...) [F](vec v, ...){return F(v.first, v.second, ...);}
#define unpackPairWeird(_F, F, vec, ...) [F](vec v, ...){return F(v.first, v.second, ...);}
#define unpackTriple(F, vec, ...) [F](vec v, ...){return F(v[0], v[1], v[2], ...);}
#define unpackTripleWeird(_F, F, vec, ...) [_F](vec v, ...){return F(v[0], v[1], v[2], ...);}

#define unpackQuad(F, vec, ...) [F](vec v, ...){return F(v[0], v[1], v[2], v[3], ...);}
#define curry(F, arg, ...) [F](arg a, ...){return F(a, ...);}
#define uncarry(F, arg, arg2, ...) [F](HOM(arg, HOM(arg2, ...)) f, arg a, ...){return f(a)(...);}
}

// #define Foo(A)_Foo33 Foo33([](float){return hom(A, glm::mat3)(t);};)
#define vec420 std::vector<float>
#define vec69 std::vector<float>
#define vec2137 std::vector<float>

#define V3CTOR$(A) std::vector<std::vector<A>>
#define MATR$X V3CTOR$(float)

template <typename A, typename B, typename C>
std::function<A(std::function<B(C)>)> curr(std::function<A(B, C)> f) {
	return [f](std::function<B(C)> g) { return [f, g](C c) { return f(g(c), c); }; };
}



#define RealFunctionRigid    HOM(TRG, float)
#define RigidBundleV3 triple_O(TRG)
#define patrickJaneBundle    HOM(float, RigidBundleV3)
#define patrickJaneCaseRigid HOM(float, RealFunctionRigid)
// #define Foo1FooRigidFoo33 hom(float, Foo(RigidBody)_Foo33 )

// #define F[vec v] [F](vec v){return F(v[0], v[1]);}
#define BigMatrix_t HOM(float, BigMatrix)
#define vec69_t HOM(float, vec69)
#define rigidMotion std::pair<vec3, vec3>

#define BigMatrix_hm std::optional<BigMatrix>
#define vec69_hm std::optional<vec69>
#define mat3_hm std::optional<glm::mat3>
#define float_hm std::optional<float>
#define R3_hm std::optional<R3>
#define HM_NO std::nullopt

namespace glm{
#define BLACK glm::vec4(0, 0, 0, 1)
#define WHITE glm::vec4(1, 1, 1, 1)
#define RED glm::vec4(1, 0, 0, 1)
#define GREEN glm::vec4(0, 1, 0, 1)
#define BLUE glm::vec4(0, 0, 1, 1)
#define YELLOW glm::vec4(1, 1, 0, 1)
#define CYAN glm::vec4(0, 1, 1, 1)
#define MAGENTA glm::vec4(1, 0, 1, 1)
#define ORANGE glm::vec4(1, .5, 0, 1)
#define PINK glm::vec4(1, .5, .5, 0)
#define LIGHT_BLUE glm::vec4(.5, .5, 1, 1)
#define LIGHT_GREEN glm::vec4(.5, 1, .5, 1)
#define LIGHT_RED glm::vec4(1, .5, .5, 1)
#define LIGHT_YELLOW glm::vec4(1, 1, .5, 1)
#define LIGHT_CYAN glm::vec4(.5, 1, 1, 1)
#define LIGHT_MAGENTA glm::vec4(1, .5, 1, 1)
#define DARK_PURPLE glm::vec4(54.f/255, 31.f/255, 39.f/255, 0)
#define PURPLE glm::vec4(82.f/255, 25.f/255, 69.f/255, 1)
#define BROWN_SUGAR glm::vec4(169.f/255, 113.f/255, 75.f/255, 1)
#define POWDER_PINK glm::vec4(255.f/255, 179.f/255, 198.f/255, 1)
#define DARKDARKBLUE glm::vec4(0, 0.016, .045, 1)
// #define com glm::vec3


#define PALETTE1 COLOR_PALETTE(glm::ivec3(255, 166, 158), glm::ivec3(170, 68, 101), glm::ivec3(70, 34, 85), glm::ivec3(147, 225, 216), glm::ivec3(221, 255, 247))
#define PALETTE2 COLOR_PALETTE(glm::ivec3(168, 32, 26), glm::ivec3(60, 107, 45), glm::ivec3(215, 68, 56), glm::ivec3(204, 88, 3), glm::ivec3(2, 2, 2))
								// red red reddy
#define PALETTE3 COLOR_PALETTE(glm::ivec3(144, 156, 194), glm::ivec3(244, 157, 55), glm::ivec3(232, 72, 85), glm::ivec3(1, 2, 2), glm::ivec3(34, 34, 34))
                               // light powder blue			light intense orange           aggressive redpink        calm darker gray      dark purple, calm

#define INTENSE_DARKER_RED_PALLETTE COLOR_PALETTE(glm::ivec3(173, 40, 49), glm::ivec3(128, 14, 19), glm::ivec3(100, 13, 20), glm::ivec3(56, 4, 14), glm::ivec3(37, 9, 2))
#define BLUE_PALLETTE COLOR_PALETTE10(glm::ivec3(226, 239, 246), glm::ivec3(183, 215, 234), glm::ivec3(141, 190, 220), glm::ivec3(100, 164, 206), glm::ivec3(72, 147, 198), glm::ivec3(45, 130, 189), glm::ivec3(36, 118, 177), glm::ivec3(25, 101, 160), glm::ivec3(15, 85, 143), glm::ivec3(0, 58, 112))
#define GREEN_PALLETTE COLOR_PALETTE10(glm::ivec3(219, 234, 215), glm::ivec3(195, 220, 188), glm::ivec3(171, 206, 161), glm::ivec3(148, 193, 134), glm::ivec3(124, 179, 107), glm::ivec3(102, 162, 83), glm::ivec3(85, 135, 69), glm::ivec3(68, 108, 55), glm::ivec3(51, 81, 42), glm::ivec3(32, 51, 26))

#define REDPINK_PALLETTE COLOR_PALETTE10(glm::ivec3(243, 219, 234), glm::ivec3(236, 197, 221), glm::ivec3(204, 154, 181), glm::ivec3(179, 116, 149), glm::ivec3(153, 77, 116), glm::ivec3(128, 39, 84), glm::ivec3(102, 0, 51), glm::ivec3(85, 6, 45), glm::ivec3(67, 12, 39), glm::ivec3(52, 9, 30))

#define R3 glm::vec3
#define R2 glm::vec2
#define R4 glm::vec4

#define R3_t HOM(float, glm::vec3)

}
