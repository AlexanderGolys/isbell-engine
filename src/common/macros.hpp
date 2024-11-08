#pragma once
#include <variant>
#include <glm/glm.hpp>

#define PI 3.14159265359f
#define TAU 6.28318530718f
#define RP1 std::optional<float>
#define maybeFloat std::optional<float>
#define PolyGroupID std::variant<int, std::string>
#define Mat2C Matrix<Complex, 2>

const PolyGroupID DEFAULT_POLY_GROUP_ID = PolyGroupID(0);
constexpr RP1 inf = std::nullopt;
constexpr RP1 unbounded = std::nullopt;
const glm::vec3 e1 = glm::vec3(1, 0, 0);
const glm::vec3 e2 = glm::vec3(0, 1, 0);
const glm::vec3 e3 = glm::vec3(0, 0, 1);
const glm::vec3 ORIGIN = glm::vec3(0, 0, 0);
const glm::vec2 PLANE_ORIGIN = glm::vec2(0, 0);


#define BLACK glm::vec4(0, 0, 0, 1)
#define WHITE glm::vec4(1, 1, 1, 1)
#define RED glm::vec4(1, 0, 0, 1)
#define GREEN glm::vec4(0, 1, 0, 1)
#define BLUE glm::vec4(0, 0, 1, 1)
#define YELLOW glm::vec4(1, 1, 0, 1)
#define CYAN glm::vec4(0, 1, 1, 1)
#define MAGENTA glm::vec4(1, 0, 1, 1)
#define ORANGE glm::vec4(1, .5, 0, 1)
#define PINK glm::vec4(1, .5, .5, 1)
#define LIGHT_BLUE glm::vec4(.5, .5, 1, 1)
#define LIGHT_GREEN glm::vec4(.5, 1, .5, 1)
#define LIGHT_RED glm::vec4(1, .5, .5, 1)
#define LIGHT_YELLOW glm::vec4(1, 1, .5, 1)
#define LIGHT_CYAN glm::vec4(.5, 1, 1, 1)
#define LIGHT_MAGENTA glm::vec4(1, .5, 1, 1)
#define DARK_PURPLE glm::vec4(54.f/255, 31.f/255, 39.f/255, 1)
#define PURPLE glm::vec4(82.f/255, 25.f/255, 69.f/255, 1)
#define BROWN_SUGAR glm::vec4(169.f/255, 113.f/255, 75.f/255, 1)
#define POWDER_PINK glm::vec4(255.f/255, 179.f/255, 198.f/255, 1)
#define DARKDARKBLUE glm::vec4(0, 0.016, .045, 1)

#define PALETTE1 COLOR_PALETTE(glm::ivec3(255, 166, 158), glm::ivec3(170, 68, 101), glm::ivec3(70, 34, 85), glm::ivec3(147, 225, 216), glm::ivec3(221, 255, 247))
#define PALETTE2 COLOR_PALETTE(glm::ivec3(168, 32, 26), glm::ivec3(60, 107, 45), glm::ivec3(215, 68, 56), glm::ivec3(204, 88, 3), glm::ivec3(2, 2, 2))
								// red red reddy
#define PALETTE3 COLOR_PALETTE(glm::ivec3(144, 156, 194), glm::ivec3(244, 157, 55), glm::ivec3(232, 72, 85), glm::ivec3(1, 2, 2), glm::ivec3(34, 34, 34))
                               // light powder blue			light intense orange           aggressive redpink        calm darker gray      dark purple, calm

#define INTENSE_DARKER_RED_PALLETTE COLOR_PALETTE(glm::ivec3(173, 40, 49), glm::ivec3(128, 14, 19), glm::ivec3(100, 13, 20), glm::ivec3(56, 4, 14), glm::ivec3(37, 9, 2))
#define BLUE_PALLETTE COLOR_PALETTE10(glm::ivec3(226, 239, 246), glm::ivec3(183, 215, 234), glm::ivec3(141, 190, 220), glm::ivec3(100, 164, 206), glm::ivec3(72, 147, 198), glm::ivec3(45, 130, 189), glm::ivec3(36, 118, 177), glm::ivec3(25, 101, 160), glm::ivec3(15, 85, 143), glm::ivec3(0, 58, 112))
#define GREEN_PALLETTE COLOR_PALETTE10(glm::ivec3(219, 234, 215), glm::ivec3(195, 220, 188), glm::ivec3(171, 206, 161), glm::ivec3(148, 193, 134), glm::ivec3(124, 179, 107), glm::ivec3(102, 162, 83), glm::ivec3(85, 135, 69), glm::ivec3(68, 108, 55), glm::ivec3(51, 81, 42), glm::ivec3(32, 51, 26))

#define REDPINK_PALLETTE COLOR_PALETTE10(glm::ivec3(243, 219, 234), glm::ivec3(236, 197, 221), glm::ivec3(204, 154, 181), glm::ivec3(179, 116, 149), glm::ivec3(153, 77, 116), glm::ivec3(128, 39, 84), glm::ivec3(102, 0, 51), glm::ivec3(85, 6, 45), glm::ivec3(67, 12, 39), glm::ivec3(52, 9, 30))

