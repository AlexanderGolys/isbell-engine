
#pragma once
#include "exceptions.hpp"

// ReSharper disable CppNonExplicitConvertingConstructor
// ReSharper disable CppNonExplicitConversionOperator
struct Color {
	vec4 color;
	Color(vec4 color) : color(color) {}
	Color(float r, float g, float b, float a=1.f) : color(vec4(r, g, b, a)) {}
	Color(const Color &other) : color(other.color) {}
	Color(vec3 color) : color(vec4(color, 1.0f)) {}
	Color(ivec3 color) : color(vec4(color.x / 255.f, color.y / 255.f, color.z / 255.f, 1.0f)) {}
	Color(int rgb) : color(rgb >> 16 & 0xFF, rgb >> 8 & 0xFF, rgb & 0xFF, 255) {}
	Color(string hexCode) {
		THROW_IF(hexCode[0] != '#', ValueError, "Hex color code must start with #.");
		THROW_IF(hexCode.length() != 7 && hexCode.length() != 9, ValueError, "Hex color code must be in format #RRGGBB or #RRGGBBAA.");
		int r = std::stoi(hexCode.substr(0, 2), nullptr, 16);
		int g = std::stoi(hexCode.substr(2, 2), nullptr, 16);
		int b = std::stoi(hexCode.substr(4, 2), nullptr, 16);
		int a = hexCode.length() == 7 ? 255 : std::stoi(hexCode.substr(6, 2), nullptr, 16);
		color = vec4(r / 255.f, g / 255.f, b / 255.f, a / 255.f);
	}
	operator vec4() const { return color; }
	float r() const { return color.r; }
	float g() const { return color.g; }
	float b() const { return color.b; }
	float a() const { return color.a; }
	float &operator[] (array_index i) {
		CHECK_OUT_OF_BOUNDS(i, 4);
		if (i == 0) return color.r;
		if (i == 1) return color.g;
		if (i == 2) return color.b;
		return color.a;
	}

};

struct COLOR_PALETTE {
	Color mainColor;
	Color second;
	Color third;
	Color accent;
	Color accent2;

	COLOR_PALETTE(Color mainColor, Color second, Color third, Color accent, Color accent2) : mainColor(mainColor), second(second), third(third), accent(accent), accent2(accent2) {}
	COLOR_PALETTE(vec3 mainColor, vec3 second, vec3 third, vec3 accent, vec3 accent2) : mainColor(mainColor), second(second), third(third), accent(accent), accent2(accent2) {}
	COLOR_PALETTE(ivec3 mainColor, ivec3 second, ivec3 third, ivec3 accent, ivec3 accent2) : mainColor(mainColor), second(second), third(third), accent(accent), accent2(accent2) {}

	Color& operator[] (array_index i) {
		CHECK_OUT_OF_BOUNDS(i, 5);
		if (i == 0) return mainColor;
		if (i == 1) return second;
		if (i == 2) return third;
		if (i == 3) return accent;
		return accent2;
	}

};

struct COLOR_PALETTE10 {
	array<Color, 10> cls;

	explicit COLOR_PALETTE10(const array<Color, 10> &colors) : cls(colors) {}
	COLOR_PALETTE10(COLOR_PALETTE p1, COLOR_PALETTE p2) : cls({p1[0], p1[1], p1[2], p1[3], p1[4], p2[0], p2[1], p2[2], p2[3], p2[4]}) {}
	COLOR_PALETTE10(ivec3 c1, ivec3 c2, ivec3 c3, ivec3 c4, ivec3 c5, ivec3 c6, ivec3 c7, ivec3 c8, ivec3 c9, ivec3 c10)
	: cls({Color(c1), Color(c2), Color(c3), Color(c4), Color(c5), Color(c6), Color(c7), Color(c8), Color(c9), Color(c10)}) {}

	Color& operator[] (int i) { return cls[i]; }
	array<Color, 10>::iterator begin() { return cls.begin(); }
	array<Color, 10>::iterator end() { return cls.end(); }
};

#define BLACK Color(0, 0, 0)
#define WHITE Color(1, 1, 1)
#define RED Color(1, 0, 0)
#define GREEN Color(0, 1, 0)
#define BLUE Color(0, 0, 1)
#define YELLOW Color(1, 1, 0)
#define CYAN Color(0, 1, 1)
#define MAGENTA Color(1, 0, 1)
#define ORANGE Color(1, .5, 0)
#define PINK Color(1, .5, .5)
#define LIGHT_BLUE Color(.5, .5, 1)
#define LIGHT_GREEN Color(.5, 1, .5)
#define LIGHT_RED Color(1, .5, .5)
#define LIGHT_YELLOW Color(1, 1, .5)
#define LIGHT_CYAN Color(.5, 1, 1)
#define LIGHT_MAGENTA Color(1, .5, 1)
#define DARK_PURPLE Color(54, 31, 39)
#define PURPLE Color(82, 25, 69)
#define BROWN_SUGAR Color(169, 113, 7)
#define POWDER_PINK Color(255, 179, 198)
#define DARKDARKBLUE Color(0, 0.016, .045)

#define INTENSE_DARKER_RED_PALLETTE COLOR_PALETTE(ivec3(173, 40, 49), ivec3(128, 14, 19), ivec3(100, 13, 20), ivec3(56, 4, 14), ivec3(37, 9, 2))
#define BLUE_PALLETTE COLOR_PALETTE10(ivec3(226, 239, 246), ivec3(183, 215, 234), ivec3(141, 190, 220), ivec3(100, 164, 206), ivec3(72, 147, 198), ivec3(45, 130, 189), ivec3(36, 118, 177), ivec3(25, 101, 160), ivec3(15, 85, 143), ivec3(0, 58, 112))
#define GREEN_PALLETTE COLOR_PALETTE10(ivec3(219, 234, 215), ivec3(195, 220, 188), ivec3(171, 206, 161), ivec3(148, 193, 134), ivec3(124, 179, 107), ivec3(102, 162, 83), ivec3(85, 135, 69), ivec3(68, 108, 55), ivec3(51, 81, 42), ivec3(32, 51, 26))
#define REDPINK_PALLETTE COLOR_PALETTE10(ivec3(243, 219, 234), ivec3(236, 197, 221), ivec3(204, 154, 181), ivec3(179, 116, 149), ivec3(153, 77, 116), ivec3(128, 39, 84), ivec3(102, 0, 51), ivec3(85, 6, 45), ivec3(67, 12, 39), ivec3(52, 9, 30))
#define GRAY_PALLETTE COLOR_PALETTE10(ivec3(248, 249, 250), ivec3(233, 236, 239), ivec3(222, 226, 230), ivec3(206, 212, 218), ivec3(173, 181, 189), ivec3(108, 117, 125), ivec3(73, 80, 87), ivec3(52, 58, 64), ivec3(33, 37, 41), ivec3(22, 26, 29))
