#ifndef SPECIFIC_HPP
#define SPECIFIC_HPP

#include "geometry.hpp"

class BoxMesh : public TriangularMesh {
public:
	BoxMesh(glm::vec3 minCorner, glm::vec3 maxCorner);
	BoxMesh(glm::vec3 minCorner, glm::vec3 maxCorner, glm::vec4 color);
	BoxMesh(glm::vec3 minCorner, glm::vec3 maxCorner, std::array<glm::vec4, 6> faceColors);


	void assignFaceColors(std::array<glm::vec4, 6> faceColors);
	void randomizeFaceColors();
};

class BallMesh : public TriangularMesh {
public:
	glm::vec3 center;
	float radius;
	std::vector<glm::vec3> bdVertices;
	BallMesh(glm::vec3 center, float radius, int radialSegments, int verticalSegments, glm::vec3 normal);
	BallMesh(glm::vec3 center, float radius, int radialSegments, int verticalSegments=3);
	BallMesh(int radialSegments, int verticalSegments=3);
	BallMesh();
	TriangularMesh extrudeAlongNormal(float h);
};

class CircularRing : public TriangularMesh {
public:
	glm::vec3 center;
	float radiusBig;
	float radiusSmall;
	std::vector<glm::vec3> bdVerticesInner;
	std::vector<glm::vec3> bdVerticesOuter;


	CircularRing(glm::vec3 center, float radiusBig, float radiusSmall, int radialSegments, int verticalSegments, glm::vec3 normal);
	CircularRing(glm::vec3 center, float radiusBig, float radiusSmall, int radialSegments, int verticalSegments = 3);
	CircularRing(int radialSegments, int verticalSegments = 3);
	CircularRing();
	TriangularMesh extrudeAlongNormal(float h);
};


class PlanarUnitDisk : public PlanarMeshWithBoundary {
public:
	PlanarUnitDisk(int radial_res, int vertical_res);
};

class PlanarConvexPolygon: public PlanarMeshWithBoundary{
public:
	std::vector<glm::vec2> vertices;
	PlanarConvexPolygon(std::vector<glm::vec2> verts);
};

class PlanarRing : public PlanarMeshWithBoundary {
public:
	glm::vec2 center;
	float radiusBig;
	float radiusSmall;
	PlanarRing(int radial_res, int vertical_res, glm::vec2 center, float radiusBig, float radiusSmall);
	void addRotationalField(float power);
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

	COLOR_PALETTE10(std::array<glm::vec4, 10> colors) : cls(colors) {}
	COLOR_PALETTE10(COLOR_PALETTE p1, COLOR_PALETTE p2) : cls({p1[0], p1[1], p1[2], p1[3], p1[4], p2[0], p2[1], p2[2], p2[3], p2[4]}) {}
	COLOR_PALETTE10(glm::ivec3 c1, glm::ivec3 c2, glm::ivec3 c3, glm::ivec3 c4, glm::ivec3 c5, glm::ivec3 c6, glm::ivec3 c7, glm::ivec3 c8, glm::ivec3 c9, glm::ivec3 c10);
	std::vector<glm::vec4> colors() const { return std::vector<glm::vec4>(cls.begin(), cls.end()); }
	glm::vec4 operator[] (int i) const { return cls[i]; }
};

SuperCurve circle(glm::vec3 v1, glm::vec3 v2, glm::vec3 center, float r, std::function<float(float)> w, std::function<MaterialPhong(float)> mat, int n, float eps=.01);
SuperPencilCurve circlePencil(glm::vec3 v1, glm::vec3 v2, glm::vec3 center, float r, std::function<float(float)> w, std::function<MaterialPhong(float)> mat, int n, float eps=.01);
SuperPencilCurve ellipse(float a, float b, glm::vec3 v1, glm::vec3 v2, glm::vec3 center, std::function<float(float)> w, std::function<MaterialPhong(float)> mat, int n, float eps=.01);
SuperPencilCurve cycloid(float a, float b, float gamma, glm::vec3 v1, glm::vec3 v2, glm::vec3 center, std::function<float(float)> w, std::function<MaterialPhong(float)> mat, int n, float eps=.01);

#endif

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