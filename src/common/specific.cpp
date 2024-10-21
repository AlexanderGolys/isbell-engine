
#include "specific.hpp"
#include <random>
#include <chrono>

using namespace glm;
using std::vector, std::array;


BoxMesh::BoxMesh(glm::vec3 minCorner, glm::vec3 maxCorner)
{
	this->triangles = std::vector<TriangleR3>();

	vector<vec3> vert = {
		vec3(minCorner.x, minCorner.y, minCorner.z),
		vec3(minCorner.x, minCorner.y, maxCorner.z),
		vec3(minCorner.x, maxCorner.y, minCorner.z),
		vec3(minCorner.x, maxCorner.y, maxCorner.z),
		vec3(maxCorner.x, minCorner.y, minCorner.z),
		vec3(maxCorner.x, minCorner.y, maxCorner.z),
		vec3(maxCorner.x, maxCorner.y, minCorner.z),
		vec3(maxCorner.x, maxCorner.y, maxCorner.z) };

	vec4 color = vec4(0.0f, 0.0f, 0.0f, 1.0f);

	this->triangles.push_back(TriangleR3({ vert[0], vert[1], vert[2] }, color));
	this->triangles.push_back(TriangleR3({vert[1], vert[3], vert[2]}, color));
	this->triangles.push_back(TriangleR3({vert[4], vert[5], vert[6]}, color));
	this->triangles.push_back(TriangleR3({vert[5], vert[7], vert[6]}, color));
	this->triangles.push_back(TriangleR3({vert[0], vert[2], vert[4]}, color));
	this->triangles.push_back(TriangleR3({vert[2], vert[6], vert[4]}, color));
	this->triangles.push_back(TriangleR3({vert[1], vert[5], vert[3]}, color));
	this->triangles.push_back(TriangleR3({vert[5], vert[7], vert[3]}, color));
	this->triangles.push_back(TriangleR3({vert[0], vert[4], vert[1]}, color));
	this->triangles.push_back(TriangleR3({vert[4], vert[5], vert[1]}, color));
	this->triangles.push_back(TriangleR3({vert[2], vert[3], vert[6]}, color));
	this->triangles.push_back(TriangleR3({vert[3], vert[7], vert[6]}, color));
										 
	for (TriangleR3 t : this->triangles) {
		t.recalculateNormal();
	}
}





void BoxMesh::assignFaceColors(array<vec4, 6> faceColors)
{
	for (int i = 0; i < 6; i++) {
		this->triangles.at(i * 2).vertexColors = { faceColors[i] , faceColors[i], faceColors[i] };
		this->triangles.at(i * 2 + 1).vertexColors = { faceColors[i] , faceColors[i], faceColors[i] };
	}
}

BoxMesh::BoxMesh(vec3 minCorner, vec3 maxCorner, array<vec4, 6> faceColors) : BoxMesh(minCorner, maxCorner)
{
	assignFaceColors(faceColors);
}

BoxMesh::BoxMesh(vec3 minCorner, vec3 maxCorner, vec4 color) : BoxMesh(minCorner, maxCorner)
{
	assignFaceColors({ color, color, color, color, color, color });
}

void BoxMesh::randomizeFaceColors()
{
	auto seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine generator(seed);
	std::uniform_real_distribution<double> distribution(0.0, 1.0);
	array<vec4, 6> colors = array<vec4, 6>();
	for (int i = 0; i < 6; i++) {
		colors[i] = vec4(distribution(generator), distribution(generator), distribution(generator), 1.0f);
	}
	assignFaceColors(colors);

}


BallMesh::BallMesh(vec3 center, float radius, int radialSegments, int verticalSegments, glm::vec3 normal)
{
	this->center = center;
	this->radius = radius;
	this->bdVertices = vector<vec3>();
	auto tangents = orthogonalComplementBasis(normal);
	vec3 t1 = tangents.first;
	vec3 t2 = tangents.second;
	auto trs = vector<TriangleR3>();
	for (int i = 0; i < radialSegments; i++)
		for (int j = 0; j < verticalSegments; j++)
		{
			float r1 = radius * j * 1.f / verticalSegments;
			float r2 = radius * (j + 1) * 1.f / verticalSegments;
			float theta1 = TAU * i * 1.f / radialSegments;
			float theta2 = TAU * (i + 1)*1.f / radialSegments;
			vec3 p1 = center + r1 * t1*cos(theta1)  + t2*r1 * sin(theta1);
			vec3 p2 = center + r1 * t1*cos(theta2) + t2*r1 * sin(theta2);
			vec3 p3 = center + r2 * t1*cos(theta1) + t2*r2 * sin(theta1);
			vec3 p4 = center + r2 * t1*cos(theta2) + t2*r2 * sin(theta2);
			auto tr1 = TriangleR3({ p1, p2, p3 }, normal);
			auto tr2 = TriangleR3({ p2, p4, p3 }, normal);
			trs.push_back(tr1);
			trs.push_back(tr2);
			if (j == verticalSegments - 1)
			{
				this->bdVertices.push_back(p3);
			}
		}
	this->triangles = trs;
}

BallMesh::BallMesh(glm::vec3 center, float radius, int radialSegments, int verticalSegments) : BallMesh(center, radius, radialSegments, verticalSegments, vec3(0.0f, 0.0f, 1.0f))
{}

BallMesh::BallMesh(int radialSegments, int verticalSegments) : BallMesh(vec3(0.0f, 0.0f, 0.0f), 1.0f, radialSegments, verticalSegments) {}


BallMesh::BallMesh() : BallMesh(20, 3) { }


TriangularMesh BallMesh::extrudeAlongNormal(float h)
{
	vec3 normal = this->triangles[0].normals[0];

	auto trs = vector<TriangleR3>();
	for (TriangleR3 tr : triangles)
	{
		TriangleR3 tri = TriangleR3(tr.vertices, tr.normals, tr.vertexColors, tr.uvs);
		trs.push_back(tri);
		
		TriangleR3 tri2 = TriangleR3(tr.vertices, tr.normals, tr.vertexColors, tr.uvs) + (normal * h);
		tri2.flipNormals();
		trs.push_back(tri2);
	}
	
	for (int i = 0; i < bdVertices.size(); i++)
	{
		auto v = bdVertices[i];
		auto w = bdVertices[(i + 1) % bdVertices.size()];
		auto tr = TriangleR3({ v, w, w +normal*h });
		auto tr2 = TriangleR3({ v, w + normal*h, v + normal*h });
		tr.recalculateNormal();
		tr2.recalculateNormal();
		trs.push_back(tr);
		trs.push_back(tr2);
	}
	
	return TriangularMesh(trs);
	
}

TriangularMesh CircularRing::extrudeAlongNormal(float h)
{
	vec3 normal = this->triangles[0].normals[0];

	auto trs = vector<TriangleR3>();
	for (TriangleR3 tr : triangles)
	{
		TriangleR3 tri = TriangleR3(tr.vertices, tr.normals, tr.vertexColors, tr.uvs);
		trs.push_back(tri);

		TriangleR3 tri2 = TriangleR3(tr.vertices, tr.normals, tr.vertexColors, tr.uvs) + (normal * h);
		tri2.flipNormals();
		trs.push_back(tri2);
	}

	for (int i = 0; i < bdVerticesInner.size(); i++)
	{
		auto v = bdVerticesInner[i];
		auto w = bdVerticesInner[(i + 1) % bdVerticesInner.size()];
		auto tr = TriangleR3({ v, w, w + normal * h });
		auto tr2 = TriangleR3({ v, w + normal * h, v + normal * h });
		tr.recalculateNormal();
		tr2.recalculateNormal();
		trs.push_back(tr);
		trs.push_back(tr2);
	}

	for (int i = 0; i < bdVerticesOuter.size(); i++)
	{
		auto v = bdVerticesOuter[i];
		auto w = bdVerticesOuter[(i + 1) % bdVerticesOuter.size()];
		auto tr = TriangleR3({ v, w, w + normal * h });
		auto tr2 = TriangleR3({ v, w + normal * h, v + normal * h });
		tr.recalculateNormal();
		tr2.recalculateNormal();
		trs.push_back(tr);
		trs.push_back(tr2);
	}

	return TriangularMesh(trs);
}



CircularRing::CircularRing(vec3 center, float radiusBig, float radiusSmall, int radialSegments, int verticalSegments, glm::vec3 normal)
{
	this->center = center;
	this->radiusBig = radiusBig;
	this->radiusSmall = radiusSmall;
	auto tangents = orthogonalComplementBasis(normal);
	vec3 t1 = tangents.first;
	vec3 t2 = tangents.second;
	auto trs = vector<TriangleR3>();
	for (int i = 0; i < radialSegments; i++)
	{
		for (int j = 0; j < verticalSegments; j++)
		{
			float r1 = lerp(radiusSmall, radiusBig, j * 1.f / verticalSegments);
			float r2 = lerp(radiusSmall, radiusBig, (j + 1) * 1.f / verticalSegments);
			float theta1 = TAU * i / radialSegments;
			float theta2 = TAU * (i + 1) / radialSegments;
			vec3 p1 = center + r1 * cos(theta1) * t1 + r1 * sin(theta1) * t2;
			vec3 p2 = center + r1 * cos(theta2) * t1 + r1 * sin(theta2) * t2;
			vec3 p3 = center + r2 * cos(theta1) * t1 + r2 * sin(theta1) * t2;
			vec3 p4 = center + r2 * cos(theta2) * t1 + r2 * sin(theta2) * t2;
			auto tr1 = TriangleR3({ p1, p2, p3 }, normal);
			auto tr2 = TriangleR3({ p2, p4, p3 }, normal);
			trs.push_back(tr1);
			trs.push_back(tr2);
			if (j == verticalSegments - 1)
			{
				this->bdVerticesOuter.push_back(p3);
			}
			if (j == 0)
			{
				this->bdVerticesInner.push_back(p1);
			}
		}
	}
	this->triangles = trs;
}

CircularRing::CircularRing(glm::vec3 center, float radiusBig, float radiusSmall, int radialSegments, int verticalSegments) : CircularRing(center, radiusBig, radiusSmall, radialSegments, verticalSegments, vec3(0.0f, 0.0f, 1.0f)){}

CircularRing::CircularRing(int radialSegments, int verticalSegments) : CircularRing(vec3(0.0f, 0.0f, 0.0f), 1.0f, 3.f, radialSegments, verticalSegments){}


CircularRing::CircularRing() : CircularRing(20, 3){}


PlanarUnitDisk::PlanarUnitDisk(int radial_res, int vertical_res)
	{
		auto trng = vector<TriangleR2>();
		vector<vec2> bd = vector<vec2>();
		for (int i = 0; i < radial_res; i++)
		{
			float theta1 = TAU * i / radial_res;
			float theta2 = TAU * (i + 1) / radial_res;
			vec2 p1 = vec2(cos(theta1) / vertical_res, sin(theta1) / vertical_res);
			vec2 p2 = vec2(cos(theta2) / vertical_res, sin(theta2) / vertical_res);
			trng.push_back(TriangleR2({ vec2(0, 0), p1, p2 }, { vec2(0, 0), p1, p2 }));
		}
		for (int i = 0; i < radial_res; i++)
		{
			for (int h = 1; h < vertical_res; h++)
			{
				float theta1 = TAU * i / radial_res;
				float theta2 = TAU * (i + 1) / radial_res;
				float h1 = h * 1.f / vertical_res;
				float h2 = (h + 1) * 1.f / vertical_res;
				vec2 p1 = vec2(cos(theta1) * h1, sin(theta1) * h1);
				vec2 p2 = vec2(cos(theta2) * h1, sin(theta2) * h1);
				vec2 p3 = vec2(cos(theta1) * h2, sin(theta1) * h2);
				vec2 p4 = vec2(cos(theta2) * h2, sin(theta2) * h2);

				trng.push_back(TriangleR2(p1, p2, p3));
				trng.push_back(TriangleR2(p4, p3, p2));
				if (h == vertical_res - 1)
				{
					bd.push_back(p3);
				}
			}
			
			
		}
		this->boundaries = { bd };
		this->boundaryCyclic = { true };
		this->triangles = trng;
	}





PlanarRing::PlanarRing(int radial_res, int vertical_res,vec2 center, float radiusBig, float radiusSmall)
	{
		auto trng = vector<TriangleR2>();
		trng.reserve(2 * radial_res * vertical_res);
		this->center = center;
		this->radiusBig = radiusBig;
		this->radiusSmall = radiusSmall;
		vector<vec2> bd = vector<vec2>();
		vector<vec2> bd2 = vector<vec2>();

		for (int i = 0; i < radial_res; i++)
		{
			for (int h = 0; h < vertical_res; h++)
			{
				float theta1 = TAU * i / radial_res;
				float theta2 = TAU * (i + 1) / radial_res;
				float h1 = lerp(radiusSmall, radiusBig, h * 1.f / vertical_res);
				float h2 = lerp(radiusSmall, radiusBig, (h + 1) * 1.f / vertical_res);
				vec2 p1 = vec2(cos(theta1) * h1, sin(theta1) * h1) + center;
				vec2 p2 = vec2(cos(theta2) * h1, sin(theta2) * h1) + center;
				vec2 p3 = vec2(cos(theta1) * h2, sin(theta1) * h2) + center;
				vec2 p4 = vec2(cos(theta2) * h2, sin(theta2) * h2) + center;

				trng.push_back(TriangleR2({ p1, p2, p3 }, { p1, p2, p3 }));
				trng.push_back(TriangleR2({ p2, p4, p3 }, { p2, p4, p3 }));
				if (h == vertical_res - 1)
				{
					bd.push_back(p3);
				}
				if (h == 0)
				{
					bd2.push_back(p1);
				}
			}


		}
		this->boundaries = { bd2, bd };
		this->boundaryCyclic = { true , true};
		this->triangles = trng;
	}

	

void PlanarRing::addRotationalField(float power)
	{
		VectorFieldR2 vf = VectorFieldR2(std::function<vec2(vec2)>([power, this](vec2 p) {return vec2(-(p - this->center).y, (p - this->center).x) * power; }));
		this->addVectorField(vf);
	}

    
PlanarConvexPolygon::PlanarConvexPolygon(vector<vec2> verts)
{
	auto trng = vector<TriangleR2>();
	float maxRadius = 0;
	vector<vec2> bd = vector<vec2>();
	vec2 center = mean(verts);
	for (int i = 0; i < verts.size(); i++)
	{
		if (norm(verts[i] - center) > maxRadius)
		{
			maxRadius = norm(verts[i] - center);
		}
	}
	for (int i = 0; i < verts.size(); i++)
	{
		vec2 p1 = verts[i];
		vec2 p2 = verts[(i + 1) % verts.size()];

		trng.push_back(TriangleR2({ center, p1, p2 }, { vec2(0, 0), (p1 - center) / (2.f * maxRadius), (p2 - center) / (2.f * maxRadius) }));
		bd.push_back(p1);
	}
	this->vertices = verts;
	this->boundaries = { bd };
	this->boundaryCyclic = { true };
	this->triangles = trng;
}

COLOR_PALETTE::COLOR_PALETTE(glm::vec4 mainColor, glm::vec4 second, glm::vec4 third, glm::vec4 accent, glm::vec4 accent2)
{
	this->mainColor = mainColor;
	this->second = second;
	this->third = third;
	this->accent = accent;
	this->accent2 = accent2;
}

COLOR_PALETTE::COLOR_PALETTE(glm::vec3 mainColor, glm::vec3 second, glm::vec3 third, glm::vec3 accent, glm::vec3 accent2) : COLOR_PALETTE(vec4(mainColor, 1.0f), vec4(second, 1.0f), vec4(third, 1.0f), vec4(accent, 1.0f), vec4(accent2, 1.0f) ) {}

COLOR_PALETTE::COLOR_PALETTE(glm::ivec3 mainColor, glm::ivec3 second, glm::ivec3 third, glm::ivec3 accent, glm::ivec3 accent2)
{
	this->mainColor = vec4(mainColor.x / 255.f, mainColor.y / 255.f, mainColor.z / 255.f, 1.0f);
	this->second = vec4(second.x / 255.f, second.y / 255.f, second.z / 255.f, 1.0f);
	this->third = vec4(third.x / 255.f, third.y / 255.f, third.z / 255.f, 1.0f);
	this->accent = vec4(accent.x / 255.f, accent.y / 255.f, accent.z / 255.f, 1.0f);
	this->accent2 = vec4(accent2.x / 255.f, accent2.y / 255.f, accent2.z / 255.f, 1.0f);
}

std::vector<glm::vec4> COLOR_PALETTE::colors()
{
    return std::vector<glm::vec4>({ mainColor, second, third, accent, accent2 });
}

glm::vec4 COLOR_PALETTE::operator[](int i)
{
    return colors()[i];
}

COLOR_PALETTE10::COLOR_PALETTE10(glm::ivec3 c1, glm::ivec3 c2, glm::ivec3 c3, glm::ivec3 c4, glm::ivec3 c5,
	glm::ivec3 c6, glm::ivec3 c7, glm::ivec3 c8, glm::ivec3 c9, glm::ivec3 c10) {
	cls = {vec4(c1.x/255.f, c1.y/255.f, c1.z/255.f, 1.0f), vec4(c2.x/255.f, c2.y/255.f, c2.z/255.f, 1.0f),
		vec4(c3.x/255.f, c3.y/255.f, c3.z/255.f, 1.0f), vec4(c4.x/255.f, c4.y/255.f, c4.z/255.f, 1.0f),
		vec4(c5.x/255.f, c5.y/255.f, c5.z/255.f, 1.0f), vec4(c6.x/255.f, c6.y/255.f, c6.z/255.f, 1.0f),
		vec4(c7.x/255.f, c7.y/255.f, c7.z/255.f, 1.0f), vec4(c8.x/255.f, c8.y/255.f, c8.z/255.f, 1.0f),
		vec4(c9.x/255.f, c9.y/255.f, c9.z/255.f, 1.0f), vec4(c10.x/255.f, c10.y/255.f, c10.z/255.f, 1.0f)};
}

SuperCurve circle(vec3 v1, vec3 v2, vec3 center, float r, std::function<float(float)> w, std::function<MaterialPhong(float)> mat, int n, float eps) {
	vec3 v1n = normalize(v1);
	vec3 v2n = normalize(v2);
	SmoothParametricCurve curve = SmoothParametricCurve([v1n, v2n, center, r](float t) {return center + r * cos(t) * v1n + r * sin(t) * v2n; },
	[v1n, v2n, r](float t) {return -r * sin(t) * v1n + r * cos(t) * v2n; },
	[v1n, v2n, r](float t) {return -r * cos(t) * v1n - r * sin(t) * v2n; }, eps);
	return SuperCurve(curve, w, mat, 0, TAU, n);
}

SuperPencilCurve circlePencil(vec3 v1, vec3 v2, vec3 center, float r, std::function<float(float)> w, std::function<MaterialPhong(float)> mat, int n, float eps) {
	vec3 v1n = normalize(v1);
	vec3 v2n = normalize(v2);
	SmoothParametricCurve curve = SmoothParametricCurve([v1n, v2n, center, r](float t) {return center + r * cos(t) * v1n + r * sin(t) * v2n; },
	[v1n, v2n, r](float t) {return -r * sin(t) * v1n + r * cos(t) * v2n; },
	[v1n, v2n, r](float t) {return -r * cos(t) * v1n - r * sin(t) * v2n; }, eps);
	return SuperPencilCurve(curve, w, mat, 0, TAU, n);
}

