
#include "specific.hpp"
#include <random>
#include <chrono>

using namespace glm;
using std::vector, std::array;


BoxMesh::BoxMesh(vec3 minCorner, vec3 maxCorner)
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


BallMesh::BallMesh(vec3 center, float radius, int radialSegments, int verticalSegments, vec3 normal)
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

BallMesh::BallMesh(vec3 center, float radius, int radialSegments, int verticalSegments) : BallMesh(center, radius, radialSegments, verticalSegments, vec3(0.0f, 0.0f, 1.0f))
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



CircularRing::CircularRing(vec3 center, float radiusBig, float radiusSmall, int radialSegments, int verticalSegments, vec3 normal)
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

CircularRing::CircularRing(vec3 center, float radiusBig, float radiusSmall, int radialSegments, int verticalSegments) : CircularRing(center, radiusBig, radiusSmall, radialSegments, verticalSegments, vec3(0.0f, 0.0f, 1.0f)){}

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

COLOR_PALETTE::COLOR_PALETTE(vec4 mainColor, vec4 second, vec4 third, vec4 accent, vec4 accent2)
{
	this->mainColor = mainColor;
	this->second = second;
	this->third = third;
	this->accent = accent;
	this->accent2 = accent2;
}

COLOR_PALETTE::COLOR_PALETTE(vec3 mainColor, vec3 second, vec3 third, vec3 accent, vec3 accent2) : COLOR_PALETTE(vec4(mainColor, 1.0f), vec4(second, 1.0f), vec4(third, 1.0f), vec4(accent, 1.0f), vec4(accent2, 1.0f) ) {}

COLOR_PALETTE::COLOR_PALETTE(ivec3 mainColor, ivec3 second, ivec3 third, ivec3 accent, ivec3 accent2)
{
	this->mainColor = vec4(mainColor.x / 255.f, mainColor.y / 255.f, mainColor.z / 255.f, 1.0f);
	this->second = vec4(second.x / 255.f, second.y / 255.f, second.z / 255.f, 1.0f);
	this->third = vec4(third.x / 255.f, third.y / 255.f, third.z / 255.f, 1.0f);
	this->accent = vec4(accent.x / 255.f, accent.y / 255.f, accent.z / 255.f, 1.0f);
	this->accent2 = vec4(accent2.x / 255.f, accent2.y / 255.f, accent2.z / 255.f, 1.0f);
}

std::vector<vec4> COLOR_PALETTE::colors()
{
    return std::vector({ mainColor, second, third, accent, accent2 });
}

vec4 COLOR_PALETTE::operator[](int i)
{
    return colors()[i];
}

COLOR_PALETTE10::COLOR_PALETTE10(ivec3 c1, ivec3 c2, ivec3 c3,
                                 ivec3 c4, ivec3 c5, ivec3 c6,
                                 ivec3 c7, ivec3 c8, ivec3 c9,
                                 ivec3 c10) {
  cls = {vec4(c1.x / 255.f, c1.y / 255.f, c1.z / 255.f, 1.0f),
         vec4(c2.x / 255.f, c2.y / 255.f, c2.z / 255.f, 1.0f),
         vec4(c3.x / 255.f, c3.y / 255.f, c3.z / 255.f, 1.0f),
         vec4(c4.x / 255.f, c4.y / 255.f, c4.z / 255.f, 1.0f),
         vec4(c5.x / 255.f, c5.y / 255.f, c5.z / 255.f, 1.0f),
         vec4(c6.x / 255.f, c6.y / 255.f, c6.z / 255.f, 1.0f),
         vec4(c7.x / 255.f, c7.y / 255.f, c7.z / 255.f, 1.0f),
         vec4(c8.x / 255.f, c8.y / 255.f, c8.z / 255.f, 1.0f),
         vec4(c9.x / 255.f, c9.y / 255.f, c9.z / 255.f, 1.0f),
         vec4(c10.x / 255.f, c10.y / 255.f, c10.z / 255.f, 1.0f)};

}
SmoothParametricPlaneCurve circle(float r, vec2 center, float eps) {
        return SmoothParametricPlaneCurve(
                [center, r](float t) {return center + r * vec2(cos(t), sin(t)); },
                [r](float t) {return r * vec2(-sin(t), cos(t)); },
                [r](float t) {return r * vec2(-cos(t), -sin(t)); },
                0, TAU, true, eps);
}


SmoothParametricPlaneCurve ellipse(float a, float b, vec2 center, float eps) {
    return SmoothParametricPlaneCurve(
                [center, a, b](float t) {return center + vec2(a*cos(t), b*sin(t)); },
                [a, b](float t) {return vec2(-a*sin(t),  b*cos(t)); },
                [a, b](float t) {return vec2(-a*cos(t), -b*sin(t)); },
                0, TAU, true, eps);
}
SmoothParametricPlaneCurve epitrochoid(float r, float R, float d, float eps) {
    return SmoothParametricPlaneCurve(
                [r, R, d](float t) {return vec2((r+R)*cos(t) - d*cos((r+R)/r*t), (r+R)*sin(t) - d*sin((r+R)/r*t)); },
                [r, R, d](float t) {return vec2(-(r+R)*sin(t) + (r+R)/r*d*sin((r+R)/r*t), (r+R)*cos(t) - (r+R)/r*d*cos((r+R)/r*t)); },
                [r, R, d](float t) {return vec2(-(r+R)*cos(t) + (r+R)*(r+R)/r/r*d*cos((r+R)/r*t), -(r+R)*sin(t) + (r+R)*(r+R)/r/r*d*sin((r+R)/r*t)); },
                0, TAU, true, eps);
}

SmoothParametricPlaneCurve hypotrochoid(float r, float R, float d, float eps) {
    return SmoothParametricPlaneCurve(
                [r, R, d](float t) {return vec2((R-r)*cos(t) + d*cos((R-r)/r*t), (R-r)*sin(t) - d*sin((R-r)/r*t)); },
                [r, R, d](float t) {return vec2(-(R-r)*sin(t) - d*(R-r)/r*sin((R-r)/r*t), (R-r)*cos(t) - d*(R-r)/r*cos((R-r)/r*t)); },
                [r, R, d](float t) {return vec2(-(R-r)*cos(t) - d*(R-r)/r*(R-r)/r*cos((R-r)/r*t), -(R-r)*sin(t) + d*(R-r)/r*(R-r)/r*sin((R-r)/r*t)); },
                0, TAU, true, eps);
}

SmoothParametricCurve circle(float r, vec3 center, vec3 v1, vec3 v2, float eps) {
    return circle(r, PLANE_ORIGIN, eps).embedding(v1, v2, center);
}
SmoothParametricCurve VivaniCurve(float r, float eps) {
    return SmoothParametricCurve(
                [r](float t) {return vec3(r*(1+cos(t)), r*sin(t), 2*r*sin(t/2)); },
                [r](float t) {return vec3(-r*sin(t), r*cos(t), r*cos(t/2)); },
                [r](float t) {return vec3(-r*cos(t), -r*sin(t), -r/2*sin(t/2)); },
                -TAU, TAU, true, eps);
}

SmoothParametricPlaneCurve LissajousCurve(float a, float b, float delta, float r1, float r2, float eps) {
    return SmoothParametricPlaneCurve(
                [r1, r2, a, b, delta](float t) {return vec2(r1*sin(a*t+delta), r2*cos(b*t)); },
                [r1, r2, a, b, delta](float t) {return vec2(a*r1*cos(a*t+delta), -b*r2*sin(b*t)); },
                [r1, r2, a, b, delta](float t) {return vec2(-a*a*r1*sin(a*t+delta), -b*b*r2*cos(b*t)); },
                0, TAU, true, eps);
}

SuperCurve circle(float r, std::function<float(float)> w, std::function<MaterialPhong(float)> mat, int n, vec3 center, vec3 v1, vec3 v2, float eps) {
    return SuperCurve(circle(r, center, v1, v2, eps), w, mat, n, 0, TAU, true);
}



SmoothParametricCurve sphericalSpiral(float a, float t_max, PolyGroupID id, float eps) {
    return SmoothParametricCurve([a](float t) {return vec3(cos(t), sin(t), -a*t)/sqrt(1+a*a*t*t); },
                                        -t_max, t_max, false, id, eps);
}

SmoothParametricCurve sphericalSpiral(float a, float r, float t_max, PolyGroupID id, float eps) {
    return SmoothParametricCurve([a, r](float t) {return vec3(cos(t), sin(t), -a*t)*r/sqrt(1+a*a*t*t); },
                                        -t_max, t_max, false, id, eps);
}

WeakSuperMesh singleTrig(vec3 v0, vec3 v1, vec3 v2, MaterialPhong &material, PolyGroupID id) {
    vec3 n = normalize(cross(v1 - v0, v2 - v0));
    vector nodes = {Vertex(v0, vec2(0, 0), n, BLACK, material),
                    Vertex(v1, vec2(0, 1), n, BLACK, material),
                    Vertex(v2, vec2(1, 1), n, BLACK, material)};
    return WeakSuperMesh(nodes, {ivec3(0, 1, 2)}, id);
}


WeakSuperMesh singleTrig(vec3 v0, vec3 v1, vec3 v2, MaterialPhong &material1, MaterialPhong &material2, MaterialPhong &material3,
                                PolyGroupID id) {
    vec3 n = normalize(cross(v1 - v0, v2 - v0));
    vector nodes = {Vertex(v0, vec2(0, 0), n, BLACK, material1),
                    Vertex(v1, vec2(0, 1), n, BLACK, material2),
                    Vertex(v2, vec2(1, 1), n, BLACK, material3)};
    return WeakSuperMesh(nodes, {ivec3(0, 1, 2)}, id);
}
WeakSuperMesh singleQuadShadeSmooth(vec3 outer1, vec3 inner1, vec3 inner2, vec3 outer2, MaterialPhong &material,
                                    std::variant<int, std::string> id) {
    vec3 n1out = normalize(cross(outer1 - inner1, outer1 - inner2));
    vec3 n2out = normalize(cross(outer2 - inner2, outer2 - inner1));
    float w1 = norm(cross(inner1-outer1, inner1-inner2));
    float w2 = norm(cross(inner2-outer2, inner2-inner1));
    vec3 nin = normalize(w1*inner1 + w2*inner2);

    vector nodes = {Vertex(outer1, vec2(0, 0), n1out, BLACK, material),
                    Vertex(inner1, vec2(0, 1), nin, BLACK, material),
                    Vertex(inner2, vec2(1, 0), nin, BLACK, material),
                    Vertex(outer2, vec2(1, 1), n2out, BLACK, material)};
    return WeakSuperMesh(nodes, {ivec3(0, 1, 2), ivec3(3, 1, 2)}, id);
}

WeakSuperMesh singleQuadShadeFlat(glm::vec3 outer1, glm::vec3 inner1, glm::vec3 inner2, glm::vec3 outer2, MaterialPhong &material,
                                  std::variant<int, std::string> id) {
    vec3 n1 = normalize(cross(outer1 - inner1, outer1 - inner2));
    vec3 n2 = normalize(cross(outer2 - inner2, outer2 - inner1));

    vector nodes = {Vertex(outer1, vec2(0, 0), n1, BLACK, material),
                    Vertex(inner1, vec2(0, 1), n1, BLACK, material),
                    Vertex(inner2, vec2(1, 1), n1, BLACK, material),
                    Vertex(inner1, vec2(0, 1), n2, BLACK, material),
                    Vertex(inner2, vec2(1, 1), n2, BLACK, material),
                    Vertex(outer2, vec2(1, 1), n2, BLACK, material)};
    return WeakSuperMesh(nodes, {ivec3(0, 1, 2), ivec3(3, 4, 5)}, id);
}

WeakSuperMesh singleQuadShadeFlat(glm::vec3 outer1, glm::vec3 inner1, glm::vec3 inner2, glm::vec3 outer2, MaterialPhong &material1,
                                  MaterialPhong &material2, std::variant<int, std::string> id) {
    vec3 n1 = normalize(cross(outer1 - inner1, outer1 - inner2));
    vec3 n2 = normalize(cross(outer2 - inner2, outer2 - inner1));

    vector nodes = {Vertex(outer1, vec2(0, 0), n1, BLACK, material1), Vertex(inner1, vec2(0, 1), n1, BLACK, material1),
                    Vertex(inner2, vec2(1, 1), n1, BLACK, material1), Vertex(inner1, vec2(0, 1), n2, BLACK, material2),
                    Vertex(inner2, vec2(1, 1), n2, BLACK, material2), Vertex(outer2, vec2(1, 1), n2, BLACK, material2)};
    return WeakSuperMesh(nodes, {ivec3(0, 1, 2), ivec3(3, 4, 5)}, id);
}

SmoothParametricSurface sphere(float r, glm::vec3 center, float eps) {
    return SmoothParametricSurface([r, center](float t, float u) {
        return center + vec3(cos(t) * sin(u), sin(t) * sin(u), cos(u))*r;
    }, vec2(0, TAU*1.1), vec2(0, PI), true, false, .01);
}





