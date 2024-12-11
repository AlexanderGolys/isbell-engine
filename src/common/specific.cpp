
#include "specific.hpp"

#include <chrono>
#include <chrono>
#include <chrono>
#include <chrono>
#include <chrono>
#include <chrono>
#include <random>
#include <chrono>
#include <iosfwd>
#include <iosfwd>
#include <iosfwd>
#include <iosfwd>
#include <iosfwd>
#include <iosfwd>
#include <vector>
#include <vector>
#include <vector>
#include <vector>
#include <vector>
#include <vector>
#include <glm/detail/_vectorize.hpp>
#include <glm/detail/_vectorize.hpp>
#include <glm/detail/_vectorize.hpp>
#include <glm/detail/_vectorize.hpp>
#include <glm/detail/_vectorize.hpp>
#include <glm/detail/_vectorize.hpp>
#include <glm/detail/_vectorize.hpp>
#include <glm/detail/_vectorize.hpp>
#include <glm/detail/_vectorize.hpp>
#include <glm/detail/_vectorize.hpp>

using namespace glm;
using std::vector, std::array;


PlanarMeshWithBoundary PlanarUnitDisk(int radial_res, int vertical_res)
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

                trng.emplace_back(p1, p2, p3);
                trng.push_back(TriangleR2(p4, p3, p2));
                if (h == vertical_res - 1)
                    bd.push_back(p3);
            }
        return PlanarMeshWithBoundary(trng, { bd }, { true });
	}





PlanarMeshWithBoundary PlanarRing(int radial_res, int vertical_res,vec2 center, float radiusBig, float radiusSmall)
	{
		auto trng = vector<TriangleR2>();
		trng.reserve(2 * radial_res * vertical_res);
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
		return PlanarMeshWithBoundary(trng, { bd, bd2 }, { true, true });
	}


    
PlanarMeshWithBoundary PlanarConvexPolygon(const vector<vec2> &verts)
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
    return PlanarMeshWithBoundary(trng, { bd }, { true });
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

SmoothParametricPlaneCurve epicycloid(float r, float R, float eps) { return epitrochoid(r, R, r, eps); }

SmoothParametricPlaneCurve cardioid(float r, float eps) { return epicycloid(r, r, eps); }

SmoothParametricPlaneCurve nephroid(float r, float eps) { return epicycloid(r, 2*r, eps); }

SmoothParametricPlaneCurve trefoiloid(float r, float eps) { return epicycloid(r, 3*r, eps); }

SmoothParametricPlaneCurve quatrefoloid(float r, float eps) { return epicycloid(r, 4*r, eps); }

SmoothParametricPlaneCurve hypotrochoid(float r, float R, float d, float eps) {
    return SmoothParametricPlaneCurve(
                [r, R, d](float t) {return vec2((R-r)*cos(t) + d*cos((R-r)/r*t), (R-r)*sin(t) - d*sin((R-r)/r*t)); },
                [r, R, d](float t) {return vec2(-(R-r)*sin(t) - d*(R-r)/r*sin((R-r)/r*t), (R-r)*cos(t) - d*(R-r)/r*cos((R-r)/r*t)); },
                [r, R, d](float t) {return vec2(-(R-r)*cos(t) - d*(R-r)/r*(R-r)/r*cos((R-r)/r*t), -(R-r)*sin(t) + d*(R-r)/r*(R-r)/r*sin((R-r)/r*t)); },
                0, TAU, true, eps);
}

SmoothParametricPlaneCurve hypocycloid(float r, float R, float eps) { return hypotrochoid(r, R, r, eps); }

SmoothParametricPlaneCurve astroid(float r, float eps) { return hypocycloid(r, 4*r, eps); }

SmoothParametricPlaneCurve deltoid(float r, float eps) { return hypocycloid(r, 3*r, eps); }

SmoothParametricPlaneCurve pentoid(float r, float eps) { return hypocycloid(r, 5*r, eps); }

SmoothParametricPlaneCurve exoid(float r, float eps) { return hypocycloid(r, 6*r, eps); }

SmoothParametricCurve circle(float r, vec3 center, vec3 v1, vec3 v2, float eps) {
    return circle(r, PLANE_ORIGIN, eps).embedding(v1, v2, center);
}

WeakSuperMesh singleTrig(vec3 v0, vec3 v1, vec3 v2, std::variant<int, std::string> id) {
	vec3 normal          = normalize(cross(v1 - v0, v2 - v0));
	vector<Vertex> nodes = {Vertex(v0, vec2(0, 0), normal),
							Vertex(v1, vec2(1, 0), normal),
							Vertex(v2, vec2(0, 1), normal)};
	return WeakSuperMesh(nodes, {{0, 1, 2}}, id);
}

SmoothParametricCurve VivaniCurve(float r, float eps) {
    return SmoothParametricCurve(
                [r](float t) {return vec3(r*(1+cos(t)), r*sin(t), 2*r*sin(t/2)); },
                [r](float t) {return vec3(-r*sin(t), r*cos(t), r*cos(t/2)); },
                [r](float t) {return vec3(-r*cos(t), -r*sin(t), -r/2*sin(t/2)); }, 420,
                -TAU, TAU, true, eps);
}

SmoothParametricPlaneCurve LissajousCurve(float a, float b, float delta, float r1, float r2, float eps) {
    return SmoothParametricPlaneCurve(
                [r1, r2, a, b, delta](float t) {return vec2(r1*sin(a*t+delta), r2*cos(b*t)); },
                [r1, r2, a, b, delta](float t) {return vec2(a*r1*cos(a*t+delta), -b*r2*sin(b*t)); },
                [r1, r2, a, b, delta](float t) {return vec2(-a*a*r1*sin(a*t+delta), -b*b*r2*cos(b*t)); },
                0, TAU, true, eps);
}

SuperCurve circle(float r, std::function<float(float)> w, const std::function<MaterialPhongConstColor(float)> &mat, int n, vec3 center, vec3 v1, vec3 v2, float eps) {
    return SuperCurve(circle(r, center, v1, v2, eps), w, mat, n, 0, TAU, true);
}



SmoothParametricCurve sphericalSpiral(float a, float t_max, PolyGroupID id, float eps) {
    return SmoothParametricCurve([a](float t) {return vec3(cos(t), sin(t), -a*t)/sqrt(1+a*a*t*t); }, id,  -t_max, t_max, false, eps);
}

SmoothParametricCurve sphericalSpiral(float a, float r, float t_max, PolyGroupID id, float eps) {
    return SmoothParametricCurve([a, r](float t) {return vec3(cos(t), sin(t), -a*t)*r/sqrt(1+a*a*t*t); },id,  -t_max, t_max, false, eps);
}

WeakSuperMesh singleTrig(vec3 v0, vec3 v1, vec3 v2, MaterialPhongConstColor &material, PolyGroupID id) {
    vec3 n = normalize(cross(v1 - v0, v2 - v0));
    vector nodes = {Vertex(v0, vec2(0, 0), n, BLACK, material),
                    Vertex(v1, vec2(0, 1), n, BLACK, material),
                    Vertex(v2, vec2(1, 1), n, BLACK, material)};
    return WeakSuperMesh(nodes, {ivec3(0, 1, 2)}, id);
}


WeakSuperMesh singleTrig(vec3 v0, vec3 v1, vec3 v2, MaterialPhongConstColor &material1, MaterialPhongConstColor &material2, MaterialPhongConstColor &material3, PolyGroupID id) {
    vec3 n = normalize(cross(v1 - v0, v2 - v0));
    vector nodes = {Vertex(v0, vec2(0, 0), n, BLACK, material1),
                    Vertex(v1, vec2(0, 1), n, BLACK, material2),
                    Vertex(v2, vec2(1, 1), n, BLACK, material3)};
    return WeakSuperMesh(nodes, {ivec3(0, 1, 2)}, id);
}
WeakSuperMesh singleQuadShadeSmooth(vec3 outer1, vec3 inner1, vec3 inner2, vec3 outer2, MaterialPhongConstColor &material,
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

WeakSuperMesh singleQuadShadeFlat(vec3 outer1, vec3 inner1, vec3 inner2, vec3 outer2, MaterialPhongConstColor &material,
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

WeakSuperMesh singleQuadShadeFlat(vec3 outer1, vec3 inner1, vec3 inner2, vec3 outer2, MaterialPhongConstColor &material1,
                                  MaterialPhongConstColor &material2, std::variant<int, std::string> id) {
    vec3 n1 = normalize(cross(outer1 - inner1, outer1 - inner2));
    vec3 n2 = normalize(cross(outer2 - inner2, outer2 - inner1));

    vector nodes = {Vertex(outer1, vec2(0, 0), n1, BLACK, material1), Vertex(inner1, vec2(0, 1), n1, BLACK, material1),
                    Vertex(inner2, vec2(1, 1), n1, BLACK, material1), Vertex(inner1, vec2(0, 1), n2, BLACK, material2),
                    Vertex(inner2, vec2(1, 1), n2, BLACK, material2), Vertex(outer2, vec2(1, 1), n2, BLACK, material2)};
    return WeakSuperMesh(nodes, {ivec3(0, 1, 2), ivec3(3, 4, 5)}, id);
}



SmoothParametricSurface sphere(float r, vec3 center, float cutdown, float eps) {
    return SmoothParametricSurface([r, center](float t, float u) {
        return center + vec3(cos(t) * sin(u), sin(t) * sin(u), cos(u))*r;
    }, vec2(0, TAU), vec2(cutdown, PI-0.01), true, false, .01);
}

SmoothParametricSurface DupinCyclide(float a, float b, float d, float eps) {
	return ellipse(a, b, PLANE_ORIGIN, eps).embedding().canal([a, b, d](float t){return d - cos(t)*sqrt(abs(a*a-b*b));});
}

SmoothParametricSurface disk(float r, vec3 center, vec3 v1, vec3 v2,float eps) {
	auto w1 = normalise(v1);
	auto w2 = normalise(v2);
	return SmoothParametricSurface([r, center, w1, w2](float u, float v) { return center + (w1*cos(u) + w2*sin(u))*r; },  vec2(0, TAU), vec2(0, 1), false,true,  eps);
}

SmoothParametricSurface cylinder(float r, vec3 c1, vec3 c2, vec3 v1, vec3 v2, float eps) {
	return  SmoothParametricSurface([r, c1, c2, v1, v2](float u, float v) { return c1 + (cos(u)*v1 + sin(u)*v2)*r + v*(c2 - c1); }, vec2(0, TAU), vec2(0, 1), true, false, eps);
}

SmoothParametricSurface hyperbolic_helicoid(float a, float eps) {
	return SmoothParametricSurface([a](float t, float u) {
		return vec3(std::sin(a*t)* std::sin(u)* std::exp(t),
					std::cos(a*t)* std::sin(u)* std::exp(t),
					std::cos(u)* std::exp(t));

	}, vec2(0, TAU), vec2(0, TAU), false, true, eps);
}

SmoothParametricSurface LawsonTwist(float alpha, Quaternion q, vec2 range_u, float eps) {
	auto L = [alpha](float t, float u) {
		return vec4(::cos(u)*cos(t),
					cos(t)*sin(u),
					sin(t)*cos(alpha*u),
					sin(t)*sin(alpha*u)); };
	return SmoothParametricSurface([L, q](float t, float u) {
									   return stereoProjection(q.normalise().rotateS3(L(t, u))); },
								   vec2(0, TAU), range_u, false, false, eps);
}

SmoothParametricSurface coolLawson(float eps) {
	return LawsonTwist(2, Quaternion(.5, .5, .5, .5), vec2(-PI/4, PI/4), eps);
}

SmoothParametricSurface sudaneseMobius(float eps) {
	return LawsonTwist(.5, Quaternion(.5, .5, .5, .5), vec2(-PI/2, PI/2), eps);
}

SmoothParametricSurface twistedTorus(float a, float m, float n, int dommul1, int dommul2, float eps) {
	return SmoothParametricSurface([a, m, n](float u, float v) {
		return vec3((a + std::cos(n * u / 2) * std::sin(v) - std::sin(n * u / 2) * std::sin(2 * v)) * std::cos(m * u / 2),
			(a + std::cos(n * u / 2) * std::sin(v) - std::sin(n * u / 2) * std::sin(2 * v)) * std::sin(m * u / 2),
			std::sin(n * u / 2) * std::sin(v) + std::cos(n * u / 2) * std::sin(2 * v));
	}, vec2(0, TAU*dommul1), vec2(0, TAU*dommul2), true, true, eps);
}

inline WeakSuperMesh icosahedron(float r, vec3 center, std::variant<int, std::string> id) {
    float phi = (1.f + sqrt(5)) / 2;
    vector verts = {
        Vertex(vec3(1, phi, 0),  vec2(1, 0),  normalise(vec3(1, phi, 0)),     BLACK     ),
        Vertex(vec3(-1, phi, 0), vec2(-1, 0), normalise(vec3(-1, phi, 0)),   BLACK   ),
        Vertex(vec3(1, -phi, 0), vec2(1, 0), normalise(vec3(1, -phi, 0)),   BLACK   ),
        Vertex(vec3(-1, -phi, 0),vec2(-1, 0), normalise(vec3(-1, -phi, 0)), BLACK ),
        Vertex(vec3(0, 1, phi),  vec2(0, phi), normalise(vec3(0, 1, phi)),     BLACK     ),
        Vertex(vec3(0, -1, phi), vec2(0, phi), normalise(vec3(0, -1, phi)),   BLACK   ),
        Vertex(vec3(0, 1, -phi), vec2(0, -phi), normalise(vec3(0, 1, -phi)),   BLACK   ),
        Vertex(vec3(0, -1, -phi),vec2(0, -phi), normalise(vec3(0, -1, -phi)), BLACK ),
        Vertex(vec3(phi, 0, 1),  vec2(phi, 1), normalise(vec3(phi, 0, 1)),     BLACK     ),
        Vertex(vec3(-phi, 0, 1), vec2(-phi, 1), normalise(vec3(-phi, 0, 1)),   BLACK   ),
        Vertex(vec3(phi, 0, -1), vec2(phi, -1), normalise(vec3(phi, 0, -1)),   BLACK   ),
        Vertex(vec3(-phi, 0, -1),vec2(-phi, -1), normalise(vec3(-phi, 0, -1)), BLACK )
        };
    vector<ivec3> faceInds = {};
    for (int i = 0; i < verts.size()-2; i++) {
        for (int j = i+1; j < verts.size()-1; j++) {
            for (int k = j+1; k < verts.size(); k++) {
                vec3 a = verts[i].getPosition();
                vec3 b = verts[j].getPosition();
                vec3 c = verts[k].getPosition();
                if (norm(a-b) < 2.1 && norm(b-c) < 2.1 && norm(c-a) < 2.1) {
                    faceInds.push_back(ivec3(i, j, k));
                }
            }
        }
    }
    for (int i = 0; i < verts.size(); i++) {
        verts[i].setPosition(normalise(verts[i].getPosition())*r + center);
        verts[i].setUV(vec2(verts[i].getPosition().y, verts[i].getPosition().z));
    }
    return WeakSuperMesh(verts, faceInds, id);
}

WeakSuperMesh icosphere(float r, int n, vec3 center, PolyGroupID id, vec4 color) {
    WeakSuperMesh unitSphere = icosahedron(1, vec3(0, 0, 0), id);
    while (n > 0) {
        unitSphere = unitSphere.subdivideEdgecentric(id);
        n--;
    }

    auto normaliseVertices = [r, center, color](BufferedVertex &v) {
        v.setNormal(normalise(v.getPosition()));
    	v.setUV(vec2(v.getPosition().y - v.getPosition().z, v.getPosition().x+ v.getPosition().z));
        v.setPosition(normalise(v.getPosition())*r + center);
    	v.setColor(color);
    };
    unitSphere.deformPerVertex(id, normaliseVertices);

    return unitSphere;
}

WeakSuperMesh disk3d(float r, vec3 center, vec3 v1, vec3 v2, int radial_res, int vertical_res, const std::variant<int, std::string> &id) {
    auto vert = vector<Vertex>();
    vector<ivec3> inds = {};
    vec3 n = normalise(cross(v1, v2));

    vert.emplace_back(center, vec2(0, 0), n, vec4(0, 0, 0, 0));
    for (int h = 1; h <= vertical_res; h++)
        for (int i = 0; i < radial_res; i++)
            {
                float theta = TAU * i / radial_res;
                vert.emplace_back( center + v1 * (cos(theta) / vertical_res * r * h)  + v2 * (sin(theta) / vertical_res * r * h),
                    vec2(i * 1.f / radial_res, h * 1.f / vertical_res), n, vec4(theta, 1.0 * h /vertical_res, 0, 0));
            }
    for (int i = 0; i < radial_res; i++)
        inds.emplace_back(0, i + 1, (i + 1) % radial_res + 1);
    for (int h = 1; h <= vertical_res - 1; h++)
        for (int i = 0; i < radial_res; i++)
        {
            inds.emplace_back(i + 1 + (h - 1) * radial_res, i + 1 + h * radial_res, (i + 1) % radial_res + 1 + h * radial_res);
            inds.emplace_back(i + 1 + (h - 1) * radial_res, (i + 1) % radial_res + 1 + h * radial_res, (i + 1) % radial_res + 1 + (h - 1) * radial_res);
        }
    return WeakSuperMesh(vert, inds, id);
}

WeakSuperMesh singleQuadShadeFlat(vec3 outer1, vec3 inner1, vec3 inner2, vec3 outer2, std::variant<int, std::string> id) {
	vec3 normal          = normalize(cross(outer1 - outer2, outer1 - inner1));
	vector<Vertex> nodes = {Vertex(inner1, vec2(0, 0), normal),
							Vertex(inner2, vec2(1, 0), normal),
							Vertex(outer1, vec2(0, 1), normal),
							Vertex(outer2, vec2(1, 1), normal)};
	return WeakSuperMesh(nodes, {{0, 1, 2}, {3, 1, 2}}, id);
}

Disk3D::Disk3D(const std::vector<Vertex> &nodes, const std::vector<ivec3> &faceInds, vec3 center, vec3 forward, vec3 down, std::variant<int, std::string> id) :
    WeakSuperMesh(nodes, faceInds, id), center(center), forward(forward), down(down), normal(normalise(cross(forward, down))), radius(0) , id(id){

    setEmpiricalRadius();
    setColorInfo();
    }

Disk3D::Disk3D(const char *filename, vec3 center, vec3 forward, vec3 down, std::variant<int, std::string> id)
    : WeakSuperMesh(filename, id), center(center), forward(forward), down(down), normal(normalise(cross(forward, down))), radius(0), id(id) {

    setEmpiricalRadius();
    setColorInfo();

}

Disk3D::Disk3D(float r, vec3 center, vec3 forward, vec3 down, int radial_res, int vertical_res, const std::variant<int, std::string> &id) : WeakSuperMesh(disk3d(r, center, forward, down, radial_res, vertical_res, id)) {
    this->center = center;
    this->forward = forward;
    this->down = down;
    this->normal = normalise(cross(forward, down));
    radius = r;
    this->id = id;
}

void Disk3D::move(vec3 center, vec3 forward, vec3 down, bool scaleWidth) {
    vec3 delta = center - this->center;

    vec3 n = normalise(cross(forward, down));
    for (BufferedVertex &v: getBufferedVertices(id)) {
        vec3 normalComponent = scaleWidth ? normal*widthNormalised(v)*radius : normal*width(v);
        v.setPosition(center + forward*cos(angle(v))*rParam(v)*radius + down*sin(angle(v))*rParam(v)*radius +  normalComponent);
        v.setNormal(normalise(dot(v.getNormal(), this->normal)*n + dot(v.getNormal(), this->forward)*forward + dot(v.getNormal(), this->down)*down));
    }
    this->center = center;
    this->forward = forward;
    this->down = down;
    this->normal = n;
}

float Disk3D::moveRotate(vec3 center, vec3 forward, vec3 down) {
    float distance = norm(center + down- this->center-this->down);
    float angle = distance / radius - asin(dot(normal, cross(this->down, down)));
    for (BufferedVertex &v: getBufferedVertices(id))
        v.setColor(v.getColor() + vec4(angle, 0, 0, 0));

    move(center, forward, down, false);
    return angle;
}

void Disk3D::rotate(float angle) {
    for (BufferedVertex &v: getBufferedVertices(id))
        v.setColor(v.getColor() + vec4(angle, 0, 0, 0));

    move(center, forward, down, false);
}

float Disk3D::rReal(const BufferedVertex &v) {
    return norm(projectVectorToPlane(v.getPosition() - this->center, this->normal));
}

void Disk3D::scaleR(float r, bool scaleWidth) {
    vec3 scaleFactors = scaleWidth ? vec3(r/this->radius) : vec3(r/this->radius, r/this->radius, 1);
    radius = r;
    SpaceAutomorphism scaling = SpaceAutomorphism::scaling(scaleFactors).applyWithBasis(forward, down, normal).applyWithShift(center);
    for (BufferedVertex &v: getBufferedVertices(id)) {
        v.applyFunction(scaling);
        if (scaleWidth) setAbsoluteWidth(v, dot(v.getPosition() - center, normal));
        else setRelativeWidth(v, dot(v.getPosition() - center, normal)/radius);
    }
}

void Disk3D::setR(float r) {
    radius = r;
}

void Disk3D::setEmpiricalRadius() {
    for ( auto &v: getBufferedVertices(id))
        radius = std::max(radius, norm(projectVectorToPlane(v.getPosition() - this->center, this->normal)));
}

void Disk3D::setColorInfo() {
    for ( auto &v: getBufferedVertices(id))
        v.setColor(vec4(   atan2(dot(v.getPosition() - center, forward), dot(v.getPosition() - center, down)),
                                rReal(v)/radius,
                                dot(v.getPosition() - center, normal),
                                dot(v.getPosition() - center, normal)/radius));
}

SmoothParametricSurface cone(const SmoothParametricCurve &base, vec3 apex, float eps) {
	return SmoothParametricSurface([base, apex](float u, float v) {return base(u) + v*(apex - base(u)); }, [base, apex](float u, float v) {return cross(base.tangent(u), apex - base(u)); }, [base, apex](float u, float v) {return cross(base.tangent(u), apex - base(u)); }, base.bounds(), vec2(0, 1), base.isPeriodic(), false, eps);
}



SmoothParametricSurface coneSide(float r, float h, vec3 center, vec3 v1, vec3 v2, vec3 dir, float eps) {
	return cone(circle(r, center, v1, v2, eps), center + normalise(dir)*h, eps);
}







SmoothImplicitSurface sphereImplicit(float r, vec3 center, float eps) {
	return SmoothImplicitSurface([r, center](vec3 p) { return norm2(p - center) - r*r; }, eps);
}

SmoothImplicitSurface torusImplicit(float r, float R, vec3 center, float eps) {
	return SmoothImplicitSurface([r, R, center](vec3 p) {
		vec3 q = p - center;
		return ::pow(norm2(q) + R*R - r*r, 2) - 4*R*R*(q.x*q.x + q.y*q.y);
	}, eps);
}

SmoothImplicitSurface genus2Implicit(float eps) {
	return SmoothImplicitSurface([](vec3 p) {
		float x = p.x, y = p.y, z = p.z;
		return 2.f*y*(y*y - 3.f*x*x)*(1.f - z*z) + ::pow((x*x + y*y),2) - (9.f*z*z - 1.f)*(1.f - z*z);
	}, eps);
}

SmoothImplicitSurface wineGlass(float eps) {
	return SmoothImplicitSurface([](vec3 p) {
		return p.x*p.x + p.y*p.y - ::pow(log(p.z +3.2), 2) - .02;
	}, eps);
}

SmoothImplicitSurface equipotentialSurface(vector<vec3> points, vector<float> charges, float potential, float eps) {
	return SmoothImplicitSurface([points, charges, potential](vec3 p) {
		float res = 0;
		for (int i = 0; i < points.size(); i++)
			res += charges[i]/norm(p - points[i]);
		return res - potential;
	}, eps);
}

SmoothImplicitSurface chair(float k, float a, float b, float eps) {
	return SmoothImplicitSurface([k, a, b](vec3 p) {
		return square(norm2(p) - a*k*k) - b*(square(p.z - k) - 2.f*p.x*p.x)*(square(p.z + k) - 2.f*square(p.y));
	}, eps);
}

SmoothImplicitSurface tangleCube(float eps) {
	return SmoothImplicitSurface(
		[](float x, float y, float z) { return pow4(x) - square(x)*5 + pow4(y) - square(y)*5 + pow4(z) - square(z)*5 + 10; },
		[](float x, float y, float z) { return vec3(cube(x)*4 - x*10,  cube(y)*4 - y*10, cube(z)*4 - z*10); },
		eps);
}
SmoothImplicitSurface wineImplicit(float eps) {
	return SmoothImplicitSurface(
		[](float x, float y, float z) { return -.4*sin(x*5) -.4*sin(y*5) -.4*sin(z*5) + .1*square(x) + .3*square(y) + .2*square(z) - .5; },
		[](float x, float y, float z) { return vec3( -2*cos(x*5)+ .2*x, -2*cos(y*5)+ .6*y, -2*cos(z*5)+ .4*z ); },
		eps);
}

SmoothImplicitSurface gumdrop(float eps) {
	//return 4*(x**4 + (y**2 + z**2)**2) + 17*x**2*(y**2 + z**2) - 20*(x**2 + y**2 + z**2) + 17
	return SmoothImplicitSurface([](float x, float y, float z) {
									 return 4*pow4(x) + 4*pow4(y) + 4*pow4(z) + 8*sq(y)*sq(z) + 17*sq(x)*sq(y) + 17*sq(x)*sq(z) - 20*sq(x) - 20*sq(y) - 20*sq(z) + 17;
								 },
								 [](float x, float y, float z) { return vec3(
										 16*cube(x) + 34*x*sq(z) + 34*x*sq(y) - 40*x,
										 16*cube(y) + 16*y*sq(z) + 34*y*sq(x) - 40*y,
										 16*cube(z) + 16*z*sq(y) + 34*z*sq(x) - 40*z); }, eps);
}

SmoothImplicitSurface genus2Implicit2(float eps) {
	return SmoothImplicitSurface([](float x, float y, float z) {
									 return .04 - pow4(x) + 2*pow6(x) - pow8(x) + 2*sq(x)*sq(y) - 2*pow4(x)*sq(y) - pow4(y) - sq(z);
								 },
								 [](float x, float y, float z) { return vec3(
										 -4*pow3(x) + 12*pow5(x) - 8*pow7(x) + 8*x*sq(y) - 8*pow3(x)*sq(y),
										 4*y*sq(x) - 4*y*pow4(x) - 4*pow3(y),
										 -2*z); }, eps);

}



SmoothImplicitSurface genus2Implicit3(float c, float d, float eps) {
	return SmoothImplicitSurface([c, d](float x, float y, float z) {
		return sq(sq(x-1) + sq(y) - sq(c))*sq(sq(x+1) + sq(y) - sq(c)) + sq(z) - d;
	}, eps);

}

SmoothImplicitSurface genus2Implicit4(float d, float eps) {
	return SmoothImplicitSurface([d](float x, float y, float z) {
		return sq(sq(x) - pow4(x) - sq(y)) + sq(z) - d;
	}, eps);

}

SmoothImplicitSurface superellipsoid(float alpha1, float alpha2, float a, float b, float c, float r, float eps) {
	return SmoothImplicitSurface([alpha1, alpha2, a, b, c, r](float x, float y, float z) {
		return ::pow(::pow(abs(x/a), 2/alpha2) + ::pow(abs(y/b), 2/alpha2), alpha2/alpha1) + ::pow(abs(z/c), 2/alpha1) - r;
	}, eps);

}

SmoothImplicitSurface superQuadric(float alpha, float beta, float gamma, float a, float b, float c, float r, float eps) {
	return SmoothImplicitSurface([alpha, beta, gamma, a, b, c, r](float x, float y, float z) {
		return ::pow(abs(x/a), alpha) + ::pow(abs(y/b), beta) + ::pow(abs(z/c), gamma) - r;
	}, eps);

}

SmoothImplicitSurface K3Surface222(float eps) {
	return SmoothImplicitSurface([](float x, float y, float z) {
		return (sq(x) + 1)*(sq(y) + 1)*(sq(z) + 1) + 8*x*y*z - 2;
	}, eps);

}

WeakSuperMesh arrow(vec3 start, vec3 head, float radius, float head_len, float head_radius, int radial, int straight, float eps, std::variant<int, std::string> id) {
	vec3 direction     = normalise(head - start);
	auto w  = orthogonalComplementBasis(direction);
	vec3 w1  = w.first;
	vec3 w2  = w.second;
	auto id1 = randomID();
	WeakSuperMesh mesh = WeakSuperMesh(WeakSuperMesh(cylinder(radius, start, head, w1, w2, eps), radial, straight, id1));
	mesh.deformPerVertex(id1, [start, head](BufferedVertex &v) {
		v.setColor(vec4(start, head.z));
		v.setUV(vec2(head.x, head.y));
	});
	id1 = randomID();
	mesh.addUniformSurface(coneSide(head_radius, head_len, head, w1, w2, direction, eps), radial, straight, id1);
	mesh.deformPerVertex(id1, [start, head](BufferedVertex &v) {
		v.setColor(vec4(start, head.z));
		v.setUV(vec2(head.x, head.y));
	});
	id1 = randomID();
	mesh.addUniformSurface( disk(head_radius, head, w1, w2, eps), radial,  straight, id1);
	mesh.deformPerVertex(id1, [start, head](BufferedVertex &v) {
		v.setColor(vec4(start, head.z));
		v.setUV(vec2(head.x, head.y));
	});
	return mesh;
}

WeakSuperMesh drawArrows(const vector<vec3> &points, const vector<vec3> &directions, float radius, float head_len, float head_radius, int radial, int straight, float eps,
						 const std::variant<int, std::string> &id) {
	WeakSuperMesh mesh;
	for (int i = 0; i < points.size(); i++)
		mesh.merge(arrow(points[i], directions[i], radius, head_len, head_radius, radial, straight, eps,  randomID()));
	return mesh;
}

WeakSuperMesh drawArrows(const vector<vec3> &points, const VectorFieldR3 &field, const std::function<float(float)>& radius,const std::function<float(float)>& len, const std::function<float(float)> &head_len, const std::function<float(float)> &head_radius,
		int radial, int straight, float eps, const std::variant<int, std::string> &id) {
	WeakSuperMesh mesh;
	for (int i = 0; i < points.size(); i++) {
		float size = length(field(points[i]));
		mesh.merge(arrow(points[i], points[i]+field(points[i])*len(size), radius(size), head_len(size), head_radius(size), radial, straight, eps, randomID()));
	}
	return mesh;
}

WeakSuperMesh drawVectorFieldArrows(const VectorFieldR3 &field, const vector<vec3> &points, const std::function<float(float)>& len, const std::function<float(float)> &radius, const std::function<float(float)>& head_len,
		const std::function<float(float)>& head_radius, int radial, int straight, float eps, const std::variant<int, std::string> &id) {
	WeakSuperMesh mesh;
	for (auto point : points) {
		float s = norm(field(point));
		mesh.merge(arrow(point, point + normalise(field(point))*len(s), radius(s), head_len(s), head_radius(s), radial, straight, eps, randomID()));
	}
	return mesh;
}

vec3 getArrayStart(const BufferedVertex &v) {
	return vec3(v.getColor().x, v.getColor().y, v.getColor().z);
}

vec3 getArrayHead(const BufferedVertex &v) {
	return vec3(v.getUV().x, v.getUV().y, v.getColor().w);
}

vec3 getArrayDirection(const BufferedVertex &v) {
	return getArrayHead(v) - getArrayStart(v);
}

vec3 getArrayHead(const std::variant<int, std::string>& id, WeakSuperMesh &mesh) {return getArrayHead(mesh.getAnyVertexFromPolyGroup(id));}
vec3 getArrayStart(const std::variant<int, std::string>& id, WeakSuperMesh &mesh) {return getArrayStart(mesh.getAnyVertexFromPolyGroup(id));}
vec3 getArrayDirection(const std::variant<int, std::string>& id, WeakSuperMesh &mesh) {return getArrayDirection(mesh.getAnyVertexFromPolyGroup(id));}


SmoothParametricCurve PLCurve(std::vector<vec3> points) {
	return SmoothParametricCurve([points](float t) {
		if (t<0) return points[0];
		if (t>=points.size()-1) return points[points.size()-1];
		int i = static_cast<int>(t);
		float t1 = t - i;
		return t1*points[i+1] + (1-t1)*points[i];
	}, "PL", 0, points.size()-1, points[0]==points[points.size()-1], .01);
}

SmoothParametricPlaneCurve GernsterWave(float a, float b, float k, float c) {
	return SmoothParametricPlaneCurve([a, b, k, c](float t) {
		return vec2(a + ::exp(b*k)/k* ::sin(a*k + t*k*c),
					b- ::exp(b*k)/k* ::cos(a*k + t*k*c)); }, 0, 1, false, .01);
}

VectorFieldR3 PousevillePlanarFlow(float h, float nabla_p, float mu, float v0, float eps) {
	return VectorFieldR3([h, nabla_p, mu, v0](vec3 p) {return vec3(0,nabla_p/(2*mu)*p.z*(h-p.z) + v0*p.z/h, 0); }, eps);
}

VectorFieldR3 PousevillePipeFlow(float nabla_p, float mu, float c1, float c2, float eps) {
	return VectorFieldR3([nabla_p, mu, c1, c2](vec3 p) {
		float r = norm(p);
		return vec3(-nabla_p*r*r/(4*mu) + c1* ::log(r) + c2); }, eps);
}

vec3 freeSurfaceComponent(float amplitude, float phase, vec2 waveNumber, float angFrequency, float h, vec2 ab, float t) {
	float k     = norm(waveNumber);
	float kx    = waveNumber.x;
	float ky    = waveNumber.y;
	float theta = kx*ab.x + ky*ab.y - angFrequency*t - phase;

	return vec3(-kx/k*amplitude/tanh(k*h)*sin(theta),
				-ky/k*amplitude/tanh(k*h)*sin(theta),
				amplitude*cos(theta));
}

SurfaceParametricPencil freeSurface(vector<float> a_m, vector<float> phi_m, vector<vec2> k_m, vector<float> omega_m, float h) {
	return SurfaceParametricPencil([a_m, phi_m, k_m, omega_m, h](float t, vec2 tu) {
		vec3 sum = vec3(tu.x, tu.y, 0);
		for (int i = 0; i < a_m.size(); i++)
			sum += freeSurfaceComponent(a_m[i], phi_m[i], k_m[i], omega_m[i], h, tu, t);
		return sum;
	}, vec2(-0, PI), vec2(-0, PI));
}

WeakSuperMesh particles(int n, vec3 bound1, vec3 bound2, float radius) {
	WeakSuperMesh mesh = WeakSuperMesh();
	for (int i = 0; i < n; i++) {
		vec3 pos  = randomUniform(bound1, bound2);
		auto part = icosphere(radius, 1, pos, randomID(), vec4(pos, 0));
		mesh.merge(part);
	}
	return mesh;
}

WeakSuperMesh box(vec3 size, vec3 center, std::variant<int, std::string> id) {
	float dx        = size.x/2;
	float dy        = size.y/2;
	float dz        = size.z/2;
	WeakSuperMesh b = singleQuadShadeFlat(center+vec3(dx, -dy, -dz), center+vec3(dx, -dy, dz), center+vec3(dx, dy, dz), center+vec3(dx, dy, -dz), id);
	b.merge(singleQuadShadeFlat(center+vec3(-dx, -dy, -dz), center+vec3(-dx, -dy, dz), center+vec3(-dx, dy, dz), center+vec3(-dx, dy, -dz), id));
	b.merge(singleQuadShadeFlat(center+vec3(-dx, dy, -dz), center+vec3(-dx, dy, dz), center+vec3(dx, dy, dz), center+vec3(dx, dy, -dz), id));
	b.merge(singleQuadShadeFlat(center+vec3(-dx, -dy, -dz), center+vec3(-dx, -dy, dz), center+vec3(dx, -dy, dz), center+vec3(dx, -dy, -dz), id));
	b.merge(singleQuadShadeFlat(center+vec3(-dx, -dy, dz), center+vec3(-dx, dy, dz), center+vec3(dx, dy, dz), center+vec3(dx, -dy, dz), id));
	b.merge(singleQuadShadeFlat(center+vec3(-dx, -dy, -dz), center+vec3(-dx, dy, -dz), center+vec3(dx, dy, -dz), center+vec3(dx, -dy, -dz), id));

	return b;
}

SmoothParametricCurve trefoil(float r, float R, float eps) {
	return SmoothParametricCurve([r, R](float t) {return
		vec3(r*cos(t) - R*cos(2*t),
			r*sin(t) + R*sin(2*t),
			R*sin(3*t)); }, "Trefoil", 0, TAU, true, eps);
}
 SmoothParametricCurve torusKnot23(float scale, float R, float eps) {
	return SmoothParametricCurve([scale, R](float t)
		{return vec3((R + cos(3 * t)) * cos(2 * t),
			(R + cos(3 * t)) * sin(2 * t),
			sin(3 * t))*scale; }, "TorusKnot23", 0, TAU, true, eps);
}
