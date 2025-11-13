
#include "specific.hpp"


Permutation Permutation::operator&(const Permutation &p) const {
	if (size() != p.size())
		throw ValueError("Permutations must have the same size to be composed.", __FILE__, __LINE__);
	auto lst = p.perm;
	for (int i = 0; i < lst.size(); i++)
		lst[i] = perm[lst[i]];
	return Permutation(lst);
}

Permutation Permutation::operator~() const {
	vector<int> lst = vector<int>(perm.size());
	for (int i = 0; i < perm.size(); i++)
		lst[perm[i]] = i;
	return Permutation(lst);
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
    return circle(r, ORIGIN_R2, eps).embedding(v1, v2, center);
}

IndexedMesh3D singleTrig(vec3 v0, vec3 v1, vec3 v2, PolyGroupID id) {
	vec3 normal          = normalize(cross(v1 - v0, v2 - v0));
	vector<Vertex> nodes = {Vertex(v0, vec2(0, 0), normal, vec4(0)),
							Vertex(v1, vec2(1, 0), normal, vec4(0)),
							Vertex(v2, vec2(0, 1), normal, vec4(0))};
	return IndexedMesh3D(nodes, {{0, 1, 2}}, id);
}

IndexedMesh3D singleQuadShadeSmooth(vec3 outer1, vec3 inner1, vec3 inner2, vec3 outer2, PolyGroupID id) {
	auto n1 = normalize(cross(inner1 - outer1, inner2 - outer1));
	auto n2 = normalize(normalize(cross(inner2 - inner1, outer1 - inner1)) + normalize(cross(inner2 - inner1, outer2 - inner1)));
	auto n3 = normalize(normalize(cross(inner1 - inner2, outer1 - inner2)) + normalize(cross(inner1 - inner2, outer2 - inner2)));
	auto n4 = normalize(cross(inner2 - outer2, inner1 - outer2));
	if (dot(n1, n2) < 0) n2 = -n2;
	if (dot(n1, n3) < 0) n3 = -n3;
	if (dot(n1, n4) < 0) n4 = -n4;
	vector<Vertex> verts = {
	Vertex(outer1, vec2(0, 0), n1, BLACK),
	Vertex(inner1, vec2(1, 0), n2, BLACK),
	Vertex(inner2, vec2(0, 1), n3, BLACK),
	Vertex(outer2, vec2(1, 1), n4, BLACK),
	};
	vector<ivec3> tris = {
	ivec3(0, 1, 2),
	ivec3(0, 2, 3),
	};
	return IndexedMesh3D(verts, tris, id);
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


SmoothParametricCurve sphericalSpiral(float a, float t_max, PolyGroupID id, float eps) {
    return SmoothParametricCurve([a](float t) {return vec3(cos(t), sin(t), -a*t)/sqrt(1+a*a*t*t); }, id,  -t_max, t_max, false, eps);
}

SmoothParametricCurve sphericalSpiral(float a, float r, float t_max, PolyGroupID id, float eps) {
    return SmoothParametricCurve([a, r](float t) {return vec3(cos(t), sin(t), -a*t)*r/sqrt(1+a*a*t*t); },id,  -t_max, t_max, false, eps);
}



IndexedMesh3D pyramid(const vector<vec3> &cornersDown, vec3 apex, PolyGroupID id) {
	auto m = singleQuadShadeFlat(cornersDown, id);
	m.merge(singleTrig(cornersDown[0], cornersDown[1], apex, id));
	m.merge(singleTrig(cornersDown[1], cornersDown[2], apex, id));
	m.merge(singleTrig(cornersDown[2], cornersDown[3], apex, id));
	m.merge(singleTrig(cornersDown[3], cornersDown[0], apex, id));
	return m;
}

SmoothParametricSurface plane(vec3 corner, vec3 side1, vec3 side2, float eps) {
	return SmoothParametricSurface([corner, side1, side2](float u, float v) {
		return corner + u * side1 + v * side2;
	}, vec2(0, 1), vec2(0, 1), true, false, eps);
}

SmoothParametricSurface sphere(float r, vec3 center, float cutdown, float eps) {
    return SmoothParametricSurface([r, center](float t, float u) {
        return center + vec3(cos(t) * sin(u), sin(t) * sin(u), cos(u))*r;
    }, vec2(0, TAU), vec2(cutdown, PI-0.01), true, false, .01);
}

SmoothParametricSurface ellipsoid(float rx, float ry, float rz, vec3 center, float eps) {
	return SmoothParametricSurface([rx, ry, rz, center](float t, float u) {
		return center - vec3(cos(t) * cos(u)*rx, sin(t) * cos(u)*ry, sin(u)*rz);
	}, vec2(0, TAU), vec2(-PI/2, PI/2), false, false, eps);
}

SmoothParametricSurface DupinCyclide(float Ra, float Rb, float r, float delta_r, float eps) {
	return SmoothParametricSurface([Ra, Rb, r, d=delta_r](float u, float v) {
		return vec3(
			(r*(d - Ra*cos(u)*cos(v)) + pow2(Rb)*cos(u))/(Ra - d*cos(u)*cos(v)),
			(Ra-r*cos(v))*sin(u)*Rb/(Ra - d*cos(u)*cos(v)),
			(d*cos(u)-r)*sin(v)*Rb/(Ra - d*cos(u)*cos(v)));
	},	 vec2(0, TAU), vec2(0, TAU), true, true, eps);
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
		return vec4(cos(u)*cos(t),
					cos(t)*sin(u),
					sin(t)*cos(alpha*u),
					sin(t)*sin(alpha*u)); };
	return SmoothParametricSurface([L, q](float t, float u) {
									   return stereoProjection(q.normalise()*Quaternion(L(t, u))); },
								   vec2(0, TAU), range_u, false, false, eps);
}

SmoothParametricSurface coolLawson(float eps) {
	return LawsonTwist(2, Quaternion(.5, .5, .5, .5), vec2(-PI/4, PI/4), eps);
}

SmoothParametricSurface sudaneseMobius(float eps) {
	return LawsonTwist(.5, Quaternion(.5, .5, .5, .5f), vec2(-PI/2, PI/2), eps);
}


SmoothParametricSurface parametricPlane2ptHull(vec3 v1, vec3 v2, float normal_shift, float eps) {
	return SmoothParametricSurface([v1, v2, normal_shift](float u, float v) { return u*v1 + v*v2 + normalise(cross(v1, v2))*normal_shift; }, vec2(0, 1), vec2(0, 1), false, false, eps);
}

SmoothParametricSurface twistedTorus(float a, float m, float n, int dommul1, int dommul2, float eps) {
	return SmoothParametricSurface([a, m, n](float u, float v) {
		return vec3((a + std::cos(n * u / 2) * std::sin(v) - std::sin(n * u / 2) * std::sin(2 * v)) * std::cos(m * u / 2),
			(a + std::cos(n * u / 2) * std::sin(v) - std::sin(n * u / 2) * std::sin(2 * v)) * std::sin(m * u / 2),
			std::sin(n * u / 2) * std::sin(v) + std::cos(n * u / 2) * std::sin(2 * v));
	}, vec2(0, TAU*dommul1), vec2(0, TAU*dommul2), true, true, eps);
}

inline IndexedMesh3D icosahedron(float r, vec3 center, PolyGroupID id) {
    float phi = (1.f + sqrt(5)) / 2;
    vector verts_ = {
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
	vector<Vertex> verts = {};
    vector<ivec3> faceInds = {};
    for (int i = 0; i < verts_.size()-2; i++) {
        for (int j = i+1; j < verts_.size()-1; j++) {
            for (int k = j+1; k < verts_.size(); k++) {
                vec3 a = verts_[i].get_position();
                vec3 b = verts_[j].get_position();
                vec3 c = verts_[k].get_position();
                if (norm(a-b) < 2.1 && norm(b-c) < 2.1 && norm(c-a) < 2.1) {
                	vec3 normal = normalise(cross(b - a, c - a));
                	if (dot(normal, a) > 0) normal = -normal;
                	verts.emplace_back(a, verts_[i].get_uv(), normal, BLACK);
                	verts.emplace_back(b, verts_[j].get_uv(), normal, BLACK);
                	verts.emplace_back(c, verts_[k].get_uv(), normal, BLACK);
                	int l = verts.size();
                    faceInds.push_back(ivec3(l-3, l-2, l-1));
                }
            }
        }
    }
    for (int i = 0; i < verts.size(); i++) {
        verts[i].set_position(normalise(verts[i].get_position())*r + center);
    }
    return IndexedMesh3D(verts, faceInds, id);
}

IndexedMesh3D icosphere(float r, int n, vec3 center, PolyGroupID id, vec4 color) {
    IndexedMesh3D unitSphere = icosahedron(1, vec3(0, 0, 0), id);
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

IndexedMesh3D disk3d(float r, vec3 center, vec3 v1, vec3 v2, int radial_res, int vertical_res, const PolyGroupID &id) {
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
    return IndexedMesh3D(vert, inds, id);
}

IndexedMesh3D singleQuadShadeFlat(vec3 inner1, vec3 outer1, vec3 inner2, vec3 outer2, PolyGroupID id) {
	return singleQuadShadeFlat(inner1, outer1, inner2, outer2, BLACK, id);
}

IndexedMesh3D singleQuadShadeFlat(vec3 inner1, vec3 outer1, vec3 inner2, vec3 outer2, vec4 color, PolyGroupID id) {
	vec3 normal          = normalize(cross(outer1 - inner1, inner2 - inner1));
	vector<Vertex> nodes = {Vertex(inner1, vec2(0, 0), normal, color),
							Vertex(outer1, vec2(1, 0), normal, color),
							Vertex(inner2, vec2(1, 1), normal, color),
							Vertex(outer2, vec2(0, 1), normal, color)};
	return IndexedMesh3D(nodes, {{0, 1, 2}, {0, 3, 2}}, id);
}

IndexedMesh3D singleQuadShadeFlat(const vector<vec3> &corners, PolyGroupID id) {
	return singleQuadShadeFlat(corners[0], corners[1], corners[2], corners[3], id);
}


SmoothParametricSurface cone(const SmoothParametricCurve &base, vec3 apex, float eps) {
	return SmoothParametricSurface([base, apex](float u, float v) {return base(u) + v*(apex - base(u)); }, [base, apex](float u, float v) {return cross(base.tangent(u), apex - base(u)); }, [base, apex](float u, float v) {return cross(base.tangent(u), apex - base(u)); }, base.bounds(), vec2(0, 1), base.isPeriodic(), false, eps);
}



SmoothParametricSurface coneSide(float r, float h, vec3 center, vec3 v1, vec3 v2, vec3 dir, float eps) {
	return cone(circle(r, center, v1, v2, eps), center + normalise(dir)*h, eps);
}




SmoothParametricSurface torus(float r, float R, vec3 center, float eps) {
	return SmoothParametricSurface([r, R, center](float u, float v) {
		return center + vec3((R + r * ::cos(v)) * ::cos(u), (R + r * ::cos(v)) * ::sin(u), r * ::sin(v));
	}, vec2(0, TAU), vec2(0, TAU), true, true, eps);
}

SmoothImplicitSurface sphereImplicit(float r, vec3 center, float eps) {
	return SmoothImplicitSurface([r, center](vec3 p) { return norm2(p - center) - r*r; }, eps);
}

SmoothImplicitSurface torusImplicit(float r, float R, vec3 center, float eps) {
	return SmoothImplicitSurface([r, R, center](vec3 p) {
		vec3 q = p - center;
		return std::pow(norm2(q) + R*R - r*r, 2) - 4*R*R*(q.x*q.x + q.y*q.y);
	}, eps);
}

SmoothImplicitSurface genus2Implicit(float eps) {
	return SmoothImplicitSurface([](vec3 p) {
		float x = p.x, y = p.y, z = p.z;
		return 2.f*y*(y*y - 3.f*x*x)*(1.f - z*z) + std::pow((x*x + y*y),2) - (9.f*z*z - 1.f)*(1.f - z*z);
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
		return std::pow(::pow(abs(x/a), 2/alpha2) + std::pow(abs(y/b), 2/alpha2), alpha2/alpha1) + std::pow(abs(z/c), 2/alpha1) - r;
	}, eps);

}

SmoothImplicitSurface superQuadric(float alpha, float beta, float gamma, float a, float b, float c, float r, float eps) {
	return SmoothImplicitSurface([alpha, beta, gamma, a, b, c, r](float x, float y, float z) {
		return ::pow(abs(x/a), alpha) + std::pow(abs(y/b), beta) + ::pow(abs(z/c), gamma) - r;
	}, eps);

}

SmoothImplicitSurface K3Surface222(float eps) {
	return SmoothImplicitSurface([](float x, float y, float z) {
		return (sq(x) + 1)*(sq(y) + 1)*(sq(z) + 1) + 8*x*y*z - 2;
	}, eps);

}

SmoothImplicitSurface planeImplicit(vec3 normal, float d, vec3 center, float eps) {
	return SmoothImplicitSurface([normal, d, center](vec3 p) {
		return dot(normal, p - center) - d;
	}, eps);
}

IndexedMesh3D arrow(vec3 start, vec3 head, float radius, float head_len, float head_radius, int radial, int straight, float eps, PolyGroupID id) {
	vec3 direction     = normalise(head - start);
	auto w  = orthogonalComplementBasis(direction);
	vec3 w1  = w.first;
	vec3 w2  = w.second;
	auto id1 = randomID();
	IndexedMesh3D mesh = IndexedMesh3D(IndexedMesh3D(cylinder(radius, start, head, w1, w2, eps), radial, straight, id1));
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

IndexedMesh3D drawArrows(const vector<vec3> &points, const vector<vec3> &directions, float radius, float head_len, float head_radius, int radial, int straight, float eps,
						 const PolyGroupID &id) {
	IndexedMesh3D mesh;
	for (int i = 0; i < points.size(); i++)
		mesh.merge(arrow(points[i], directions[i], radius, head_len, head_radius, radial, straight, eps,  randomID()));
	return mesh;
}

IndexedMesh3D drawArrows(const vector<vec3> &points, const VectorField &field, const std::function<float(float)>& radius,const std::function<float(float)>& len, const std::function<float(float)> &head_len, const std::function<float(float)> &head_radius,
		int radial, int straight, float eps, const PolyGroupID &id) {
	IndexedMesh3D mesh;
	for (int i = 0; i < points.size(); i++) {
		float size = length(field(points[i]));
		mesh.merge(arrow(points[i], points[i]+field(points[i])*len(size), radius(size), head_len(size), head_radius(size), radial, straight, eps, randomID()));
	}
	return mesh;
}

IndexedMesh3D drawVectorFieldArrows(const VectorField &field, const vector<vec3> &points, const std::function<float(float)>& len, const std::function<float(float)> &radius, const std::function<float(float)>& head_len,
		const std::function<float(float)>& head_radius, int radial, int straight, float eps, const PolyGroupID &id) {
	IndexedMesh3D mesh;
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

vec3 getArrayHead(const PolyGroupID& id, IndexedMesh3D &mesh) {return getArrayHead(mesh.getAnyVertexFromPolyGroup(id));}
vec3 getArrayStart(const PolyGroupID& id, IndexedMesh3D &mesh) {return getArrayStart(mesh.getAnyVertexFromPolyGroup(id));}
vec3 getArrayDirection(const PolyGroupID& id, IndexedMesh3D &mesh) {return getArrayDirection(mesh.getAnyVertexFromPolyGroup(id));}


SmoothParametricCurve PLCurve(std::vector<vec3> points) {
	return SmoothParametricCurve([points](float t) {
		if (t<0) return points[0];
		if (t>=points.size()-1) return points[points.size()-1];
		int i = static_cast<int>(t);
		float t1 = t - i;
		return t1*points[i+1] + (1-t1)*points[i];
	}, "PL", 0, points.size()-1, points[0]==points[points.size()-1], .01);
}

SmoothParametricCurve segment(vec3 p0, vec3 p1, float t0, float t1) {
	return SmoothParametricCurve([p0, p1, t0, t1](float t) { return (p1*(t-t0) + p0*(t1-t))/(t1-t0); }, "seg", t0, t1, false, .01);
}

SmoothParametricPlaneCurve GernsterWave(float a, float b, float k, float c) {
	return SmoothParametricPlaneCurve([a, b, k, c](float t) {
		return vec2(a + ::exp(b*k)/k* ::sin(a*k + t*k*c),
					b- ::exp(b*k)/k* ::cos(a*k + t*k*c)); }, 0, 1, false, .01);
}

VectorField PousevillePlanarFlow(float h, float nabla_p, float mu, float v0, float eps) {
	return VectorField([h, nabla_p, mu, v0](vec3 p) {return vec3(0,nabla_p/(2*mu)*p.z*(h-p.z) + v0*p.z/h, 0); }, eps);
}

VectorField PousevillePipeFlow(float nabla_p, float mu, float c1, float c2, float eps) {
	return VectorField([nabla_p, mu, c1, c2](vec3 p) {
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

ImplicitVolume implicitBall(float r, vec3 center, float eps) {
	return ImplicitVolume(RealFunctionR3([r, center](vec3 p) { return pow2(r) - pow2(norm(p - center)); }, eps), center-vec3(r, r, r)*1.05, center+vec3(r, r, r)*1.05);
}

ImplicitVolume implicitVolumeEllipsoid(float rx, float ry, float rz, vec3 center, float eps) {
	return ImplicitVolume(RealFunctionR3([rx, ry, rz, center](vec3 p) {
		auto x = p - center;
		return 1 - pow2(x.x/rx) - pow2(x.y/ry) - pow2(x.z/rz);
	}, eps), center-vec3(rx, ry, rz)*1.5, center+vec3(rx, ry, rz)*1.5);
}

IndexedMesh3D generate_random_particle_mesh(int n, vec3 bound1, vec3 bound2, float radius) {
	IndexedMesh3D mesh = IndexedMesh3D();
	for (int i = 0; i < n; i++) {
		vec3 pos  = randomUniform(bound1, bound2);
		auto part = icosphere(radius, 1, pos, randomID(), vec4(pos, 0));
		mesh.merge(part);
	}
	return mesh;
}

IndexedMesh3D box(vec3 size, vec3 center, vec4 color, PolyGroupID id) {
	vec3 dx = vec3(size.x/2, 0, 0);
	vec3 dy = vec3(0, size.y/2, 0);
	vec3 dz = vec3(0, 0, size.z/2);
	vec3 c = center;
	vector<vec3> cornersUp = {c - dx - dy + dz, c + dx - dy + dz, c + dx + dy + dz, c - dx + dy + dz};
	vector<vec3> cornersDown = {c - dx - dy - dz, c + dx - dy - dz, c + dx + dy - dz, c - dx + dy - dz};
	IndexedMesh3D b = box(cornersUp, cornersDown, color);
	return b;
}

IndexedMesh3D box(vec3 size, vec3 center, PolyGroupID id) {
	return box(size, center, vec4(0, 0, 0, 0), id);
}

IndexedMesh3D box(vector<vec3> cornersUp, vector<vec3> cornersDown, vec4 color) {
	vec3 a000 = cornersUp[0];
	vec3 a001 = cornersUp[1];
	vec3 a010 = cornersUp[3];
	vec3 a011 = cornersUp[2];
	vec3 a100 = cornersDown[0];
	vec3 a101 = cornersDown[1];
	vec3 a110 = cornersDown[3];
	vec3 a111 = cornersDown[2];
	auto m = singleQuadShadeFlat(a000, a001, a011, a010,color,  randomID());
	m.merge(singleQuadShadeFlat(a100, a101, a111, a110, color, randomID()));
	m.merge(singleQuadShadeFlat(a000, a001, a101, a100, color, randomID()));
	m.merge(singleQuadShadeFlat(a010, a011, a111, a110, color, randomID()));
	m.merge(singleQuadShadeFlat(a000, a010, a110, a100, color, randomID()));
	m.merge(singleQuadShadeFlat(a001, a011, a111, a101, color, randomID()));
	return m;
}

IndexedMesh3D paraleblahblapid(vec3 corner, vec3 dir1, vec3 dir2, vec3 dir3) {
	vector cornersUp = {corner, corner + dir1, corner + dir2 + dir1, corner + dir2 , };
	vector cornersDown = {corner + dir3, corner + dir1 + dir3, corner + dir1 + dir2+ dir3, corner + dir2+ dir3,  };
	return box(cornersUp, cornersDown);
}

SmoothParametricCurve trefoil(float r, float R, float eps) {
	return SmoothParametricCurve([r, R](float t) {return
		vec3(r*cos(t) - R*cos(2*t),
			r*sin(t) + R*sin(2*t),
			R*sin(3*t)); }, "Trefoil", 0, TAU, true, eps);
}
 SmoothParametricCurve torusKnot23(float scale, float R, float eps) {
	return SmoothParametricCurve([scale, R](float t)
		{
		float r = R + cos(3 * t);
		return vec3(r*cos(2 * t),
			r*sin(2 * t),
			-sin(3 * t))*scale; }, "TorusKnot23", 0, TAU, true, eps);
}

SmoothParametricCurve torusKnot_pq(int p, int q, float R, float scale, float eps) {
	return SmoothParametricCurve([scale, R, p, q](float t)
			{
			float r = R + cos(q * t);
			return vec3(r*cos(p * t),
				r*sin(p * t),
				-sin(q * t))*scale; }, "TorusKnotpq", 0, TAU, true, eps);}



RealFunction smoothenReLu(float c0, float x_change_to_id) {
	return RealFunction([c0, x_change_to_id](float t) {
		if (t>=x_change_to_id) return t;
		float a = 2.0*c0 - x_change_to_id;
		float b = 2.0*x_change_to_id - 3.0*c0;
		float x = t/x_change_to_id;
		return (a*x + b)*x*x + c0;
	});
}

RealFunction expImpulse(float peak, float decay) {
	auto F = [](float x) {
		if (abs(x) >= 1) return 0.f;
		return std::exp(1.f-1.f/(1.f-x*x));
	};
	return RealFunction([F, peak, decay](float t) {
		return F((t-peak)/decay);
	});
}

RealFunction expImpulse_k(float peak, float decay, float k) {
	auto F = [k](float x) {
		if (abs(x) >= 1) return 0.f;
		return 1/(std::exp(k*((1-abs(2*x))/(x*x-abs(x))))+1);
	};
	return RealFunction([F, peak, decay](float t) {
		return F((t-peak)/decay);
	});
}

RealFunction cubicShroom(float center, float margin) {
	return RealFunction([center, margin](float t) {
			if (abs(t-center) > margin) return 0.f;
			return 1.f - sq(t-center)*(3.f-2.f*abs(t-center)/margin)/sq(margin);
		});
}

RealFunction powerShroom(float begin, float end, float zeroExponent, float oneExponent) {
	if (zeroExponent <= 0) throw IllegalArgumentError("zeroExponent must be positive.", __FILE__, __LINE__);
	if (oneExponent <= 0) throw IllegalArgumentError("oneExponent must be positive.", __FILE__, __LINE__);
	return RealFunction([begin, end, zeroExponent, oneExponent](float t) {
		float x = (t-begin)/(end-begin);
		if (x < 0) return 0.f;
		if (x > 1) return 0.f;
		float c = std::pow(zeroExponent + oneExponent, zeroExponent + oneExponent)/(std::pow(zeroExponent, zeroExponent)*std::pow(oneExponent, oneExponent));
		return c*std::pow(x, zeroExponent)*std::pow(1.f - x, oneExponent);
	});
}

RealFunction toneMap(float k) {
	return RealFunction([k](float t) {
		return (k*t)/(1+k*t);
	});
}

RealFunction rationalInfiniteShroom(float steepFactor, float center) {
	if (steepFactor <= 0) throw IllegalArgumentError("steepFactor must be positive.", __FILE__, __LINE__);
	return RealFunction([steepFactor, center](float t) {
			return 1.f/(steepFactor*sq(t-center) + 1);
		});
}
