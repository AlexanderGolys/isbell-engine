#include "smoothImplicit.hpp"
#include "src/fundamentals/func.hpp"
#include "src/common/indexedRendering.hpp"

using std::vector, std::string, std::shared_ptr, std::unique_ptr, std::pair, std::make_unique, std::make_shared, std::function;

using namespace glm;

SmoothImplicitSurface::SmoothImplicitSurface(const RealFunctionR3 &F) {
    _F = F;
}

float SmoothImplicitSurface::operator()(vec3 p) const {
    return _F(p);
}

AffinePlane::AffinePlane(vec3 n, float d) :
SmoothImplicitSurface(RealFunctionR3([n, d](vec3 p) {return dot(p, n) - d; }, [n](vec3 p) {return n; })) {
    this->n = normalise(n);
    this->d = d;
    auto tangents = orthogonalComplementBasis(n);
    this->v1 = tangents.first;
    this->v2 = tangents.second;
    this->pivot = this->n * d;
}



AffinePlane::AffinePlane(vec3 pivot, vec3 v1, vec3 v2):
SmoothImplicitSurface(RealFunctionR3(
    [v1, v2, pivot](vec3 p) {return dot(p, normalise(cross(v1, v2))) - dot(pivot, normalise(cross(v1, v2))); },
    [v1, v2](vec3 p) {return normalise(cross(v1, v2)); })) {
    this->pivot = pivot;
    this->v1 = v1;
    this->v2 = v2;
    this->n = normalise(cross(v1, v2));
    this->d = dot(pivot, n);
}
    AffinePlane AffinePlane::spanOfPts(vec3 p0, vec3 p1, vec3 p2) {
    return AffinePlane(p0, p1 - p0, p2 - p0);
}

AffineLine AffinePlane::intersect(const AffinePlane &p) const {
	vec3 v = cross(n, p.normal());
	float denominator = norm2(n)*norm2(p.normal()) - dot(n, p.normal())*dot(n, p.normal());
	vec3 pivot = (d*norm2(p.normal()) - p.getD()*dot(n, p.normal()))/denominator +
				p.normal()*(p.getD()*norm2(n) - d*dot(n, p.normal()))/denominator;
	return AffineLine(pivot, v);
}


vec3 AffinePlane::intersect(const AffineLine &l) const {
	return l.pivot() + dot(pivot - l.pivot(), n) / dot(l.direction(), n) * l.direction();
}

mat3 AffinePlane::pivotAndBasis() const {
    return mat3(pivot, v1, v2);
}

vec2 AffinePlane::localCoordinates(vec3 p) const {
    return vec2(dot(p - pivot, v1), dot(p - pivot, v2));
}



SmoothParametricCurve::SmoothParametricCurve(Foo13 f,Foo13 df,  Foo13 ddf, PolyGroupID id, float t0, float t1, bool periodic, float epsilon)
{
    this->_f = f;
    this->_df = df;
    this->_ddf = ddf;
    this->eps = epsilon;
    this->t0 = t0;
    this->t1 = t1;
    this->periodic = periodic;
    this->id = id;
}

SmoothParametricCurve::SmoothParametricCurve(Foo13 f, Foo13 df, PolyGroupID id, float t0, float t1, bool periodic, float epsilon)
: SmoothParametricCurve(f, df, derivativeOperator(df, epsilon), id, t0, t1, periodic, epsilon) {}



SmoothParametricCurve::SmoothParametricCurve(const SmoothParametricCurve &other) :
    _f(other._f), _df(other._df), _ddf(other._ddf), _der_higher(other._der_higher), eps(other.eps), t0(other.t0), t1(other.t1),
    periodic(other.periodic), id(other.id){}

SmoothParametricCurve::SmoothParametricCurve(SmoothParametricCurve &&other) noexcept :
    _f(std::move(other._f)), _df(std::move(other._df)), _ddf(std::move(other._ddf)), _der_higher(std::move(other._der_higher)),
    eps(other.eps), t0(other.t0), t1(std::move(other.t1)), periodic(other.periodic), id(other.id) {}

SmoothParametricCurve &SmoothParametricCurve::operator=(const SmoothParametricCurve &other) {
    if (this == &other)
        return *this;
    _f = other._f;
    _df = other._df;
    _ddf = other._ddf;
    _der_higher = other._der_higher;
    eps = other.eps;
    t0 = other.t0;
    t1 = other.t1;
    periodic = other.periodic;
    id = other.id;
    return *this;
}
SmoothParametricCurve &SmoothParametricCurve::operator=(SmoothParametricCurve &&other) noexcept {
    if (this == &other)
        return *this;
    _f = std::move(other._f);
    _df = std::move(other._df);
    _ddf = std::move(other._ddf);
    _der_higher = std::move(other._der_higher);
    eps = other.eps;
    t0 = std::move(other.t0);
    t1 = std::move(other.t1);
    periodic = other.periodic;
    id = other.id;
    return *this;
}

SmoothParametricCurve::SmoothParametricCurve(Foo13 f, std::vector<Foo13> derivatives, PolyGroupID id, float t0, float t1, bool periodic, float epsilon) {
	this->_f = f;
	if (derivatives.size() > 0)
		this->_df = derivatives[0];
	else
		this->_df = derivativeOperator(f, epsilon);

	if (derivatives.size() > 1)
		this->_ddf = derivatives[1];
	else
		this->_ddf = derivativeOperator(_df, epsilon);

	this->_der_higher = [derivatives, f](int n) { return n==0 ? f : n <= derivatives.size() ? derivatives[n-1] : derivativeOperator(derivatives[n-1], 0.01f);} ;
	this->eps = epsilon;
    this->t0 = t0;
    this->t1 = t1;
    this->periodic = periodic;
    id = randomCurvaID();
}



SmoothParametricCurve SmoothParametricCurve::constCurve(vec3 v) {
	return SmoothParametricCurve([v](float t) {return v; },-10, 10, false, 0.01);
}

float SmoothParametricCurve::length(float t0, float t1, int n) const {
    float res = 0;
    float dt = (t1 - t0) / n;
    for (int i = 0; i < n; i++) {
        res += norm(_f(lerp(t0, t1, dt * i)) - _f(lerp(t0, t1, dt * (i + 1))));
    }
    return res;
}
SmoothParametricCurve SmoothParametricCurve::precompose(SpaceEndomorphism g_) const {
    return SmoothParametricCurve([f = this->_f, g=g_](float t) {return g(f(t)); },
                                  [f = this->_f, d = this->_df, g=g_](float t) {return g.df(f(t)) * d(t); },
                                  id, this->t0, this->t1, this->periodic, this->eps);
}
void SmoothParametricCurve::precomposeInPlace(SpaceEndomorphism g) {
    this->_f = [f = this->_f, gg=g](float t) {return gg(f(t)); };
    this->_df = [f = this->_f, d = this->_df, g](float t) {return g(f(t)) * d(t); };
    this->_ddf = [f = this->_f, d = this->_df, dd = this->_ddf, gg=g](float t) {return gg(f(t)) * dd(t) + gg.df(f(t)) * d(t); };
}

vec3 SmoothParametricCurve::operator()(float t) const { return this->_f(t); }

mat3 SmoothParametricCurve::FrenetFrame(float t) const {
	return mat3(tangent(t), normal(t), binormal(t));

}

float SmoothParametricCurve::curvature(float t) const {
	return norm(cross(_ddf(t), _df(t))) / pow(norm(_df(t)), 3);
}

float SmoothParametricCurve::torsion(float t) const {
	return dot(cross(_df(t), ddf(t)), higher_derivative(t, 3)) / norm2(cross(df(t), ddf(t)));
}

vec3 SmoothParametricCurve::curvature_vector(float t) const {
	return normal(t)/curvature(t);
}

AffineLine::AffineLine(vec3 p0, vec3 v) :
SmoothParametricCurve([p0, v](float t) {return p0 + t * normalise(v); },
					  [v](float t) {return normalise(v); },
					  [](float t) {return vec3(0, 0, 0); }) {
	this->p0 = p0;
	this->v = normalise(v);
}

AffineLine AffineLine::spanOfPts(vec3 p0, vec3 p1) {
	return AffineLine(p0, p1 - p0);
}


RealFunctionR3 AffineLine::distanceField() const {
	return RealFunctionR3([this](vec3 p) {return distance(p); });
}

vec3 AffineLine::orthogonalProjection(vec3 p) const {
	return p0 + dot(p - p0, v) * v / dot(v, v);
}

vec3 AffineLine::pivot() const {
	return p0;
}

vec3 AffineLine::direction() const {
	return v;
}

float AffineLine::distance(AffineLine &l) const {
	return dot((p0 - l.pivot()), cross(v, l.direction())) / norm(cross(v, l.direction()));
}

AffineLine AffineLine::operator+(vec3 v) const {
	return AffineLine(p0 + v, this->v);
}


float AffineLine::distance(vec3 p) const {
	return norm(cross(p - p0, v)) / norm(v);
}

bool AffineLine::contains(vec3 p, float eps) const {
	return distance(p) < eps;
}



ImplicitSurfacePoint::ImplicitSurfacePoint(vec3 p, vec3 normal, bool border) : p(p), border(border) {
	auto tangents = orthogonalComplementBasis(normal);
	orthoFrame = mat3(tangents.first, tangents.second, normal);
}

vec3 ImplicitSurfacePoint::projectOnTangentPlane(vec3 q) const {
	vec3 ortho_coords = coords_in_frame(q);
	return getTangent1() * ortho_coords.x + getTangent2() * ortho_coords.y;
}
vec3 ImplicitSurfacePoint::rotateAroundNormal(vec3 q, float angle) const {
	vec3 projected = projectOnTangentPlane(q);
	return projected*cos(angle) + cross(projected, getNormal())*sin(angle);
}

void FrontPolygon::recalculate_angle(int index) {
	if (points.size() < 3)
		return;
	vec3 p0 = points[(index - 1 + points.size()) % points.size()].getPosition();
	vec3 p1 = points[index].getPosition();
	vec3 p2 = points[(index + 1) % points.size()].getPosition();
	vec3 v1 = p0 - p1;
	vec3 v2 = p2 - p1;
	float angle = acos(dot(v1, v2) / (norm(v1)*norm(v2)));
	if (dot(cross(v1, v2), points[index].getNormal()) < 0)
		angle = TAU - angle;
	points[index].setAngle(angle);
}
void FrontPolygon::recalculate_angles() { for (int i = 0; i < points.size(); i++) recalculate_angle(i); }

int FrontPolygon::argminAngle() const {
	int argmin = 0;
	float min = points[0].getAngle();
	for (int i = 1; i < points.size(); i++) {
		if (points[i].getAngle() < min) {
			min = points[i].getAngle();
			argmin = i;
		}
	}
	return argmin;
}
float FrontPolygon::minAngle() const { return points[argminAngle()].getAngle(); }

void FrontPolygon::addPoint(ImplicitSurfacePoint p, int index) {
	points.insert(points.begin() + index, p);
	recalculate_angle((index - 1 + points.size()) % points.size());
	recalculate_angle(index);
	recalculate_angle((1 + index) % points.size());
	for (int i = 0; i < excluded_for_distance_check.size(); i++) {
		if (excluded_for_distance_check[i] >= index)
			excluded_for_distance_check[i]++;
	}
}

bool FrontPolygon::removePoint(int index) {
	if (points.size() <= 2)
		return true;
	points.erase(points.begin() + index);
	recalculate_angle((index - 1 + points.size()) % points.size());
	recalculate_angle(index % points.size());
	for (int i = 0; i < excluded_for_distance_check.size(); i++) {
		if (excluded_for_distance_check[i] == index) {
			excluded_for_distance_check.erase(excluded_for_distance_check.begin() + i);
			i--;
		}
		else if (excluded_for_distance_check[i] > index)
			excluded_for_distance_check[i]--;
	}
	return false;
}

void FrontPolygon::expandVertex(int index, TriangulatedImplicitSurface &surface) {
	if (points.size() < 3)
		return;
	vec3 p0 = points[(index - 1 + points.size()) % points.size()].getPosition();
	vec3 p1 = points[index].getPosition();
	vec3 p2 = points[(index + 1) % points.size()].getPosition();
	vec3 v1 = p0 - p1;
	vec3 v2 = p2 - p1;
	float angle = points[index].getAngle();
	int divisions = max(1, 3*(int)(angle/PI) + 1);
	float delta = angle / divisions;

	if (delta < .8 && divisions > 1) {
		divisions--;
		delta = angle / divisions;
	}
	if (divisions == 1 && delta > .8 && norm(v1-v2) > 1.25*side) {
		divisions = 2;
		delta = angle / divisions;
	}
	if (angle<3 && min(norm(v1), norm(v2)) < .5*side) {
		divisions = 1;
		delta = angle;
	}
	std::vector verts = { points[(index - 1 + points.size()) % points.size()]};
	for (int i = 1; i < divisions; i++)
		verts.push_back(surface.findNearbyPoint( p1 + normalize(points[index].rotateAroundNormal(v1, delta*i)) * side));
	verts.push_back(points[(index + 1) % points.size()]);

	for (int i = 0; i < verts.size()-1; i++)
		surface.addTriangle({ verts[i].getPosition(), verts[i+1].getPosition(), p1 });

	removePoint(index);
	for (int i = 0; i < verts.size()-1; i++) {
		addPoint(verts[i], index + i);
		index++;
	}
}
void FrontPolygon::step(TriangulatedImplicitSurface &surface) { expandVertex(argminAngle(), surface); }
bool FrontPolygon::checkForSelfIntersections() { return true; }
bool FrontPolygon::merge(FrontPolygon &other, int ind_self, int ind_other) { return true; }
bool FrontPolygon::checkForCrossIntersections(const FrontPolygon &other) { return true; }


void TriangulatedImplicitSurface::removeDegeneratePolygons() {
	for (int i = 0; i < polygons.size(); i++) {
		if (polygons[i].size() < 3) {
			polygons.erase(polygons.begin() + i);
			i--;
		}
	}
}
ImplicitSurfacePoint TriangulatedImplicitSurface::findNearbyPoint(vec3 p) const {
	for (int i = 0; i < max_iter; i++) {
		vec3 next = p - F(p)*F.df(p)/norm2(F.df(p));
		if (norm(next - p) < eps) {
			p = next;
			break;
		}
		p = next;
	}
	vec3 normal = normalise(F.df(p));
	return ImplicitSurfacePoint(p, normal, false);
}
void TriangulatedImplicitSurface::initHexagon(vec3 p, float side) {
	vec3 p0 = findNearbyPoint(p).getPosition();
	vector<ImplicitSurfacePoint> points = {};
	for (int i = 0; i < 6; i++) {
		vec3 next = p0 + side * cos(i * TAU / 6) * vec3(1, 0, 0) + side * sin(i * TAU / 6) * vec3(0, 1, 0);
		points.push_back(findNearbyPoint(next));
	}
	triangles = { {points[0].getPosition(), points[1].getPosition(), points[2].getPosition()},
				  {points[0].getPosition(), points[2].getPosition(), points[3].getPosition()},
				  {points[0].getPosition(), points[3].getPosition(), points[4].getPosition()},
				  {points[0].getPosition(), points[4].getPosition(), points[5].getPosition()},
				  {points[0].getPosition(), points[5].getPosition(), points[1].getPosition()},
				  {points[1].getPosition(), points[2].getPosition(), points[3].getPosition()},
				  {points[3].getPosition(), points[4].getPosition(), points[5].getPosition()},
				  {points[5].getPosition(), points[1].getPosition(), points[2].getPosition()} };
	polygons = {FrontPolygon(points, side)};
	front_polygon = 0;

}

bool TriangulatedImplicitSurface::step() {
	if (polygons.empty()) return false;
	for (int i = 0; i < polygons.size(); i++)
		polygons[i].noExcluded();
	polygons[front_polygon].step(*this);
	// todo

	current_step++;
	return true;
}
WeakSuperMesh TriangulatedImplicitSurface::compute(int max_iter, vec3 p0, float side) {
	initHexagon(p0, side);
	while (step() && current_step < max_iter);
	vector<Vertex> vertices = {};
	vector<ivec3> trs = {};
	for (auto &tr : triangles) {
		vertices.emplace_back(tr[0], vec2(tr[0].x, tr[0].y));
		vertices.emplace_back(tr[1], vec2(tr[1].x, tr[1].y));
		vertices.emplace_back(tr[2], vec2(tr[2].x, tr[2].y));
		trs.emplace_back(vertices.size()-3, vertices.size()-2, vertices.size()-1);
	}
	return WeakSuperMesh(vertices, trs, 0);
}

AffinePlane SmoothParametricCurve::osculatingPlane(float t) const {
	return AffinePlane(binormal(t), dot(binormal(t), _f(t)));
}
