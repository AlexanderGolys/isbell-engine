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
    for (int i = 0; i < n; i++)
		res += norm(_f(t0 + dt*i) - _f(t0 + dt*(i+1)));
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




vec3 ImplicitSurfacePoint::projectOnTangentPlane(vec3 q) const {
	vec3 ortho_coords = coords_in_frame(q);
	return getTangent1() * ortho_coords.x + getTangent2() * ortho_coords.y;
}
vec3 ImplicitSurfacePoint::rotateAroundNormal(vec3 q, float angle) const {
	vec3 projected = projectOnTangentPlane(q);
	return projected*cos(angle) + cross(projected, getNormal())*sin(angle);
}

float TriangulatedImplicitSurface::angle(int i) { return points[i].getAngle(); }


ImplicitSurfacePoint TriangulatedImplicitSurface::projectOnSurface(vec3 p) {
	vec3 last_q = p;
	for (int i = 0; i < NewtonMaxSteps; i++) {
		vec3 q = last_q - F.df(last_q)*F(last_q)/norm2(F.df(last_q));
		if (norm(q - last_q) < NewtonEps) {
			last_q = q;
			break;
		}
		last_q = q;
	}
	return constructPoint(last_q);
}

ImplicitSurfacePoint TriangulatedImplicitSurface::constructPoint(vec3 p) {
	vec3 normal = normalise(F.df(p));
	vec3 t1;
	if (normal.x > .5 || normal.y>.5)
		t1 = normalise(vec3(normal.y, -normal.x, 0));
	else
		t1 = normalise(vec3(-normal.z, 0, normal.x));
	vec3 t2 = cross(normal, t1);
	return ImplicitSurfacePoint(p, mat3(normal, t1, t2), age, borderCheck(p));
}


void TriangulatedImplicitSurface::initialiseHexagon(vec3 p0, float len) {
	points = {projectOnSurface(p0)};
	for (int i = 0; i < 6; i++) {
		vec3 p = p0 + len*cos(i*TAU / 6)*vec3(1, 0, 0) + len*sin(i*TAU / 6)*vec3(0, 1, 0);
		points.push_back(projectOnSurface(p));
	}
	polygons = {{1, 2, 3, 4, 5, 6}};
}


void TriangulatedImplicitSurface::expandFrontPolygon() { throw std::logic_error("Not implemented"); }


void TriangulatedImplicitSurface::calculateMissingAngles() {
}


void TriangulatedImplicitSurface::calculateAngle(int i) {
	int ind_p = polygons[0][i];
	vec3 p = points[ind_p].getPosition();
	vec3 r = points[polygons[0][(i + 1) % 6]].getPosition();
	vec3 q = points[polygons[0][(i - 1 + polygons[0].size()) % 6]].getPosition();
	vec3 v1 = q-p;
	vec3 v2 = r-p;
	float angle = acos(dot(v1, v2) / (norm(v1)*norm(v2)));
	if (dot(cross(v1, v2), points[ind_p].getNormal()) < 0)
		angle = TAU - angle;
	points[ind_p].setAngle(angle);
}

int TriangulatedImplicitSurface::minAngleIndex() {
	int minInd = 0;
	float minAngle = TAU;
	for (int i = 0; i < polygons[0].size(); i++) {
			int ind = polygons[0][i];
			if (points[ind].getAngle() < minAngle && !points[ind].border) {
				minAngle = points[ind].getAngle();
				minInd = ind;
			}
		}
	return minInd;
}



void TriangulatedImplicitSurface::checkFrontPolygonSelfIntersections() { throw std::logic_error("Not implemented"); }
void TriangulatedImplicitSurface::checkFrontPolygonIntersectionWithPolygon(int i) { throw std::logic_error("Not implemented"); }


void TriangulatedImplicitSurface::expandFrontPolygon(int i) { throw std::logic_error("Not implemented"); }


AffinePlane SmoothParametricCurve::osculatingPlane(float t) const {
	return AffinePlane(binormal(t), dot(binormal(t), _f(t)));
}
