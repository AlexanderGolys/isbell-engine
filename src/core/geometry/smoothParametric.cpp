#include "smoothParametric.hpp"
#include <map>
// #include <glm/gtc/constants.hpp>
#include <utility>

using namespace glm;
using std::vector, std::string, std::shared_ptr, std::unique_ptr, std::pair, std::make_unique, std::make_shared, std::function;





SmoothParametricCurve SmoothParametricSurface::precompose(const SmoothParametricPlaneCurve& c, PolyGroupID id) const {
    return SmoothParametricCurve([f=*this, g=c](float t) {return f(g(t)); }, std::move(id), c.bounds().x, c.bounds().y, c.isPeriodic(), epsilon);
}

SmoothParametricSurface SmoothParametricSurface::postcompose(const SpaceEndomorphism& g) const {
	return SmoothParametricSurface([f=_f, g](float t, float s) {return g(f(t, s)); },
//		[df=_df_t, _f=_f, g](float t, float s) {return g.df(_f(t, s))*df(t, s); },
//		[df=_df_u, _f=_f, g](float t, float s) {return g.df(_f(t, s))*df(t, s); },
		rangeT, rangeU, periodicT, periodicU, epsilon);
}

SmoothParametricSurface SmoothParametricSurface::precompose(const PlaneSmoothEndomorphism &g, vec2 t_bounds, vec2 u_bounds) const {
	return SmoothParametricSurface([f = _f, g](float t, float s) {return f(g(t, s).x, g(t, s).y); },
		t_bounds, u_bounds, periodicT, periodicU, epsilon);
}

SmoothParametricSurface SmoothParametricSurface::precompose(const PlaneAutomorphism &g) const {
	vec2 new0 = g.inv(t0(), u0());
	vec2 new1 = g.inv(t1(), u1());
	return precompose(g,vec2(new0.x, new1.x), vec2(new0.y, new1.y));
}

SmoothParametricCurve SmoothParametricSurface::restrictToInterval(vec2 p0, vec2 p1, PolyGroupID id) const {
    return SmoothParametricCurve([f=*this, p0, p1](float t) {return f(p1*t + p0*(1-t)); }, id, 0, 1,
    	nearlyEqual((*this)(p0), (*this)(p1)), epsilon);
}

SmoothParametricCurve SmoothParametricSurface::restrictToInterval(vec2 p0, vec2 p1, bool periodic,  PolyGroupID id) const {
	return SmoothParametricCurve([f=*this, p0, p1](float t) {return f(p1*t + p0*(1-t)); }, id, 0, 1,
		periodic, epsilon);
}

SmoothParametricCurve SmoothParametricSurface::constT(float ti) const {
	return SmoothParametricCurve([f=_f, ti](float u) {return f(ti, u); },
								[df=_df_u, ti](float u) {return df(ti, u); },
								randomID(), u0(), u1(), periodicU, epsilon);

}
SmoothParametricCurve SmoothParametricSurface::constU(float ui) const {
	return SmoothParametricCurve([f=_f, ui](float t) {return f(t, ui); },
								[df=_df_t, ui](float t) {return df(t, ui); },
								randomID(), t0(), t1(), periodicT, epsilon);
}




SmoothParametricSurface AffineLine::tube(float radius, float t0, float t1) const {
	auto w = orthogonalComplementBasis(direction());
	SmoothParametricCurve circle = SmoothParametricCurve([radius, w, p=(*this)(t0)](float t) {return  p + radius*cos(t)*w.first + radius*sin(t)*w.second; });
	return circle.cylinder(direction(), t1 - t0);
}

SmoothParametricSurface::SmoothParametricSurface(const Foo113& f, const Foo113& df_t, const Foo113& df_u,
	vec2 t_range, vec2 u_range,  bool t_periodic, bool u_periodic, float epsilon)  :
	_f(f), _df_t(df_t), _df_u(df_u), rangeT(t_range), rangeU(u_range), periodicT(t_periodic), periodicU(u_periodic), epsilon(epsilon) {}

SmoothParametricSurface::SmoothParametricSurface(const std::function<SmoothParametricCurve(float)>& pencil, vec2 t_range, vec2 u_range, bool t_periodic, bool u_periodic, float eps)
	:	_f([pencil](float t, float u) {return pencil(t)(u); }), _df_t([pencil](float t, float u) {return pencil(t).tangent(u); }), _df_u([pencil](float t, float u) {return pencil(t).normal(u); }),
		rangeT(t_range), rangeU(u_range), periodicT(t_periodic), periodicU(u_periodic), epsilon(eps) {}

SmoothParametricSurface::SmoothParametricSurface(RealFunctionR2 plot, vec2 t_range, vec2 u_range): SmoothParametricSurface([plot](float t, float u) { return vec3(t, u, plot(t, u)); },
																														   [plot](float t, float u) { return vec3(1, 0, plot.dx(t, u)); },[plot](float t, float u) { return vec3(0, 1, plot.dy(t, u)); },
																														   t_range, u_range) {}



mat3 SmoothParametricSurface::DarbouxFrame(float t, float s) const { throw std::logic_error("Not implemented"); }

//mat2x3 SmoothParametricSurface::tangentSpace(float t, float s) const {
//	vec3 p = _df_t(t, s);
//	vec3 q = _df_u(t, s);
//	return mat2x3(p, q);
//}

mat2 SmoothParametricSurface::changeOfTangentBasisToPrincipal(float t, float s) const { throw std::logic_error("Not implemented"); }
mat2 SmoothParametricSurface::changeOfPrincipalBasisToStandard(float t, float s) const { throw std::logic_error("Not implemented"); }
mat2x3 SmoothParametricSurface::tangentSpacePrincipalBasis(float t, float s) const { throw std::logic_error("Not implemented"); }
float SmoothParametricSurface::globalAreaIntegral(const RealFunctionPS &f) const { throw std::logic_error("Not implemented"); }
float SmoothParametricSurface::DirichletFunctional() const { throw std::logic_error("Not implemented"); }
float SmoothParametricSurface::biharmonicFunctional() const { throw std::logic_error("Not implemented"); }

void SmoothParametricSurface::changeDomain(vec2 t_range, vec2 u_range, bool t_periodic, bool u_periodic) {
	rangeT = t_range;
	rangeU = u_range;
	this->periodicT = t_periodic;
	this->periodicU = u_periodic;
}

mat2 SmoothParametricSurface::metricTensor(float t, float s) const {
	vec3 p = _df_t(t, s);
	vec3 q = _df_u(t, s);
	return mat2(dot(p, p), dot(p, q), dot(q, p), dot(q, q));
}


mat2 SmoothParametricSurface::firstFundamentalForm(float t, float s) const {
	return mat2(_E(t, s), _F(t, s), _F(t, s), _G(t, s));
}

mat2 SmoothParametricSurface::secondFundamentalForm(float t, float s) const {
	return mat2(_L(t, s), _M(t, s), _M(t, s), _N(t, s));
}

float SmoothParametricSurface::_E(float t, float s) const {
	return norm2(_df_t(t, s));
}

float SmoothParametricSurface::_F(float t, float s) const {
	return dot(_df_t(t, s), _df_u(t, s));
}

float SmoothParametricSurface::_G(float t, float s) const {
	return norm2(_df_u(t, s));
}
vec3 SmoothParametricSurface::d2f_tt(float t, float s) const {
	return (_df_t(t + epsilon, s) - _df_t(t - epsilon, s)) / (2 * epsilon);
}
vec3 SmoothParametricSurface::d2f_uu(float t, float u) const {
	return (_df_u(t, u + epsilon) - _df_u(t, u - epsilon)) / (2 * epsilon);
}

vec3 SmoothParametricSurface::d2f_tu (float t, float u) const {
	return (_df_t(t, u + epsilon) - _df_t(t, u - epsilon)) / (2 * epsilon);
}


float SmoothParametricSurface::_L(float t, float s) const {
	return dot(d2f_tt(t, s), normal(t, s));
}

float SmoothParametricSurface::_M(float t, float s) const {
	return dot(d2f_tu(t, s), normal(t, s));
}

float SmoothParametricSurface::_N(float t, float s) const {
	return dot(d2f_uu(t, s), normal(t, s));
}

vec3 SmoothParametricSurface::normalCurvature(float t, float s, vec2 v) const { throw std::logic_error("Not implemented"); }

mat2 SmoothParametricSurface::shapeOperator(float t, float s) const { throw std::logic_error("Not implemented"); }

float SmoothParametricSurface::meanCurvature(float t, float s) const {
	return (_E(t, s) * _N(t, s) + _G(t, s) * _L(t, s) - 2 * _F(t, s) * _M(t, s)) / (2 * (_E(t, s) * _G(t, s) - _F(t, s) * _F(t, s)));
}

float SmoothParametricSurface::gaussianCurvature(float t, float s) const {
	return (_L(t, s) * _N(t, s) - _M(t, s) * _M(t, s)) / (_E(t, s) * _G(t, s) - _F(t, s) * _F(t, s)); }

mat2x3 SmoothParametricSurface::principalDirections(float t, float s) const { throw std::logic_error("Not implemented"); }

std::pair<float, float> SmoothParametricSurface::principalCurvatures(float t, float s) const {
	return {meanCurvature(t, s) - sqrt(meanCurvature(t, s)*meanCurvature(t, s) - gaussianCurvature(t, s)),
			meanCurvature(t, s) + sqrt(meanCurvature(t, s)*meanCurvature(t, s) - gaussianCurvature(t, s))};
}

vec3 SmoothParametricSurface::Laplacian(float t, float s) const { throw std::logic_error("Not implemented"); }

SmoothParametricSurface SmoothParametricSurface::meanCurvatureFlow(float dt) const {
	return SmoothParametricSurface([S=*this, dt](float t, float s) {return S(t, s) + S.normal(t, s)*S.meanCurvature(t, s)*dt; },
		rangeT, rangeU, periodicT, periodicU);
}

SmoothParametricSurface::SmoothParametricSurface(const Foo113 &f, vec2 t_range, vec2 u_range, bool t_periodic, bool u_periodic, float epsilon)
	: _f(f), rangeT(t_range), rangeU(u_range), periodicT(t_periodic), periodicU(u_periodic), epsilon(epsilon)
{
	_df_t = [f, epsilon](float t, float u) {return (f(t + epsilon, u) - f(t - epsilon, u)) / (2 * epsilon); };
	_df_u = [f, epsilon](float t, float u) {return (f(t, u + epsilon) - f(t, u - epsilon)) / (2 * epsilon); };
}

vec3 SmoothParametricSurface::operator()(float t, float s) const {
	return _f(t, s);
}

vec3 SmoothParametricSurface::operator()(vec2 tu) const { return _f(tu.x, tu.y); }

SmoothParametricSurface SmoothParametricSurface::operator+(const SmoothParametricSurface &S) const {
	return SmoothParametricSurface([f=_f, g=S._f](float t, float u) {return f(t, u) + g(t, u); },
		[df_t=_df_t, dg_t=S._df_t](float t, float u) {return df_t(t, u) + dg_t(t, u); },
		[df_u=_df_u, dg_u=S._df_u](float t, float u) {return df_u(t, u) + dg_u(t, u); },
		rangeT, rangeU, periodicT&&S.periodicT, periodicU&&S.periodicU, min(epsilon, S.epsilon));
}

SmoothParametricSurface SmoothParametricSurface::operator*(float a) const {
	return SmoothParametricSurface([f=_f, a](float t, float u) {return a*f(t, u); },
		[df_t=_df_t, a](float t, float u) {return a*df_t(t, u); },
		[df_u=_df_u, a](float t, float u) {return a*df_u(t, u); },
		rangeT, rangeU, periodicT, periodicU, epsilon);
}

SmoothParametricSurface SmoothParametricSurface::normaliseParameters() const {
	return SmoothParametricSurface(
			[f=_f, t0=t0(), t1=t1(), u0=u0(), u1=u1()](float t, float u) {return f(lerp(t0, t1, t), lerp(u0, u1, u)); },
			[df=_df_t, t0=t0(), t1=t1(), u0=u0(), u1=u1()](float t, float u) {return (t1-t0)*df(lerp(t0, t1, t), lerp(u0, u1, u)); },
			[df=_df_u, t0=t0(), t1=t1(), u0=u0(), u1=u1()](float t, float u) {return (u1-u0)*df(lerp(t0, t1, t), lerp(u0, u1, u)); },
			vec2(0, 1), vec2(0, 1), periodicT, periodicU, epsilon);
}


SmoothParametricSurface SmoothParametricSurface::shift(vec3 v) const {
	auto shift = SpaceAutomorphism::translation(v);
	return postcompose(shift);
}

SmoothParametricSurface SmoothParametricSurface::scale(float a, vec3 center) const {
	auto scale = SpaceAutomorphism::scaling(a).applyWithShift(center);
	return postcompose(scale);
}

SmoothParametricSurface SmoothParametricSurface::rotate(vec3 axis, float angle, vec3 center) const {
	auto rot = SpaceAutomorphism::rotation(axis, angle).applyWithShift(center);
	return postcompose(rot);
}

vec3 SmoothParametricSurface::normal(float t, float u) const {
    vec3 vn = cross(_df_t(t, u), _df_u(t, u));
	// return norm(vn) < .01 ? vec3(0, 0, 1) : normalise(vn);

    if (norm(vn) < .0002) {
        float e_t = (t1() - t0()) / 400;
        float e_u = (u1() - u0()) / 400;
        vec3 new_n = -(cross(_df_t(t + e_t, u), _df_u(t + e_t, u)) + cross(_df_t(t, u + e_u), _df_u(t, u + e_u)) +
                      cross(_df_t(t - e_t, u), _df_u(t - e_t, u)) + cross(_df_t(t, u - e_u), _df_u(t, u - e_u))) /
                     4.f;

        if (norm(new_n) < .0003) {
            e_t = (t1() - t0()) / 200;
            e_u = (u1() - u0()) / 200;
            new_n = -(cross(_df_t(t + e_t, u), _df_u(t + e_t, u)) + cross(_df_t(t, u + e_u), _df_u(t, u + e_u)) +
                     cross(_df_t(t - e_t, u), _df_u(t - e_t, u)) + cross(_df_t(t, u - e_u), _df_u(t, u - e_u))) /
                    4.f;
            return norm(new_n) < .01 ? vec3(0, 0, 1) : normalise(new_n);
        }

        return -normalise(new_n);
    }
    return -normalize(vn);
}
void SmoothParametricSurface::normaliseDomainToI2() {
    _f = [f=_f, t0=t0(), t1=t1(), u0=u0(), u1=u1()](float t, float u) {return f(lerp(t0, t1, t), lerp(u0, u1, u)); };
    _df_t = [df=_df_t, t0=t0(), t1=t1(), u0=u0(), u1=u1()](float t, float u) {return (t1-t0)*df(lerp(t0, t1, t), lerp(u0, u1, u)); };
    _df_u = [df=_df_u, t0=t0(), t1=t1(), u0=u0(), u1=u1()](float t, float u) {return (u1-u0)*df(lerp(t0, t1, t), lerp(u0, u1, u)); };
    epsilon /= max(t1()-t0(), u1()-u0());
	rangeT = vec2(0, 1);
	rangeU = vec2(0, 1);
}

SmoothParametricCurve SmoothParametricPlaneCurve::embedding(vec3 v1, vec3 v2, vec3 pivot) const {
    auto aff = mat3(v1, v2, cross(v1, v2));
    return SmoothParametricCurve([pivot, aff, f=this->_f](float t) {return aff*vec3(f(t), 0) + pivot; },
                                           [aff, d=this->_df](float t) {return aff*vec3(d(t), 0); },
                                           [aff, dd=this->_ddf](float t) {return aff*vec3(dd(t), 0); });
}

SmoothParametricCurve::SmoothParametricCurve(const std::function<vec3(float)>& f, PolyGroupID id, float t0, float t1, bool periodic, float epsilon)
: SmoothParametricCurve(f, derivativeOperator(f, epsilon), id, t0, t1, periodic, epsilon) {}

SmoothParametricCurve::SmoothParametricCurve(const std::function<vec3(float)>& f, vec2 dom, PolyGroupID id, bool periodic, float epsilon): SmoothParametricCurve(f, derivativeOperator(f, epsilon), id, dom.x, dom.y, periodic, epsilon) {}

PolyGroupID SmoothParametricCurve::getID() const {
	return id;
}

void SmoothParametricCurve::setID(PolyGroupID id) {
	this->id = id;
}

void SmoothParametricCurve::copyID(const SmoothParametricCurve& other) {
	this->id = other.id;
}

SmoothParametricCurve SmoothParametricCurve::shift(vec3 v) const {
	return precompose(SpaceEndomorphism::translation(v));
}

SmoothParametricCurve SmoothParametricCurve::rotate(vec3 axis, float angle, vec3 center) const {
	return precompose(SpaceAutomorphism::rotation(axis, angle).applyWithShift(center));
}

SmoothParametricCurve SmoothParametricCurve::scale(float a, vec3 center) const {
	return precompose(SpaceAutomorphism::scaling(a, center));
}

vec3 SmoothParametricCurve::derivative(float t) const {
	return _df(t);
}

vec3 SmoothParametricCurve::df(float t) const {
	return derivative(t);
}

vec2 SmoothParametricCurve::bounds() const {
	return vec2(t0, t1);
}

float SmoothParametricCurve::getT0() const {
	return t0;
}

float SmoothParametricCurve::getT1() const {
	return t1;
}

vec3 SmoothParametricCurve::second_derivative(float t) const {
	return _ddf(t);
}

vec3 SmoothParametricCurve::higher_derivative(float t, int n) const {
	return _der_higher(n)(t);
}

vec3 SmoothParametricCurve::ddf(float t) const {
	return second_derivative(t);
}

float SmoothParametricCurve::curvature_radius(float t) const {
	return 1 / curvature(t);
}

float SmoothParametricCurve::speed(float t) const {
	return norm(df(t));
}

bool SmoothParametricCurve::isPeriodic() const {
	return periodic;
}

float SmoothParametricCurve::getEps() const {
	return eps;
}

vec3 SmoothParametricSurface::parametersNormalised(vec2 tu) const {
	return operator()(t0() + tu.x * (t1() - t0()), u0() + tu.y * (u1() - u0()));
}

vec3 SmoothParametricSurface::parametersNormalised(float t, float u) const {
	return operator()(t0() + t * (t1() - t0()), u0() + u * (u1() - u0()));
}

SmoothParametricSurface SmoothParametricSurface::operator-(const SmoothParametricSurface& S) const {
	return *this + S * (-1);
}

SmoothParametricCurve operator&(const SmoothParametricSurface& S, const SmoothParametricPlaneCurve& c) {
	return S.precompose(c);
}

SmoothParametricSurface operator&(const SpaceEndomorphism& f, const SmoothParametricSurface& S) {
	return S.postcompose(f);
}

SmoothParametricSurface operator&(const SmoothParametricSurface& S, const PlaneAutomorphism& c) {
	return S.precompose(c);
}



float SmoothParametricSurface::t0() const {
	return rangeT.x;
}

float SmoothParametricSurface::t1() const {
	return rangeT.y;
}

float SmoothParametricSurface::u0() const {
	return rangeU.x;
}

float SmoothParametricSurface::u1() const {
	return rangeU.y;
}


float SmoothParametricSurface::periodT() const {
	return periodicT ? t1() - t0() : 0;
}

float SmoothParametricSurface::periodU() const {
	return periodicU ? u1() - u0() : 0;
}

vec3 SmoothParametricSurface::normal(vec2 tu) const {
	return normal(tu.x, tu.y);
}

mat2x3 SmoothParametricSurface::tangentStandardBasis(float t, float s) const {
	return mat2x3(_df_t(t, s), _df_u(t, s));
}

float SmoothParametricSurface::meanCurvature(vec2 tu) const {
	return meanCurvature(tu.x, tu.y);
}

SmoothParametricSurface ruledSurfaceJoinT(const SmoothParametricCurve& c1, const SmoothParametricCurve& c2, vec2 bounds) {
	return ruledSurfaceJoinT(c1, c2, bounds.x, bounds.y);
}

SmoothParametricSurface ruledSurfaceJoinU(const SmoothParametricCurve& c1, const SmoothParametricCurve& c2, vec2 bounds) {
	return ruledSurfaceJoinU(c1, c2, bounds.x, bounds.y);
}

SmoothParametricSurface CoonsPatchDisjoint(const SmoothParametricCurve& c1, const SmoothParametricCurve& c2) {
	SmoothParametricCurve c_bd = SmoothParametricCurve::span(c1(c1.bounds().x), c2(c2.bounds().x));
	return CoonsPatch(c1, c_bd, c2, c_bd);
}

SurfaceParametricPencil::SurfaceParametricPencil(const std::function<SmoothParametricSurface(float)>& pencil): pencil(pencil) {}

SmoothParametricSurface SurfaceParametricPencil::operator()(float t) const {
	return pencil(t);
}

vec3 SurfaceParametricPencil::operator()(float t, float u, float s) const {
	return pencil(t)(u, s);
}

vec3 SurfaceParametricPencil::operator()(float t, vec2 us) const {
	return pencil(t)(us);
}

CurveParametricPencil::CurveParametricPencil(const std::function<SmoothParametricCurve(float)>& pencil): pencil(pencil) {}

SmoothParametricCurve CurveParametricPencil::operator()(float t) const {
	return pencil(t);
}

vec3 CurveParametricPencil::operator()(float t, float u) const {
	return pencil(t)(u);
}

SmoothParametricCurve ParametricSurfaceFoliation::getLeaf(float t) const {
	return pencil_of_leaves(t);
}

SmoothParametricCurve ParametricSurfaceFoliation::getSpecialLeaf(int i) const {
	return special_leaves[i];
}

vector<SmoothParametricCurve> ParametricSurfaceFoliation::getSpecialLeaves() const {
	return special_leaves;
}

vec2 ParametricSurfaceFoliation::getDomain() const {
	return pencil_domain;
}

RealFunctionPS::RealFunctionPS(const std::function<float(float, float)>& f, const shared_ptr<SmoothParametricSurface>& surface): _f(f), surface(surface) {}

RealFunctionPS::RealFunctionPS(const std::function<float(vec2)>& f, const shared_ptr<SmoothParametricSurface>& surface): surface(surface), _f(pack(f, f, vec2, float)), _df(pack(f, derivativeOperator(f, .01), vec2, float)) {}

RealFunctionPS RealFunctionPS::constant(float c, const shared_ptr<SmoothParametricSurface>& surface) {
	return RealFunctionPS([c](float, float) {
		return c;
	}, surface);
}

RealFunctionPS RealFunctionPS::constant(float c) const {
	return constant(c, surface);
}

RealFunctionPS RealFunctionPS::operator*(float a) const {
	return RealFunctionPS([f=_f, a](float t, float s) {
		return f(t, s) * a;
	}, surface);
}

RealFunctionPS RealFunctionPS::operator/(float a) const {
	return (*this) * (1 / a);
}

RealFunctionPS RealFunctionPS::operator-() const {
	return (*this) * (-1);
}

RealFunctionPS RealFunctionPS::operator+(const RealFunctionPS& f) const {
	return RealFunctionPS([f1=_f, f2=f._f](float t, float s) {
		return f1(t, s) + f2(t, s);
	}, surface);
}

RealFunctionPS RealFunctionPS::operator-(const RealFunctionPS& f) const {
	return *this + (-f);
}

RealFunctionPS RealFunctionPS::operator*(const RealFunctionPS& f) const {
	return RealFunctionPS([f1=_f, f2=f._f](float t, float s) {
		return f1(t, s) * f2(t, s);
	}, surface);
}

RealFunctionPS RealFunctionPS::operator/(const RealFunctionPS& f) const {
	return RealFunctionPS([f1=_f, f2=f._f](float t, float s) {
		return f1(t, s) / f2(t, s);
	}, surface);
}

RealFunctionPS operator/(float a, const RealFunctionPS& g) {
	return g.constant(a) / g;
}

RealFunctionPS operator+(float a, const RealFunctionPS& g) {
	return g.constant(a) + g;
}

RealFunctionPS operator-(float a, const RealFunctionPS& g) {
	return g.constant(a) - g;
}

RealFunctionPS operator*(float a, const RealFunctionPS& g) {
	return g.constant(a) * g;
}

template <typename V>
Linear1Form2D<V>::Linear1Form2D(vec2 omega, V basis1, V basis2): v1(basis1), v2(basis2), coefs(omega) {}

template <typename V>
float Linear1Form2D<V>::operator()(V v) const {
	return dot(vec2(glm::dot(v, v1), glm::dot(v, v2)), coefs);
}

template <typename V>
Linear1Form2D<V> Linear1Form2D<V>::operator*(float a) const {
	return Linear1Form2D(coefs * a, v1, v2);
}

template <typename V>
Linear1Form2D<V> Linear1Form2D<V>::operator/(float a) const {
	return Linear1Form2D(coefs / a, v1, v2);
}

template <typename V>
Linear1Form2D<V> Linear1Form2D<V>::operator+(const Linear1Form2D& other) const {
	return Linear1Form2D(coefs + other.coefs, v1, v2);
}

template <typename V>
Linear1Form2D<V> Linear1Form2D<V>::operator-(const Linear1Form2D& other) const {
	return Linear1Form2D(coefs - other.coefs, v1, v2);
}

template <typename V>
Linear1Form2D<V> Linear1Form2D<V>::dual(V v, V v1, V v2) {
	return Linear1Form2D(vec2(glm::dot(v, v1), glm::dot(v, v2)), v1, v2);
}

template <typename V>
Linear1Form2D<V> Linear1Form2D<V>::dx(V x, V y) {
	return Linear1Form2D(vec2(1, 0), x, y);
}

template <typename V>
Linear1Form2D<V> Linear1Form2D<V>::dy(V x, V y) {
	return Linear1Form2D(vec2(0, 1), x, y);
}

template <typename V>
pair<Linear1Form2D<V>, Linear1Form2D<V>> Linear1Form2D<V>::basisForms(V v1, V v2) {
	return {dx(v1, v2), dy(v1, v2)};
}

template <typename V>
vec2 Linear1Form2D<V>::localCoefs() const {
	return coefs;
}

template <typename V>
Linear2Form2D<V>::Linear2Form2D(float c, V basis1, V basis2): v1(basis1), v2(basis2), coef(c) {}

template <typename V>
float Linear2Form2D<V>::operator()(V a, V b) const {
	return coef * glm::dot(a, v1) * glm::dot(b, v2);
}

template <typename V>
Linear2Form2D<V> Linear2Form2D<V>::operator*(float a) const {
	return Linear2Form2D(coef * a, v1, v2);
}

template <typename V>
Linear2Form2D<V> Linear2Form2D<V>::operator/(float a) const {
	return Linear2Form2D(coef / a, v1, v2);
}

template <typename V>
Linear2Form2D<V> Linear2Form2D<V>::operator+(const Linear2Form2D& other) const {
	return Linear2Form2D(coef + other.coef, v1, v2);
}

template <typename V>
Linear2Form2D<V> Linear2Form2D<V>::operator-(const Linear2Form2D& other) const {
	return Linear2Form2D(coef - other.coef, v1, v2);
}

template <typename V>
float Linear2Form2D<V>::localCoef() const {
	return coef;
}

Differential1FormPS::Differential1FormPS(const std::function<Linear1Form2D<vec3>(float, float)>& omega, const shared_ptr<SmoothParametricSurface>& surface): _omega(omega), surface(surface) {}

Differential1FormPS::Differential1FormPS(const std::function<Linear1Form2D<vec3>(vec2)>& omega, const shared_ptr<SmoothParametricSurface>& surface): _omega(pack(omega, omega, vec2, float)), surface(surface) {}

Linear1Form2D<vec3> Differential1FormPS::operator()(float t, float s) const {
	return _omega(t, s);
}

Linear1Form2D<vec3> Differential1FormPS::operator()(vec2 tu) const {
	return _omega(tu.x, tu.y);
}

float Differential1FormPS::operator()(float t, float s, vec3 v) const {
	return _omega(t, s)(v);
}

float Differential1FormPS::operator()(vec2 tu, vec3 v) const {
	return unpack(w=_omega, w, vec2)(tu)(v);
}

Differential1FormPS Differential1FormPS::operator*(float a) const {
	return Differential1FormPS([w=_omega, a](float t, float s) {
		return w(t, s) * a;
	}, surface);
}

Differential1FormPS Differential1FormPS::operator/(float a) const {
	return (*this) * (1 / a);
}

Differential1FormPS Differential1FormPS::operator-() const {
	return (*this) * (-1);
}

Differential1FormPS Differential1FormPS::operator+(const Differential1FormPS& eta) const {
	return Differential1FormPS([w1=_omega, w2=eta._omega](float t, float s) {
		return w1(t, s) + w2(t, s);
	}, surface);
}

Differential1FormPS Differential1FormPS::operator-(const Differential1FormPS& eta) const {
	return *this + (-eta);
}

Differential1FormPS Differential1FormPS::operator*(const RealFunctionPS& f) const {
	return Differential1FormPS([w=_omega, f](float t, float s) {
		return w(t, s) * f(t, s);
	}, surface);
}

Differential1FormPS Differential1FormPS::operator/(const RealFunctionPS& f) const {
	return Differential1FormPS([w=_omega, f](float t, float s) {
		return w(t, s) / f(t, s);
	}, surface);
}

Differential2FormPS::Differential2FormPS(const std::function<Linear2Form2D<vec3>(float, float)>& omega, const shared_ptr<SmoothParametricSurface>& surface): _omega(omega), surface(surface) {}

Differential2FormPS::Differential2FormPS(const std::function<Linear2Form2D<vec3>(vec2)>& omega, const shared_ptr<SmoothParametricSurface>& surface): _omega([omega](float t, float s) {
	return omega(vec2(t, s));
}), surface(surface) {}

float Differential2FormPS::operator()(float t, float s, vec3 v, vec3 w) const {
	return _omega(t, s)(v, w);
}

float Differential2FormPS::operator()(vec2 tu, vec3 v, vec3 w) const {
	return [om=_omega](vec2 tu, vec3 v, vec3 w) {
		return om(tu.x, tu.y)(v, w);
	}(tu, v, w);
}

Differential2FormPS Differential2FormPS::operator*(float a) const {
	return Differential2FormPS([w=_omega, a](float t, float s) {
		return w(t, s) * a;
	}, surface);
}

Differential2FormPS Differential2FormPS::operator/(float a) const {
	return (*this) * (1 / a);
}

Differential2FormPS Differential2FormPS::operator-() const {
	return (*this) * (-1);
}

Differential2FormPS Differential2FormPS::operator+(const Differential2FormPS& eta) const {
	return Differential2FormPS([w1=_omega, w2=eta._omega](float t, float s) {
		return w1(t, s) + w2(t, s);
	}, surface);
}

Differential2FormPS Differential2FormPS::operator-(const Differential2FormPS& eta) const {
	return *this + (-eta);
}

Differential2FormPS Differential2FormPS::operator*(const RealFunctionPS& f) const {
	return Differential2FormPS([w=_omega, f](float t, float s) {
		return w(t, s) * f(t, s);
	}, surface);
}

Differential2FormPS Differential2FormPS::operator/(const RealFunctionPS& f) const {
	return Differential2FormPS([w=_omega, f](float t, float s) {
		return w(t, s) / f(t, s);
	}, surface);
}

VectorFieldPS::VectorFieldPS(const std::function<vec3(float, float)>& f_dt, const std::function<vec3(float, float)>& f_du, const shared_ptr<SmoothParametricSurface>& surface): _f_dt(f_dt), _f_ds(f_du), surface(surface) {}

vec3 VectorFieldPS::operator()(float t, float s) const {
	return surface->tangentStandardBasis(t, s)[0] * _f_dt(t, s) + surface->tangentStandardBasis(t, s)[1] * _f_ds(t, s);
}

VectorFieldPS VectorFieldPS::operator*(float a) const {
	return VectorFieldPS([f=_f_dt, a](float t, float s) {
							 return f(t, s) * a;
						 }, [f=_f_ds, a](float t, float s) {
							 return f(t, s) * a;
						 }, surface);
}

VectorFieldPS VectorFieldPS::operator/(float a) const {
	return (*this) * (1 / a);
}

VectorFieldPS VectorFieldPS::operator+(const VectorFieldPS& v) const {
	return VectorFieldPS([f1=_f_dt, f2=v._f_dt](float t, float s) {
							 return f1(t, s) + f2(t, s);
						 }, [f1=_f_ds, f2=v._f_ds](float t, float s) {
							 return f1(t, s) + f2(t, s);
						 }, surface);
}

VectorFieldPS VectorFieldPS::operator-(const VectorFieldPS& v) const {
	return *this + (-v);
}

VectorFieldPS VectorFieldPS::operator-() const {
	return *this * (-1);
}

VectorFieldPS VectorFieldPS::operator*(const RealFunctionPS& f) const {
	return VectorFieldPS([f1=_f_dt, f](float t, float s) {
							 return f(t, s) * f1(t, s);
						 }, [f2=_f_ds, f](float t, float s) {
							 return f(t, s) * f2(t, s);
						 }, surface);
}

VectorFieldPS VectorFieldPS::operator/(const RealFunctionPS& f) const {
	return VectorFieldPS([f1=_f_dt, f](float t, float s) {
							 return f1(t, s) / f(t, s);
						 }, [f2=_f_ds, f](float t, float s) {
							 return f2(t, s) / f(t, s);
						 }, surface);
}

vec2 VectorFieldPS::deform(vec2 tu, vec2 dv) const {
	return tu + dv.x * vec2(_f_dt(tu.x, tu.y)) + dv.y * vec2(_f_ds(tu.x, tu.y));
}

vec3 VectorFieldPS::shiftAmbient(vec2 tu, vec2 dv) const {
	return (*surface)(deform(tu, dv)) - (*surface)(tu);
}

FunctionalPartitionOfUnity::FunctionalPartitionOfUnity(const std::vector<std::function<float(float)>>& F_i): _F_i(F_i) {}

std::function<float(float)> FunctionalPartitionOfUnity::operator[](int i) const {
	return _F_i[i];
}

int FunctionalPartitionOfUnity::size() const {
	return _F_i.size();
}

std::function<float(float)> BernsteinPolynomial(int n, int i, float t0, float t1) {
	return [n, i, t0, t1](float t) {
		return binomial(n, i) * std::pow(t - t0, i) * std::pow(t1 - t, n - i) / std::pow(t1 - t0, n);
	};
}

std::function<float(float)> BernsteinPolynomial(int n, int i) {
	return [n, i](float t) {
		return binomial(n, i) * std::pow(t, i) * std::pow(1 - t, n - i);
	};
}

SmoothParametricCurve BezierCurve(const std::vector<vec3>& controlPoints, float t0, float t1, float eps) {
	return freeFormCurve(BernsteinBasis(controlPoints.size() - 1, t0, t1), controlPoints, vec2(t0, t1), eps);
}

SmoothParametricCurve BSplineCurve(const std::vector<vec3>& controlPoints, float t0(), float t1(), int k, float eps) {
	return BSplineCurve(controlPoints, uniformKnots(controlPoints.size() - 1, k), k, eps);
}

SmoothParametricSurface BezierSurface(const std::vector<std::vector<vec3>>& controlPoints, float t0, float t1, float u0, float u1, float eps) {
	return freeFormSurface(BernsteinBasis(controlPoints.size() - 1, t0, t1), BernsteinBasis(controlPoints[0].size() - 1, u0, u1), controlPoints, vec2(t0, t1), vec2(u0, u1), eps);
}

SmoothParametricSurface BSplineSurfaceUniform(const std::vector<std::vector<vec3>>& controlPoints, vec2 t_range, vec2 u_range, int k, float eps) {
	return BSplineSurface(controlPoints, linspace(t_range.x, t_range.y, controlPoints.size() + k + 1), linspace(u_range.x, u_range.y, controlPoints[0].size() + k + 1), k, eps);
}

SmoothParametricCurve::SmoothParametricCurve(const RealFunction& fx, const RealFunction& fy, const RealFunction& fz, std::variant<int, std::string> id, float t0, float t1, bool periodic, float epsilon)
	: SmoothParametricCurve([fx, fy, fz](float t) {return vec3(fx(t), fy(t), fz(t)); },
	[fx, fy, fz](float t) {return vec3(fx.df(t), fy.df(t), fz.df(t)); },
	[fx, fy, fz](float t) {return vec3(fx.ddf(t), fy.ddf(t), fz.ddf(t)); }, std::move(id), t0, t1, periodic, epsilon) {}


vec3 SmoothParametricCurve::tangent(float t) const {
	vec3 dxdt = _df(t);
	if (norm2(dxdt) < eps*eps) {
		dxdt = _f(t + 2*eps) - _f(t - 2*eps);
		if (norm(dxdt) < eps*eps) {
			dxdt = _f(t + 5*eps) - _f(t - 3*eps);
			if (norm(dxdt) < eps*eps)
				return e1;
		}
	}
		return normalise(dxdt);
}

vec3 SmoothParametricCurve::binormal(float t) const {
	return normalise(cross(tangent(t), normal(t)));
}

vec3 SmoothParametricCurve::normal(float t) const {
	vec3 v = tangent(t+eps) - tangent(t - eps);
	if (norm(v) < eps*eps) {
		vec3 vv = tangent(t+2*eps) - tangent(t-2*eps);

		if (norm(vv) < eps*eps) {
			vec3 T = tangent(t);
			if (norm(T) < eps)
				T = tangent(t+eps*4);
			if (norm(T) < eps)
				T = tangent(t-eps*4);
			if (norm(T) < eps)
				return e3;

			vec3 N = projectVectorToPlane(e3, T);
			if (norm(N) < .1)
				N = projectVectorToPlane(e2, T);
			if (norm(N) < .1)
				N = projectVectorToPlane(e1, T);

			return normalise(N);
		}
		return normalise(vv);
	}
	return normalise(v);
}

SmoothParametricSurface SmoothParametricCurve::surfaceOfRevolution(const AffineLine& axis) const
{
	Foo113 param = [f = _f, axis](float t, float s) {
		vec3 p0 = axis.orthogonalProjection(f(t));
		vec3 v1 = f(t) - p0;
		vec3 v2 = cross(axis.direction(), v1);
		return p0 + v1*cos(s) + v2*sin(s);
	};
	return SmoothParametricSurface(param, vec2(t0, t1), vec2(0, TAU), periodic, true, eps);
}

SmoothParametricSurface SmoothParametricCurve::screwMotion(float speed, int iterations) const {
	Foo113 param = [f = _f, speed](float t, float s) {
		return vec3(f(t).x*cos(s) - f(t).y*sin(s), f(t).x*sin(s) + f(t).y*cos(s), f(t).z + speed*s);
	};
	return SmoothParametricSurface(param, vec2(t0, t1), vec2(0, TAU*iterations), periodic, false, eps);
}

SmoothParametricSurface SmoothParametricCurve::cylinder(vec3 direction, float h) const {
	Foo113 param = [f = _f, direction](float t, float s) {
		return f(t) + normalise(direction)*s;
	};
	return SmoothParametricSurface(param, vec2(t0, t1), vec2(0, h), periodic, false, eps);
}

SmoothParametricSurface SmoothParametricCurve::pipe(float radius, float eps, bool useFrenetFrame) const {
	Foo113 param = [c=*this, radius](float t, float s) {
		vec3 n = c.normal(t);
		vec3 b = c.binormal(t);
		if (dot(n, e2) < 0) {
			n = -n;
			b = -b;
		}

		// if (dot(b, e3) < 0)
		// 	b = -b;
		return c(t) + (n*cos(s) + b*sin(s))*radius;

	};
	return SmoothParametricSurface(param, vec2(t0, t1), vec2(0, TAU), periodic, true, this->eps/3);
}
SmoothParametricSurface SmoothParametricCurve::pipe(float radius, bool useFrenetFrame) const {
	return canal([radius](float t) {return radius; });
}


SmoothParametricSurface SmoothParametricCurve::canal(const function<float(float)>& radius) const {
	Foo113 param = [c=*this, r=radius, dr=derivativeOperator(radius, eps)](float t, float s) {
		return c(t) - r(t)*dr(t)/norm2(c.df(t))*c.df(t) + (c.normal(t)*cos(s) + c.binormal(t)*sin(s))*r(t)*sqrt(1 - dr(t)*dr(t)/norm2(c.df(t)));
	};
	return SmoothParametricSurface(param, vec2(t0, t1), vec2(0, TAU), periodic, true, eps);
}

SmoothParametricCurve SmoothParametricCurve::span(vec3 p1, vec3 p2) { return SmoothParametricCurve([p1, p2](float t) {return p1*(1-t) + p2*t; },
																								   [p1, p2](float t) {return p2-p1; }, [p1, p2](float t) {return vec3(0); }, randomID(), 0, 1, false, .01); }

FunctionalPartitionOfUnity BernsteinBasis(int n) {
	vector<Fooo> basis = {};
	for (int i = 0; i <= n; i++)
		basis.push_back(BernsteinPolynomial(n, i));
	return FunctionalPartitionOfUnity(basis);
}
FunctionalPartitionOfUnity BernsteinBasis(int n, float t0, float t1) {
	vector<Fooo> basis = {};
	for (int i = 0; i <= n; i++)
		basis.push_back(BernsteinPolynomial(n, i, t0, t1));
	return FunctionalPartitionOfUnity(basis);
}

function<float(float)> BSpline(int i, int k, const std::vector<float> &knots) {
	return [i, k, &knots](float t) {
		if (t <= knots[i] || t >= knots[i + k + 1])
			return 0.f;
		if (k==1)
			return knots[i] <= t && t < knots[i+1] ? 1.f : 0.f;
		return (t-knots[i])/(knots[i+k-1]-knots[i])*BSpline(i, k-1, knots)(t) + (knots[i+k]-t)/(knots[i+k]-knots[i+1])*BSpline(i+1, k-1, knots)(t);
	};
}


FunctionalPartitionOfUnity BSplineBasis(int n, int k, const std::vector<float>& knots) {
	vector<Fooo> basis = {};
	for (int i = 0; i <= n-1; i++)
		basis.push_back(BSpline(i, k, knots));
	return FunctionalPartitionOfUnity(basis);
}

SmoothParametricCurve freeFormCurve(const FunctionalPartitionOfUnity& family, const std::vector<vec3>& controlPts,
	vec2 domain, float eps) {
	return SmoothParametricCurve([controlPts, family](float t) {
		vec3 res = vec3(0);
		for (int i = 0; i < family.size(); i++)
			res += family[i](t)*controlPts[i];
		return res;
	}, randomID(), domain.x, domain.y, false, eps);
}

int findKnot(const std::vector<float> &knots, float t) {
	for (int i = 0; i < knots.size()-1; i++)
		if (t >= knots[i] && t < knots[i+1])
			return i;
	return -1;
}

vector<float> uniformKnots(int n, int k) {
	vector<float> knots = {};
	knots.reserve(n+k+1);
	for (int i = 0; i < n+k+1; i++)
		knots.push_back(i);
	return knots;
}

SmoothParametricCurve BSplineCurve(const std::vector<vec3> &controlPoints, const std::vector<float> &knots, int k, float eps) {
	auto basis = BSplineBasis(controlPoints.size()-1, k, knots);
	return SmoothParametricCurve([controlPoints, basis, knots, k](float t) {
		vec3 res = vec3(0);
		int l = findKnot(knots, t);
		for (int i = max(0, l-k+1); i <=l; i++)
			res += controlPoints[i]*basis[i](t);
		return res;
	}, 0, knots[0], knots[knots.size()-1], controlPoints[0] == controlPoints[controlPoints.size()-1], eps);
}


SmoothParametricSurface BSplineSurface(const std::vector<std::vector<vec3>> &controlPoints, const std::vector<float> &knots_t, const std::vector<float> &knots_u, int k, float eps) {
	auto basis_t = BSplineBasis(controlPoints.size()-1, k, knots_t);
	auto basis_u = BSplineBasis(controlPoints[0].size()-1, k, knots_u);
	return SmoothParametricSurface([controlPoints, basis_t, basis_u, knots_t, knots_u, k](float t, float u) {
		vec3 res = vec3(0);
		int l_t = findKnot(knots_t, t);
		int l_u = findKnot(knots_u, u);
		for (int i = max(0, l_t-k+1); i <=l_t; i++)
			for (int j = max(0, l_u-k+1); j <=l_u; j++)
				res += basis_t[i](t)*basis_u[j](u)*controlPoints[i][j];
		return res;
	},  vec2(knots_t[0], knots_t[knots_t.size()-1]), vec2(knots_u[0], knots_u[knots_u.size()-1]), false, false, eps);
}


SmoothParametricSurface freeFormSurface(const FunctionalPartitionOfUnity &F_i,
										const FunctionalPartitionOfUnity &G_i,
										const std::vector<std::vector<vec3>> &controlPts,
										vec2 range_t, vec2 range_u, float eps) {
	return SmoothParametricSurface([F_i, G_i, controlPts](float t, float s) {
		vec3 res = vec3(0);
		for (int i = 0; i < F_i.size(); i++)
			for (int j = 0; j < G_i.size(); j++)
				res += controlPts[i][j]* F_i[i](t)*G_i[j](s);
		return res;
	}, range_t, range_u, false, false, eps);
}


SmoothParametricSurface ruledSurfaceJoinT(const SmoothParametricCurve &c1, const SmoothParametricCurve &c2, float u0, float u1) {
	return SmoothParametricSurface([c1, c2, u0, u1](float t, float u) {
		float param = (u-u0)/(u1-u0);
		return (1-param)*c1(t) + param*c2(t);
	},
	c1.bounds(), vec2(u0, u1), c1.isPeriodic(), false, c1.getEps());
}


SmoothParametricSurface ruledSurfaceJoinU(const SmoothParametricCurve &c1, const SmoothParametricCurve &c2, float t0, float t1) {
	return SmoothParametricSurface([c1, c2, t0, t1](float t, float u) {
		float param = (t-t0)/(t1-t0);
		return (1-param)*c1(u) + param*c2(u);
	},
	vec2(t0, t1), c1.bounds(), false, c1.isPeriodic(), c1.getEps());
}

SmoothParametricSurface bilinearSurface(const vec3 &p00, const vec3 &p01, const vec3 &p10, const vec3 &p11, vec2 t_range, vec2 u_range) {
	return SmoothParametricSurface(
		[p00, p01, p10, p11, t_range, u_range](float t_, float u_) {
		float t = (t_-t_range.x)/(t_range.y-t_range.x);
		float u = (u_-u_range.x)/(u_range.y-u_range.x);
		return (1-t)*(1-u)*p00 + (1-t)*u*p01 + t*(1-u)*p10 + t*u*p11; },
		[p00, p01, p10, p11, t_range, u_range](float t_, float u_) {
				float t = 1.f/(t_range.y-t_range.x);
				float u = (u_-u_range.x)/(u_range.y-u_range.x);
				return -t*(1-u)*p00 + -t*u*p01 + t*(1-u)*p10 + t*u*p11; },
				[p00, p01, p10, p11, t_range, u_range](float t_, float u_) {
		float t = (t_-t_range.x)/(t_range.y-t_range.x);
		float u = 1.f/(u_range.y-u_range.x);
		return -(1-t)*u*p00 + (1-t)*u*p01 - t*u*p10 + t*u*p11; },
		t_range, u_range, false, false);
}

SmoothParametricSurface CoonsPatch(const SmoothParametricCurve &cDown, const SmoothParametricCurve &cLeft, const SmoothParametricCurve &cUp, const SmoothParametricCurve &cRight) {
	return ruledSurfaceJoinU(cDown, cUp, cLeft.bounds()) + ruledSurfaceJoinT(cLeft, cRight, cDown.bounds()) - bilinearSurface(cDown(cDown.bounds().x), cDown(cDown.bounds().y), cUp(cUp.bounds().x), cUp(cUp.bounds().y), cLeft.bounds(), cDown.bounds());
}

SmoothParametricSurface polarCone(const SmoothParametricCurve &r, vec3 center) {
	return SmoothParametricSurface([r, center](float t, float s) { return (r(t)-center)*s+center; }, r.bounds(), vec2(0, 1), true, false, r.getEps());
}

SurfaceParametricPencil::SurfaceParametricPencil(const BIHOM(float, vec2, vec3) &foo, vec2 range_t, vec2 range_u, float eps)
: pencil([foo, eps, range_t, range_u](float t) {
	return SmoothParametricSurface([foo=foo, t](float u, float s) {
		return foo(t, vec2(u, s));
	}, range_t, range_u, false, false, eps);
}) {}

CurveParametricPencil::CurveParametricPencil(std::function<vec3(float, float)> foo, vec2 bounds, float eps)
: pencil([foo, eps, bounds](float t) {
	return SmoothParametricCurve([foo, t](float u ) {
		return foo(t, u);
	}, randomID(), bounds.x, bounds.y, false, eps);
}) {}

ParametricSurfaceFoliation::ParametricSurfaceFoliation( std::function<SmoothParametricCurve(float)> pencil, vec2 pencil_domain, bool periodic, const vector<SmoothParametricCurve> &special_leaves)
: pencil_of_leaves(pencil), pencil_domain(pencil_domain), pencil_periodic(periodic), special_leaves(special_leaves) {}

SmoothParametricSurface ParametricSurfaceFoliation::getFoliatedSurface() const {
	auto leaf = pencil_of_leaves(pencil_domain.x);
	return SmoothParametricSurface([p=pencil_of_leaves](float t, float u) { return p(t)(u); }, pencil_domain, leaf.bounds(), pencil_periodic, leaf.isPeriodic(), leaf.getEps());
}

vector<SmoothParametricCurve> ParametricSurfaceFoliation::sampleLeaves(int res) const {
	vector<SmoothParametricCurve> leaves;
	leaves.reserve(res);
	for (float t: linspace(pencil_domain.x, pencil_domain.y, res))
		leaves.emplace_back(pencil_of_leaves(t));
	return leaves;
}
