#pragma once

#include "mat.hpp"
#include "../geometry/smoothImplicit.hpp"
//#include "src/geometry/smoothImplicit.hpp"




template <typename V>
class ODESolver {
protected:
	BIHOM(float, V, V) _f;
	float _t0;
	V initial;
	std::vector<V> solutionSequence = {initial};
	std::vector<float> times = {_t0};

public:
	virtual ~ODESolver() = default;
	ODESolver( BIHOM(float, V, V) &f, float t0, const V &initial) : _f(f), _t0(t0), initial(initial) {}
	virtual float getStep() { throw std::format_error("not implemented"); }
	virtual void computeStep() { throw std::format_error("not implemented"); }
	virtual V operator[](int i) const { return solutionSequence[i]; }
	virtual float timeAtStep(int i) const { return times[i]; }
	virtual float timeReached() const { return times.back(); }
	virtual bool solutionValidAtTime(float t) const { return t >= _t0 && t <= times.back(); }
	virtual int lowerIndexOfTime(float t) const { return binSearch(times, t); }
	virtual std::vector<V> solution() const { return solutionSequence; }
};

//template <VectorSpaceConcept<float> V>
class RK4 : public ODESolver<vec3> {
	vec69 c = {0, 0, 1.f /3, 2.f / 3, 1};
	MATR$X a = {{0, 0, 0, 0},
		{0, 0, 0, 0},
		{0, 1.f/3, 0, 0},
		{0, -1.f/3, 1.f, 0},
		{0, 1, -1, 1}};
	vec69 b = {0, 1.f/8, 3.f/8, 3.f/8, 1.f/8};
	float h;
public:
	RK4( BIHOM(float, vec3, vec3) &f, float t0, const vec3 &initial, float h) : ODESolver<vec3>(f, t0, initial), h(h) {}
	float getStep() override { return h; }
	void computeStep() override;
	void solveUpTo(float t1);
	void solveNSteps(int n);
	SmoothParametricCurve integralCurveBezier(float t0, float t1);
	std::vector<vec3> solution(int n)  {solveNSteps(n); return solutionSequence;}
	std::vector<vec3> solution() const override {return solutionSequence;}
};






SmoothParametricSurface BezierSurfaceFromSolution(const vector<vector<vec3>> &solutions, vec2 spatial_range, float t0, float t1);
SmoothParametricSurface SplineSurfaceFromSolution(const vector<vector<vec3>> &solutions, vec2 spatial_range, float t0, float t1, int k=3, float eps=.01);
vector<SmoothParametricCurve> streamLinesFromSolution(const vector<vector<vec3>> &solutions, float t0, float t1, int k=3, float eps=.01, int step=1);




vector<vector<vec3>> developSolutionsAlongCurve( BIHOM(float, vec3, vec3) &f, float t0,
													   const vector<vec3> &initials, int iters, float h) {
	vector<vector<vec3>> solutions = {};

	for (vec3 initial : initials) {
		RK4 solver(f, t0, initial, h);
		solver.solveNSteps(iters);
		solutions.push_back(solver.solution());
	}
	return solutions;
}



vector<vector<vec3>> developSolutionsAlongCurve( BIHOM(float, vec3, vec3) &f, float t0, float t1,
													   const vector<vec3> &initials, int iters) {
	vector<vector<vec3>> solutions = {};
	float h = (t1 - t0) / iters;

	for (vec3 initial : initials) {
		RK4 solver(f, t0, initial, h);
		solver.solveUpTo(t1);
		solutions.push_back(solver.solution());
	}
	return solutions;
}

SmoothParametricSurface BezierSurfaceFromSolution(const vector<vector<vec3>> &solutions, float t0, float t1) {
	return BezierSurface(solutions, 0, 1, t0, t1);
}
SmoothParametricSurface SplineSurfaceFromSolution(const vector<vector<vec3>> &solutions, vec2 spatial_range, float t0, float t1, int k, float eps) {
	return BSplineSurfaceUniform(solutions, spatial_range, vec2(t0, t1), 3, eps);
}


vector<SmoothParametricCurve> streamLinesFromSolution(const vector<vector<vec3>> &solutions, float t0, float t1, int k, float eps, int step) {
	vector<SmoothParametricCurve> curves = {};
	for (const vector<vec3>& solution : rangeStep(solutions, step))
		curves.push_back(BezierCurve(solution, t0, t1, eps));
	return curves;
}







//template<VectorSpaceConcept<float> V>
void RK4::computeStep() {
	float tn = this->times.back();
	vec3 yn = this->solutionSequence.back();
	float h  = getStep();

	auto k1 = _f(tn, yn);
	auto k2 = _f(tn + c[2]*h, yn + k1* a[2][1]*h);
	auto k3 = _f(tn + c[3]*h, yn + k1*a[3][1]*h + k2 * a[3][2]*h);
	auto k4 = _f(tn + c[4]*h, yn + k1*a[4][1]*h + k2*a[4][2]*h + k3*a[4][3]*h);

	vec3 y = yn + (k1*b[1] + k2*b[2] + k3*b[3] + k4*b[4]) * h;
	this->solutionSequence.push_back(y);
	this->times.push_back(tn + h);
}

void RK4::solveUpTo(float t1) {
	while (timeReached() < t1)
		computeStep();
}

void RK4::solveNSteps(int n) {
	for (int i = 0; i < n; i++)
		computeStep();
}

SmoothParametricCurve RK4::integralCurveBezier(float t0, float t1) {
	while (!solutionValidAtTime(t1))
		computeStep();
	vector<vec3> points = {};
	for (int i = 0; i < solutionSequence.size(); i++)
		if (times[i] >= t0 && times[i] <= t1)
			points.push_back(solutionSequence[i]);
	return BezierCurve(points, t0, t1);
};
