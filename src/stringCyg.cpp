#include "common/specific.hpp"
#include "fundamentals/macros.hpp"
#include "fundamentals/prob.hpp"
#include "geometry/pde_dicrete.hpp"
#include "common/interface.hpp"

using namespace glm;
using std::vector, std::string, std::shared_ptr, std::unique_ptr, std::pair, std::make_unique, std::make_shared;





int main() {
	Renderer renderer = Renderer(0.2, vec4(.0, .0, 0.0, 1.0f));
	renderer.initMainWindow(FHD, "flows");

	vec2 dom = vec2(0, 10);

	auto u0 = expImpulse_k(3, .4, .5)*.3;

	auto mat = MaterialPhong(REDPINK_PALLETTE[1], .037, .556, 0.248, 2.0);
	auto floormat = MaterialPhong(BLUE_PALLETTE[5], .02, .536, .2348, 2.0);

	float c = .1;
	float r = .02;

	// DiscreteRealFunctionR2 snt = DiscreteRealFunctionR2(
	// 	[](float t, float u) {
	// 		return sin(t*u);
	// 	}, vec2(-2, 2), vec2(0, 5), 100, 100);

	// auto fn = DiscreteRealFunction(max(2-abs(3*X_R-3), 0), vec2(0, 9), 400);
	auto fn = DiscreteRealFunction(1-pow(1-EXP_R(-pow(X_R-2, 2)*16), 16) + max(1.5-abs(4*(X_R-4)), 0), vec2(0, 10), 500);
	// auto cn = fn.convolve(kn)/100;
	float k = 0.0101;
	auto heat = HeatRealLineHomoDiscrete(fn, k, 15, 150).solution();
	 // fn = fn.fft().ifft().re();
	// auto mesh = make_shared<PipeCurveVertexShader>(heat[1], r, 21);
	// auto mesh = make_shared<PipeCurveVertexShader>(cn, r, 21);

	auto pip = make_shared<PipeCurveVertexShader>(fn, r, 21);
	auto mesh = make_shared<SurfacePlotDiscretisedMesh>(heat);

	// mesh->flipNormals();
	// mesh->recalculateNormals();
	// mesh->shift(vec3(-dom.y/2, 0, 0));



	PointLight light1 = PointLight(vec3(-1,3, 3), .03, .032);
	PointLight light2 = PointLight(vec3(1, -1, 6), .09, .0132);
	PointLight light3 = PointLight(vec3(-1,-2, -2), .01, .02);
	auto lights = vector<Light>({light1, light2, light3});

	WeakSuperMesh b0 = disk3d(7, vec3(0, 0, -.6), e1, e2, 7, 4, randomID());

	auto shader = ShaderProgram(
		R"(C:\Users\PC\Desktop\ogl-master\src\shaders2\curve.vert)",
		R"(C:\Users\PC\Desktop\ogl-master\src\shaders2\heat.frag)");
	auto shader_floor = ShaderProgram(
	R"(C:\Users\PC\Desktop\ogl-master\src\shaders2\guitar.vert)",
	R"(C:\Users\PC\Desktop\ogl-master\src\shaders2\heat.frag)");

	shared_ptr<Camera> camera = make_shared<Camera>(
		make_shared<SmoothParametricCurve>([](float t) { return vec3(3+5*sin(-t), +5*cos(-t), 4); }),
		vec3(3, 0, 0),  vec3(0, 0, 1),  PI/4);
	renderer.setCamera(camera);
	renderer.setLights(lights);
	// renderer.addMeshStep(shader_floor, make_shared<WeakSuperMesh>(b0), floormat);
	// renderer.addMeshStep(shader_floor, mesh, floormat);
	renderer.addMeshStep(shader, pip, floormat);


	int freq = 3;
	float speed_wave = 1;
	float r0 = 0.01;
	float r_amp = 0.01;

	// auto rad = [freq, speed_wave, r0, r_amp](float t) {
	// 	return [t, freq, speed_wave, r0, r_amp](float x) {
	// 		return r0+r_amp*sin2(freq*x + speed_wave*t);
	// 	};
	// };

	renderer.addCustomAction([&](float t) {
		pip->updateCurve(heat(t));
	});

	renderer.addTimeUniform();
	return renderer.mainLoop();
}
