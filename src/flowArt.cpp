#include "physics/rigid.hpp"

using namespace glm;
using std::vector, std::string, std::shared_ptr, std::unique_ptr, std::pair, std::make_unique, std::make_shared;

float w(float t, float t0, float speed, vec2 x, vec2 x0) {
	float total_dt = 3;
	float max_w = .004;
	return powerShroom(0, 0+total_dt, 2, 6)(t)*max_w;
}

vec4 color(float t, float t0, float speed, vec2 x, vec2 x0) {
//	return lerp(vec4(.8, .6, 0, 1), vec4(.9, .7, 0, 1), sin(PI*t));
	return lerp(vec4(.9, .5, 0, 1), vec4(.8, .3, 0, 1), sin(TAU*t+speed*5));

}

int main() {
	Renderer renderer = Renderer(.05f, vec4(.089, .089, 0.089, 1.0f));
	renderer.initMainWindow(UHD, "flows");


	float camSpeed = .00;


	PointLight light1 = PointLight(vec3(-1,3, 1), .01, .032);
	PointLight light2 = PointLight(vec3(3, 26, 20), vec4(1, 1, .8, 1), .009, .00132);
	PointLight light3 = PointLight(vec3(1,10, 2), .001, .0302);
	auto lights = vector<Light>({light1, light2, light3});
	auto bdmat = MaterialPhong(BLUE_PALLETTE[5], .07, .6, .48, 15.0);
	auto floormat = MaterialPhong(GRAY_PALLETTE[7], .13, .1, .01, 155.0);
	auto graymat = MaterialPhong(GRAY_PALLETTE[4], .13, .1, .01, 155.0);


	auto X = VectorFieldR2([](vec2 x) { return vec2(1.2*x.y+cos(3*x.y*x.x-x.y), cos(2*x.y*x.x)+x.x*(1.1+sin(3*x.x + 4*x.y*x.y))); });
	auto lines = PlanarFlowLines(X, .03, 300, w, color);
	lines.generateGrid(vec2(-2, -1), vec2(2, 1), ivec2(60, 50));
	lines.generateStartTimesAll0();
	lines.generateLines();



	auto shader = Shader( R"(C:\Users\PC\Desktop\ogl-master\src\shaders2\flowart.vert)",
		R"(C:\Users\PC\Desktop\ogl-master\src\shaders2\flowart.frag)");



	renderer.setLights(lights);
	renderer.addMeshStep(shader, make_shared<WeakSuperMesh>(lines), bdmat);

	shared_ptr<Camera> camera = make_shared<Camera>(vec3(0, .1, 4), vec3(0), vec3(0, 0, 1), PI/8);
	renderer.setCamera(camera);
	renderer.addTimeUniform();


	return renderer.mainLoop();

}
