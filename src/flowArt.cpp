#include "shading/glsl_utils.hpp"

using namespace glm;
using std::vector, std::string, std::shared_ptr, std::unique_ptr, std::pair, std::make_unique, std::make_shared;

float w(float t, float t0, float speed, vec2 x, vec2 x0) {
	float total_dt = 1.3;
	float max_w = .007;
	return powerShroom(0, 0+total_dt, 2, 15)(t)*max_w;
}

vec4 color(float t, float t0, float speed, vec2 x, vec2 x0) {
//	return lerp(vec4(.8, .6, 0, 1), vec4(.9, .7, 0, 1), sin(PI*t));
	return lerp(vec4(.9, .4, 0, 1), vec4(.8, .6, 0, 1), .5 + .5*sin(4*TAU*t+speed));

}


int main() {
	Renderer renderer = Renderer(.03f, vec4(.089, .089, 0.089, 1.0f));
	renderer.initMainWindow(UHD, "flows");


	float camSpeed = .00;


	PointLight light1 = PointLight(vec3(-1,3, 1), .01, .032);
	PointLight light2 = PointLight(vec3(3, 26, 20), vec4(1, 1, .8, 1), .009, .00132);
	PointLight light3 = PointLight(vec3(1,10, 2), .001, .0302);
	auto lights = vector<Light>({light1, light2, light3});
	auto bdmat = MaterialPhong(BLUE_PALLETTE[5], .07, .6, .48, 15.0);
	auto floormat = MaterialPhong(GRAY_PALLETTE[7], .13, .1, .01, 155.0);
	auto graymat = MaterialPhong(GRAY_PALLETTE[4], .13, .1, .01, 155.0);


	auto X = VectorFieldR2([](vec2 x) { return vec2(0.4*x.y*(x.y-0.2)-sin(7*x.y*x.x-11*x.y+1)-.2, cos(5*x.y*x.x)+x.x*(1.3+sin(13*x.x + 4*x.y*x.y))); });
	auto lines = PlanarFlowLines(X, .02, 400, w, color);
	lines.generateGrid(vec2(-1.5, -1), vec2(1.5, 1), ivec2(80, 70));
//	lines.generateRandomUniform(vec2(-1, -1), vec2(1, 1), 900);

//	lines.generateStartTimesUniform(.6);
	lines.generateStartTimesAll0();
	lines.generateLines();



	auto shader = ShaderProgram( R"(C:\Users\PC\Desktop\ogl-master\src\shaders2\flowart.vert)",
		R"(C:\Users\PC\Desktop\ogl-master\src\shaders2\flowart.frag)");



	renderer.setLights(lights);
	renderer.addMeshStep(shader, make_shared<WeakSuperMesh>(lines), bdmat);

	shared_ptr<Camera> camera = make_shared<Camera>(vec3(0, .1, 4), vec3(0), vec3(0, 0, 1), PI/8);
	renderer.setCamera(camera);
	renderer.addTimeUniform();


	return renderer.mainLoop();

}
