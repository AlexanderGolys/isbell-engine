#include "common/glsl_utils.hpp"
#include "common/specific.hpp"
#include <memory>
#include <cmath>
#include <stdlib.h>

using namespace glm;
using std::vector, std::string, std::map, std::shared_ptr, std::unique_ptr, std::pair, std::make_unique, std::make_shared;

float debugSin(float t)
{
	return -cos(t)/2 + 0.5;
}

float width1(float t)
{
    return .03f + .00f*sin(4*t)*sin(4*t);
}

float width2(float t)
{
    return .03f + .00f*sin(4*t)*sin(4*t);
}

int main(void) {
	Renderer renderer = Renderer(.15f, vec4(.016f, .01509f, 0.05585f, 1.0f));
	renderer.initMainWindow(UHD, "flows");
	float camSpeed = 5.f;

	PolyGroupID curveID = PolyGroupID(420);
	PolyGroupID axisId = PolyGroupID(69);

	SmoothParametricCurve curve = epitrochoid(.2, 1, .15).embedding(e1, e3, vec3(4, 0, 0));
	AffineLine axis = AffineLine(vec3(1, -1, 0), vec3(0, 0, 1));
	float speedScrew = .0;
	//    SmoothParametricSurface surface = curve.screwMotion(speedScrew, 4);
	SmoothParametricSurface axisSurface = axis.tube(.02, 0, 8);

	SmoothParametricSurface surface = DupinCyclide(1, .8, .8, .0001);



	auto camCurve = make_shared<SmoothParametricCurve>([camSpeed, speedScrew](float t) { return vec3(-5*cos(t*camSpeed), -5*sin(t*camSpeed),  -3); }, curveID, 0, TAU, true, .01);
	auto camCurve2 = make_shared<SmoothParametricCurve>([speedScrew, axis](float t) { return vec3(0); }, curveID, 0, TAU, true, .01);

	shared_ptr<Camera> camera = make_shared<Camera>(camCurve, camCurve2, vec3(0, 0, -1));
	auto lights = vector({std::make_shared<PointLight>(vec3(1, 4, 0), vec4(.90952795, .785579, .94197285861, 1), 9.0f),
						  std::make_shared<PointLight>(vec3(1, 0, -3.1), vec4(.898769, .75369864, .903, 1), 9.0f),
						  std::make_shared<PointLight>(vec3(-2, -.1, -4), vec4(.9698, .9292598, .9399785938, 1), 5.0f)});

	auto midmat =  MaterialPhong(BLUE_PALLETTE[5], 0.031031423, .5962453956641656 , 0.21931160145739731, 50.0, "");


	auto shader = make_shared<Shader>(
			"C:\\Users\\PC\\Desktop\\ogl-master\\src\\shaders\\rev.vert",
			"C:\\Users\\PC\\Desktop\\ogl-master\\src\\shaders\\rev.frag");

	auto step = make_shared<RenderingStep>(shader);
	//	auto step2 = make_shared<RenderingStep>(shader);

	PolyGroupID sphereID = PolyGroupID(421);


	auto s = make_shared<WeakSuperMesh>(surface, 200, 400, curveID);
//	auto F = RealFunctionR3([](vec3 p) {return pow(p.x, 4) + pow(p.y, 4) + pow(p.z, 4) - 9;});
//	auto ts = TriangulatedImplicitSurface(F);
//	auto s = make_shared<WeakSuperMesh>(ts.compute(10, vec3(1), .05));

//	auto axis_s = make_shared<WeakSuperMesh>(axisSurface, 15, 100, axisId);

//	s->addUniformSurface(axisSurface, 15, 100, axisId);

    s->addGlobalMaterial(midmat);
//	axis_s->addGlobalMaterial(axmat);
    //
    // oscCircMid->scaleR(oscCircBand->getR()*2, false);
    // oscCircMid->scaleR(oscCircBand->getR()/2, true);


    // supeer->addUniformSurface(curva1(0), 150, 10, curveID);

    step->setWeakSuperMesh(s);
	step->addConstVec4("bdColor", vec4(85.f/255, 193.f/255, 255.f/255, 1.0f));

//	step2->setWeakSuperMesh(axis_s);


    // step->addUniform("intencities", MAT4, spec_uni([camera](float t, const shared_ptr<Shader> &shader) {shader->setUniform("mvp", camera->mvp(t));}));
    renderer.addRenderingStep(step);
//	renderer.addRenderingStep(step2);
//	renderer.nonlinearSpeed([speedScrew](float t) { return .4;} );
	renderer.addTimeUniform();

	renderer.setCamera(camera);
	renderer.setLights(lights);
	// renderer.addConstUniform("intencities", VEC4

	return renderer.mainLoop();

}
